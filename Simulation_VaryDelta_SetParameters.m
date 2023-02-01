%% This file runs and saves the data for the results in Figure 1
% Parameters D, beta are held constant and Delta is randomly chosen.

%% Machine and simulation parameters
N = 1000; % Number of simulations
NDirs = 64;
SNR = [5, 10, 20]; % Signal-to-Noise Ratio
Sigma = 1./(sqrt(NDirs)*SNR); % St Dev of noise
numDelta = [1,2,3,4,5,2.1]; % 2.1 represents the special case of two Deltas that are the same
delta = 8*1E-3; % pulse duration (s)
gyroRatio = 267.522E6; % gyromagnetic ratio in units of rad/s/T
G = [0, 31, 68, 105, 142, 179, 216, 253, 290]' * 1E-3; % Gradient strengths in units of T/m
qvec = gyroRatio*G*delta/1E3; % units of rad/mm
DeltaMin = delta;
DeltaMax = delta + 50*1E-3;
diffDelta = 0.1*1E-3;
DeltaSep = 30*1E-3;

%% Set parameters for simulated data
Dvec = [3e-4,5e-4];
betaVec = [0.75,0.85];

%% Set Range of Delta values
Delta = DeltaMin:diffDelta:DeltaMax;
diffTime = Delta-delta/3;

%% Parameter estimation limits and initial guesses
D_guess = 1e3; Dmin = 0; Dmax = 100e3; % diffusion parameter, wide range used here, so same limits used for both sub-diffusion and DKI fitting, parameters are defined in mm^2/s^beta, but fitting is done in micro m^2/s^beta
beta_guess = 0.8; betamin = 0; betamax = 1; % time fractional index
% K_guess = 0; Kmin = 0; Kmax = 3; % K limits for DKI fitting

options = optimset('display','off','TolFun',1e-4,'TolX',1e-6);

%% Compute Noiseless data for all Delta and Q values
[Qvec, DiffTime] = meshgrid(qvec, diffTime);
[m,n] = size(Qvec);
numSim = length(Dvec);
NLData = zeros(m,n,numSim);
for ii = 1:numSim
    NLData(:,:,ii) = reshape(ml(-Dvec(ii).*Qvec(:).^2.*DiffTime(:).^betaVec(ii),betaVec(ii)), [m,n]);
end

for ii = 1:length(SNR) % for each SNR
    sigma = Sigma(ii);
    for ij = 1:length(Dvec) % for each parameter set
        D = Dvec(ij)*1e6;
        beta = betaVec(ij);
        K = 6*(gamma(1+beta)).^2./ gamma(1+2*beta) - 3;

        Dfit = zeros(N, length(numDelta));
        betaFit = zeros(N, length(numDelta));

        for ik = 1:length(numDelta) % for each set of number of Deltas
            nDelta = numDelta(ik);
            parfor il = 1:N % for the number of simulations per Delta and parameter set. These simulations are independent and can be performed in parallel

                %% choose Delta values and extract corresponding data and add noise
                if mod(nDelta,1)==0
                    nDelta2 = nDelta;
                    idx = randi(m, [nDelta,1]);
                    Deltas = Delta(idx);

                    while min(diff(sort(Deltas))) < DeltaSep/(nDelta-1) % Enforce minimum separation of Delta Values
                        idx = randi(m, [nDelta,1]);
                        Deltas = Delta(idx);
                    end

                    data = (squeeze(NLData(idx, :, ij)) + sigma*randn(nDelta,n))';
                    diff_time = DiffTime(idx,:)';
                    
                elseif nDelta == 2.1 % Incorporate special case of 2 Deltas that are the same value, i.e. repeated scan with same machine parameters
                    nDelta2 = 2;
                    idx = randi(m, [1,1]);
                    Deltas = Delta([idx;idx]);

                    data = (squeeze(NLData([idx;idx], :, ij)) + sigma*randn(2,n))';
                    diff_time = DiffTime([idx;idx],:)';
                else
                    error('code only set up to handle two Delta special case')
                end
                data = data(:);
                diff_time = diff_time(:);

                % Estimate parameters
                SUB = @(params, q) ml(-params(1) .* q.^2 .* diff_time .^params(2), params(2));

                [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], repmat(qvec, nDelta2,1), data, [Dmin betamin], [Dmax betamax], options);
                
                % Save parameter values from fitting
                Dfit(il, ik) = params(1)/1e6;
                betaFit(il, ik) = params(2);

            end
            % Display algorithm progress
            percentage = floor(( (ik + length(numDelta)*(ij-1+length(Dvec)*(ii-1))))/(length(SNR)*length(Dvec)*length(numDelta))*100);
            disp(['Simulation ', num2str(percentage), '% complete, SNR = ', num2str(SNR(ii)), ', parameter set ', num2str(ij), ', nDelta = ', num2str(nDelta),])
        end
        % Compute K values from beta.
        Kfit = 6*(gamma(1+betaFit)).^2./gamma(1+2*betaFit) - 3;
        % Save fitting results
        save(['VaryDeltaSim_SNR_', num2str(SNR(ii)), '_D_', num2str(D*1e4),'_beta_',num2str(beta*100)], 'D', 'beta', 'K', 'Dfit', 'betaFit', 'Kfit', 'numDelta')
    end
end
