%% This file runs and saves the results used in Figure 2.
% for each set of Delta (2500 sets), 1000 simulations of randomly chosen
% D,beta pairs are used to test the R^2 values of simulations.

%% Machine and simulation parameters
N = 1000; % Number of simulations
NDirs = 64;
SNR = [5,10,20]; % Signal-to-Noise Ratio
Sigma = 1./(sqrt(NDirs)*SNR); % St Dev of noise
numDelta = [1,2,3,4,5,2.1]; % 2.1 represents the special case of two Deltas that are the same
delta = 8*1E-3; % pulse duration (s)
gyroRatio = 267.522E6; % gyromagnetic ratio in units of rad/s/T
G = [0, 31, 68, 105, 142, 179, 216, 253, 290]' * 1E-3; % Gradient strengths in units of T/m
qvec = gyroRatio*G*delta/1E3 /1E3; % units of rad/micro m
Qvec = [qvec; qvec];

%% Define range of Delta values for figure 2
deltarange = (9:58)*1e-3;
diffRange = (1:50)*1e-3;
[Delta1, dDelta] = meshgrid(deltarange, diffRange);
Delta2 = Delta1 + dDelta;

diffTime1 = Delta1 - delta/3;
diffTime2 = Delta2 - delta/3;

%% Parameter estimation limits and initial guesses
D_guess = 1e3; Dmin = 0; Dmax = 100e3; % diffusion parameter, wide range used here, so same limits used for both sub-diffusion and DKI fitting
beta_guess = 0.8; betamin = 0; betamax = 1; % time fractional index
% K_guess = 0; Kmin = 0; Kmax = 3; % K limits for DKI fitting

options = optimset('display','off','TolFun',1e-4,'TolX',1e-6);



%% Randomise the values of beta and D
betaVec = rand(N,numel(Delta1))*0.5+0.5; % 0.5 <= beta <= 1
Dvec = rand(N,numel(Delta1))*(1e-3-1e-4)+1e-4; % 1e-4 <= D <= 1e-3 (mm^2/s^beta), fitting done in (micro m^2/s^beta)
KVec = 6*(gamma(1+betaVec).^2)./ gamma(1+2*betaVec) - 3; % Kurtosis value computed from beta


for ii = 1:length(Sigma) % For every Sigma/SNR value
    sigma = Sigma(ii);

    %% Initialise results matrices
    betaFitVec = zeros(N,numel(Delta1));
    DFitVec = zeros(N,numel(Delta1));
    ErrorVec = zeros(N,numel(Delta1));
    for ij = 1:numel(Delta1) % For every (Delta1, Delta2) pair

        DeltaA = Delta1(ij);
        DeltaB = Delta2(ij);
        diffTimeA = diffTime1(ij);
        diffTimeB = diffTime2(ij);

        parfor ik = 1:N % for the number of simulations to be performed with each pair
            % parameter values for simulation
            beta = betaVec(ik,ij);
            D = Dvec(ik,ij)*1e6;

            % generate simulated data
            dataA = ml(-D.*qvec.^2.*diffTimeA.^beta,beta) + sigma*randn(size(qvec));
            dataB = ml(-D.*qvec.^2.*diffTimeB.^beta,beta) + sigma*randn(size(qvec));
            data = [dataA; dataB];

            % Estimate parameters
            SUB = @(params, q) ml(-params(1) .* q.^2 .* [diffTimeA*ones(size(dataA)); diffTimeB*ones(size(dataB))].^params(2), params(2));
            [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], Qvec, data, [Dmin betamin], [Dmax betamax], options);

            % Store results of parameter fitting
            DFitVec(ik,ij)= params(1)/1e6;
            betaFitVec(ik,ij)  = params(2);
            ErrorVec(ik,ij) = resnorm/length(Qvec);
        end

        if mod(ij,100)==0
        percentage = floor(((ii-1)*numel(Delta1)+ij)/(length(Sigma)*numel(Delta1))*100);
        disp([num2str(percentage), ' % complete'])
        end
    end
    % Compute K values from beta.
    KFitVec = 6*(gamma(1+betaFitVec).^2)./ gamma(1+2*betaFitVec) - 3;

    % Compute R^2 values
    R_2_D = reshape(1 - sum((Dvec - DFitVec).^2)./((N-1)*var(DFitVec)), size(Delta1));
    R_2_beta = reshape(1 - sum((betaVec - betaFitVec).^2)./((N-1)*var(betaFitVec)), size(Delta1));
    R_2_K = reshape(1 - sum((KVec - KFitVec).^2)./((N-1)*var(KFitVec)), size(Delta1));
    
    % Save results
    save(['DeltaCombinations', num2str(SNR(ii)), '_', num2str(N),'.mat'], 'diffRange', 'deltarange', 'R_2_D', 'R_2_beta', 'R_2_K', 'Dvec', 'DFitVec', 'betaVec', 'betaFitVec', 'KVec', 'KFitVec')
end

