%% This file performs a simulation and saves the data used in figure 4 and table 1.
% This simulation computes R^2 values for all combinations of reduce b values 
% to determine the optimal b-values to use in the b-value reduction.

%% Machine and simulation parameters
N = 1000; % Number of simulations
NDirs = 64;
SNR = [5,10,20, inf]; % Signal-to-Noise Ratio
Sigma = 1./(sqrt(NDirs)*SNR); % St Dev of noise
numDelta = [1,2,3,4,5,2.1]; % 2.1 represents the special case of two Deltas that are the same
delta = 8*1E-3; % pulse duration (s)
gyroRatio = 267.522E6; % gyromagnetic ratio in units of rad/s/T
G = [31, 68, 105, 142, 179, 216, 253, 290]' * 1E-3; % Gradient strengths in units of T/m
qvec = gyroRatio*G*delta/1E3 /1E3; % units of rad/micro m
numQ = length(qvec);
numB = [3,4,5];

%%

%% Define Delta values for figure 3
DeltaA = 19*1E-3; % Units of seconds.
DeltaB = 49*1E-3; % Units of seconds.

diffTimeA = DeltaA-delta/3; % Units of seconds
diffTimeB = DeltaB-delta/3; % Units of seconds

Qvec = [0; qvec; qvec];
diffTime = [0;diffTimeA*ones(size(qvec)); diffTimeB*ones(size(qvec))];

%% Define all combinations of q values
qIdx = 2:17;
% choose3 = nchoosek(qIdx,2);
% type3 = zeros(size(choose3,1),1);
% type3(all(choose3>9,2)) = 3; % Classify type as number of Delta2 + 1
% type3(sum(choose3<=9,2)==1) = 2;
% type3(all(choose3<=9,2)) = 1;
% choose3 = [ones(size(choose3,1),1), choose3];
% choose4 = nchoosek(qIdx,3);
% type4 = zeros(size(choose4,1),1);
% type4(all(choose4>9,2)) = 4;
% type4(sum(choose4<=9,2)==1) = 3;
% type4(sum(choose4<=9,2)==2) = 2;
% type4(all(choose4<=9,2)) = 1;
% choose4 = [ones(size(choose4,1),1), choose4];
% choose5 = nchoosek(qIdx,4);
% type5 = zeros(size(choose5,1),1);
% type5(all(choose5>9,2)) = 5;
% type5(sum(choose5<=9,2)==1) = 4;
% type5(sum(choose5<=9,2)==2) = 3;
% type5(sum(choose5<=9,2)==3) = 2;
% type5(all(choose5<=9,2)) = 1;
% choose5 = [ones(size(choose5,1),1), choose5];


%% Parameter estimation limits and initial guesses
D_guess = 1e3; Dmin = 0; Dmax = 100e3; % diffusion parameter, wide range used here, so same limits used for both sub-diffusion and DKI fitting
beta_guess = 0.8; betamin = 0; betamax = 1; % time fractional index
K_guess = 0; Kmin = 0; Kmax = 3; % K limits for DKI fitting

options = optimset('display','off','TolFun',1e-4,'TolX',1e-6);

%% Randomise the values of beta and D and generate data
minBeta = 0.5;
maxD = 1e-3; %in units micro m^2 / second
minD = 1e-4;

Dvec = rand(N,1)*(maxD-minD)+minD;
betaVec = rand(N,1)*(1-minBeta) + minBeta;
KVec = 6*(gamma(1+betaVec)).^2 ./ gamma(1+2*betaVec) - 3;



for ii = 1:length(SNR)
    sigma = Sigma(ii);
    data = zeros(N,length(Qvec));
    for ij = 1:N
        D = Dvec(ij)*1e6;
        data(ij,:) = ml(-D.*Qvec.^2.*diffTime.^betaVec(ij),betaVec(ij)) + sigma*randn(size(Qvec));
    end

    for ij = 1:length(numB)
        combsIdx = nchoosek(2:(numQ*2+1), numB(ij)-1);
        R_2_D = zeros(size(combsIdx,1),1);
        R_2_beta = zeros(size(combsIdx,1),1);
        R_2_K = zeros(size(combsIdx,1),1);
        for ik = 1:size(combsIdx,1)
            idx = [1, combsIdx(ik,:)];
            data_ij = data(:,idx);
            Qvec_ij = Qvec(idx);
            diffTime_ij = diffTime(idx);
            DFitVec = zeros(N,1);
            betaFitVec = zeros(N,1);

            parfor il = 1:N
                data19_49 = data_ij(il,:)';
                SUB19_49 = @(params, q) ml(-params(1) .* q.^2 .* diffTime_ij.^params(2), params(2));
                [params,resnorm] = lsqcurvefit(SUB19_49, [D_guess beta_guess], Qvec_ij, data19_49, [Dmin betamin], [Dmax betamax], options);
                DFitVec(il) = params(1)/1e6;
                betaFitVec(il) = params(2);
            end

            KFitVec = 6*(gamma(1+betaFitVec)).^2 ./ gamma(1+2*betaFitVec) - 3;

            R_2_D(ik) = 1 - sum( (Dvec - DFitVec).^2)./ sum( (Dvec - mean(Dvec)).^2 );
            R_2_beta(ik) = 1 - sum( (betaVec - betaFitVec).^2)./ sum( (betaVec - mean(betaVec)).^2 );
            R_2_K(ik) = 1 - sum( (KVec - KFitVec).^2)./ sum( (KVec - mean(KVec)).^2 );
        end
        save(['b_red_data_', num2str(numB(ij)), '_', num2str(SNR(ii))],'combsIdx', 'Dvec', 'DFitVec', 'betaVec', 'betaFitVec', 'KVec', 'KFitVec', 'R_2_D', 'R_2_beta', 'R_2_K')
    end
end