%% This file runs and saves the data used for Figure 3.
% 1000 pairs of D, beta are simulated and fitted for Delta = 19, Delta =
% 49, Delta = 19,49 and Delta = 19,34,49.

%% Machine and simulation parameters
N = 1000; % Number of simulations
NDirs = 64;
SNR = [5,10,20, inf]; % Signal-to-Noise Ratio
Sigma = 1./(sqrt(NDirs)*SNR); % St Dev of noise
numDelta = [1,2,3,4,5,2.1]; % 2.1 represents the special case of two Deltas that are the same
delta = 8*1E-3; % pulse duration (s)
gyroRatio = 267.522E6; % gyromagnetic ratio in units of rad/s/T
G = [0, 31, 68, 105, 142, 179, 216, 253, 290]' * 1E-3; % Gradient strengths in units of T/m
qvec = gyroRatio*G*delta/1E3 /1E3; % units of rad/micro m

%% Define Delta values for figure 3
DeltaA = 19*1E-3; % Units of seconds.
DeltaB = 49*1E-3; % Units of seconds.
DeltaC = 34*1E-3; % Units of seconds.

diffTimeA = DeltaA-delta/3; % Units of seconds
diffTimeB = DeltaB-delta/3; % Units of seconds
diffTimeC = DeltaC-delta/3; % Units of seconds

bvecA = (qvec*1E3).^2 * diffTimeA; % Units of s/mm^2
DKI_idxA = bvecA<2500;
bvecA = bvecA(DKI_idxA);
bvecB = (qvec*1E3).^2 * diffTimeB; % Units of s/mm^2
DKI_idxB = bvecB<2500;
bvecB = bvecB(DKI_idxB);

%% Parameter estimation limits and initial guesses
D_guess = 1e3; Dmin = 0; Dmax = 100e3; % diffusion parameter, wide range used here, so same limits used for both sub-diffusion and DKI fitting
beta_guess = 0.8; betamin = 0; betamax = 1; % time fractional index
K_guess = 0; Kmin = 0; Kmax = 3; % K limits for DKI fitting

options = optimset('display','off','TolFun',1e-4,'TolX',1e-6);

%% Randomise the values of beta and D and generate data
min_beta = 0.5;
betaVec = rand(N,1)*(1-min_beta)+min_beta; % 0.5 <= beta <= 1
Dvec = rand(N,1)*(1e-3-1e-4)+1e-4; % 1e-4 <= D <= 1e-3 (mm^2/s^beta), fitting done in (micro m^2/s^beta)
KVec = 6*(gamma(1+betaVec).^2)./ gamma(1+2*betaVec) - 3; % Kurtosis value computed from beta

dataA = zeros(length(qvec),N);
dataB = zeros(length(qvec),N);
dataC = zeros(length(qvec),N);

for ii = 1:N
    dataA(:,ii) = ml(-Dvec(ii)*1E6.*qvec.^2.*diffTimeA.^betaVec(ii), betaVec(ii));
    dataB(:,ii) = ml(-Dvec(ii)*1E6.*qvec.^2.*diffTimeB.^betaVec(ii), betaVec(ii));
    dataC(:,ii) = ml(-Dvec(ii)*1E6.*qvec.^2.*diffTimeC.^betaVec(ii), betaVec(ii));
end

for ii = 1:length(Sigma) % For every Sigma/SNR value
    sigma = Sigma(ii);

    %% Initialise results matrices
    betaFitVec = zeros(N,4);
    DFitVec = zeros(N,4);

    DFitVec_DKI19 = zeros(N,1);
    KFitVec_DKI19 = zeros(N,1);
    DFitVec_DKI49 = zeros(N,1);
    KFitVec_DKI49 = zeros(N,1);

    parfor ij = 1:N
        %% Fit subdiffusion model for Delta = 19, Delta = 49, Delta = 19,49 and Delta = 19,34,49
        data19 = dataA(:,ij) + sigma*randn(length(qvec),1);
        data49 = dataB(:,ij) + sigma*randn(length(qvec),1);

        data19_49 = [dataA(:,ij);dataB(:,ij)] + sigma*randn(length(qvec)*2,1);
        data19_34_49 = [dataA(:,ij);dataB(:,ij);dataC(:,ij)] + sigma*randn(length(qvec)*3,1);

        SUB19 = @(params, q) ml(-params(1) .* q.^2 .* diffTimeA.^params(2), params(2));
        SUB49 = @(params, q) ml(-params(1) .* q.^2 .* diffTimeB.^params(2), params(2));
        SUB19_49 = @(params, q) ml(-params(1) .* q.^2 .* [diffTimeA*ones(size(qvec)); diffTimeB*ones(size(qvec))].^params(2), params(2));
        SUB19_34_49 = @(params, q) ml(-params(1) .* q.^2 .* [diffTimeA*ones(size(qvec)); diffTimeB*ones(size(qvec)); diffTimeC*ones(size(qvec))].^params(2), params(2));
        
        Dvals = zeros(1,4);
        betaVals = zeros(1,4);
        [params,resnorm] = lsqcurvefit(SUB19, [D_guess beta_guess], qvec, data19, [Dmin betamin], [Dmax betamax], options);
        Dvals(1)= params(1)/1E6;
        betaVals(1)  = params(2);

        [params,resnorm] = lsqcurvefit(SUB49, [D_guess beta_guess], qvec, data49, [Dmin betamin], [Dmax betamax], options);
        Dvals(2)= params(1)/1E6;
        betaVals(2)  = params(2);

        [params,resnorm] = lsqcurvefit(SUB19_49, [D_guess beta_guess], [qvec;qvec], data19_49, [Dmin betamin], [Dmax betamax], options);
        Dvals(3)= params(1)/1E6;
        betaVals(3)  = params(2);

        [params,resnorm] = lsqcurvefit(SUB19_34_49, [D_guess beta_guess], [qvec;qvec;qvec], data19_34_49, [Dmin betamin], [Dmax betamax], options);
        Dvals(4)= params(1)/1E6;
        betaVals(4)  = params(2);

        DFitVec(ij,:) = Dvals;
        betaFitVec(ij,:) = betaVals;

        %% Fit DKI model for Delta = 19 and Delta = 49 for comparison
        DKI = @(params,b) exp(-b*params(1) + 1/6*b.^2.*params(1)^2*params(2));
        [params, resnorm] = lsqcurvefit(DKI, [D_guess/1e6 K_guess], bvecA, data19(DKI_idxA), [Dmin/1e6 Kmin], [Dmax/1e6 Kmax], options);
        DFitVec_DKI19(ij) = params(1);
        KFitVec_DKI19(ij) = params(2);

        [params, resnorm] = lsqcurvefit(DKI, [D_guess/1e6 K_guess], bvecB, data49(DKI_idxB), [Dmin/1e6 Kmin], [Dmax/1e6 Kmax], options);
        DFitVec_DKI49(ij) = params(1);
        KFitVec_DKI49(ij) = params(2);


    end
    
save(['BetaDfit_', num2str(SNR(ii)), '_', num2str(min_beta*10),'.mat'], 'betaVec', 'betaFitVec', 'Dvec', 'DFitVec', 'DFitVec_DKI19', 'KFitVec_DKI19', 'DFitVec_DKI49', 'KFitVec_DKI49')

end
