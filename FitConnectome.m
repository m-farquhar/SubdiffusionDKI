%% This file produces and saves the data for all of the connectome data fitting results.

%% Machine and simulation parameters
delta = 8*1E-3; % pulse duration (s)
gyroRatio = 267.522E6; % gyromagnetic ratio in units of rad/s/T
G = [31, 68, 105, 142, 179, 216, 253, 290]' * 1E-3; % Gradient strengths in units of T/m
qvec = gyroRatio*G*delta/1E3 /1E3; % units of rad/micro m
qvecAll = [0;qvec; qvec];
Delta19 = 19*1e-3;
Delta49 = 49*1e-3;

diff_time19 = Delta19 - delta/3;
diff_time49 = Delta49 - delta/3;


bvec19 = [0; (qvec(1:5)*1e3).^2*diff_time19];
qvec19 = [0; qvec];

bvec49 = [0; (qvec(1:3)*1e3).^2*diff_time49];
qvec49 = [0; qvec];

diff_t = [0; repelem(diff_time19,8)'; repelem(diff_time49,8)'];

%% Parameter estimation limits and initial guesses
D_guess = 1e3; Dmin = 0; Dmax = 100e3; % diffusion parameter, wide range used here, so same limits used for both sub-diffusion and DKI fitting
beta_guess = 0.8; betamin = 0; betamax = 1; % time fractional index
K_guess = 0; Kmin = 0; Kmax = 3; % K limits for DKI fitting
D_guess_DKI = 1e-3; Dmin_DKI = 0; Dmax_DKI = 100e-3; % D limits for DKI fitting

options = optimset('display','off','TolFun',1e-4,'TolX',1e-6);

%% Specify which direction reductions to run, 0 is using all directions, 
RED_Dir = [true, true,false]; % are you reducing directions
Dir_32 = [16, 4, 0];
Dir_64 = [16, 4, 0];


%% Initialise which indices correspond to which Delta
idx19 = [1, 2:9];
idx49 = [1, 10:17];

%% loading trace data
subjects = {'sub_001', 'sub_003', 'sub_004', 'sub_005', 'sub_006', 'sub_007', 'sub_001_rescan', 'sub_003_rescan', 'sub_004_rescan', 'sub_005_rescan', 'sub_006_rescan', 'sub_007_rescan'};
for sub = 1:length(subjects) % Loop over all subjects
    SUBJECT = subjects{sub};
    for diridx = 1:length(RED_Dir) % Loop over all reduced directions you wish to do
        PATH   = pwd;
        FOLDER = [PATH '\' SUBJECT '\' 'dwi_real\'];

        % Load brain mask to only fit data within the brain
        brainmask = double(niftiread([FOLDER SUBJECT '_dwi_real_brainmask.nii.gz']));

        RedDir = RED_Dir(diridx);
        numDir32 = Dir_32(diridx);
        numDir64 = Dir_64(diridx);

        % Define Indices for Reduced b-values
        if numDir32 == 16
            %SNR = 10
            idx5a = [1,3,7,11,15]; % 1820 R2 = 0.9123
            idx4a = [1,3,6,12]; % 560 R2 = 0.8485
            idx5b = [1,5,7,14,16]; % 326 R2 = 0.4486
            idx4b = [1,6,8,12]; % 260 R2 = 0.4493
        elseif numDir32 == 4
            % SNR = 5
            idx5a = [1,3,6,11,13]; % 1820 R2 = 0.6309
            idx4a = [1,3,11,14]; % 560 R2 = 0.4879
            idx5b = [1,4,5,10,12];% 1263 R2 = 0.3000
            idx4b = [1,3,4,14]; % 514 R2 = 0.2961
        elseif numDir32 == 0
            % SNR = 20
            idx5a = [1,3,8,11,14]; % 1820 R2 = 0.9618
            idx4a = [1,2,6,13]; % 560 R3 = 0.9002
            idx5b = [1,7,8,12,13]; % 116 R2 = 0.4995
            idx4b = [1,6,8,12]; % 147 R2 = 0.4979
        end

        % Load in data
        if ~RedDir
            S0 = double(niftiread([FOLDER, SUBJECT,'_dwi_real_b0_delta0_image.nii']));

            % dwi_real with Delta = 19 ms
            S1 = double(niftiread([FOLDER SUBJECT '_dwi_real_b50_delta19_image.nii']));
            S2 = double(niftiread([FOLDER SUBJECT '_dwi_real_b350_delta19_image.nii']));
            S3 = double(niftiread([FOLDER SUBJECT '_dwi_real_b800_delta19_image.nii']));
            S4 = double(niftiread([FOLDER SUBJECT '_dwi_real_b1500_delta19_image.nii']));
            S5 = double(niftiread([FOLDER SUBJECT '_dwi_real_b2400_delta19_image.nii']));
            S6 = double(niftiread([FOLDER SUBJECT '_dwi_real_b3450_delta19_image.nii']));
            S7 = double(niftiread([FOLDER SUBJECT '_dwi_real_b4750_delta19_image.nii']));
            S8 = double(niftiread([FOLDER SUBJECT '_dwi_real_b6000_delta19_image.nii']));

            % dwi_real with Delta = 49 ms
            S9  = double(niftiread([FOLDER SUBJECT '_dwi_real_b200_delta49_image.nii']));
            S10 = double(niftiread([FOLDER SUBJECT '_dwi_real_b950_delta49_image.nii']));
            S11 = double(niftiread([FOLDER SUBJECT '_dwi_real_b2300_delta49_image.nii']));
            S12 = double(niftiread([FOLDER SUBJECT '_dwi_real_b4250_delta49_image.nii']));
            S13 = double(niftiread([FOLDER SUBJECT '_dwi_real_b6750_delta49_image.nii']));
            S14 = double(niftiread([FOLDER SUBJECT '_dwi_real_b9850_delta49_image.nii']));
            S15 = double(niftiread([FOLDER SUBJECT '_dwi_real_b13500_delta49_image.nii']));
            S16 = double(niftiread([FOLDER SUBJECT '_dwi_real_b17800_delta49_image.nii']));
        else
            S0 = double(niftiread([FOLDER, SUBJECT,'_dwi_real_b0_delta0_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));

            % dwi_real with Delta = 19 ms
            S1 = double(niftiread([FOLDER SUBJECT '_dwi_real_b50_delta19_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S2 = double(niftiread([FOLDER SUBJECT '_dwi_real_b350_delta19_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S3 = double(niftiread([FOLDER SUBJECT '_dwi_real_b800_delta19_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S4 = double(niftiread([FOLDER SUBJECT '_dwi_real_b1500_delta19_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S5 = double(niftiread([FOLDER SUBJECT '_dwi_real_b2400_delta19_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S6 = double(niftiread([FOLDER SUBJECT '_dwi_real_b3450_delta19_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S7 = double(niftiread([FOLDER SUBJECT '_dwi_real_b4750_delta19_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S8 = double(niftiread([FOLDER SUBJECT '_dwi_real_b6000_delta19_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));

            % dwi_real with Delta = 49 ms
            S9 = double(niftiread([FOLDER SUBJECT '_dwi_real_b200_delta49_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S10 = double(niftiread([FOLDER SUBJECT '_dwi_real_b950_delta49_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S11 = double(niftiread([FOLDER SUBJECT '_dwi_real_b2300_delta49_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S12 = double(niftiread([FOLDER SUBJECT '_dwi_real_b4250_delta49_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S13 = double(niftiread([FOLDER SUBJECT '_dwi_real_b6750_delta49_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S14 = double(niftiread([FOLDER SUBJECT '_dwi_real_b9850_delta49_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S15 = double(niftiread([FOLDER SUBJECT '_dwi_real_b13500_delta49_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
            S16 = double(niftiread([FOLDER SUBJECT '_dwi_real_b17800_delta49_image_dirRed_P_D_MM_n1_', num2str(numDir32), '_n2_', num2str(numDir64),'.nii']));
        end


        %% Initialise result arrays

        [X,Y,L] = size(S0);

        D_SUB = nan(X,Y,L);
        beta_SUB = nan(X,Y,L);
        error_SUB = nan(X,Y,L);


        D_SUB5a = nan(X,Y,L);
        beta_SUB5a = nan(X,Y,L);
        error_SUB5a = nan(X,Y,L);

        D_SUB5b = nan(X,Y,L);
        beta_SUB5b = nan(X,Y,L);
        error_SUB5b = nan(X,Y,L);

        D_SUB4a = nan(X,Y,L);
        beta_SUB4a = nan(X,Y,L);
        error_SUB4a = nan(X,Y,L);

        D_SUB4b = nan(X,Y,L);
        beta_SUB4b = nan(X,Y,L);
        error_SUB4b = nan(X,Y,L);

        if ~RedDir % Only fit DKI and single Delta models for full dataset
            D_SUB19 = nan(X,Y,L);
            beta_SUB19 = nan(X,Y,L);
            error_SUB19 = nan(X,Y,L);

            D_SUB49 = nan(X,Y,L);
            beta_SUB49 = nan(X,Y,L);
            error_SUB49 = nan(X,Y,L);

            D_DKI19 = nan(X,Y,L);
            K_DKI19 = nan(X,Y,L);
            error_DKI19 = nan(X,Y,L);

            D_DKI49 = nan(X,Y,L);
            K_DKI49 = nan(X,Y,L);
            error_DKI49 = nan(X,Y,L);
        end

        
% Loop over every voxel within the brain
        tic
        for l = 1:L
            for x = 1:X
                parfor y = 1:Y
                    if brainmask(x,y,l) && ~isnan(S0(x,y,l)) && S0(x,y,l)~=0
                        Svec = zeros(17,1);
                        Svec(1)= S0(x,y,l);
                        Svec(2)= S1(x,y,l);
                        Svec(3)= S2(x,y,l);
                        Svec(4)= S3(x,y,l);
                        Svec(5)= S4(x,y,l);
                        Svec(6)= S5(x,y,l);
                        Svec(7)= S6(x,y,l);
                        Svec(8)= S7(x,y,l);
                        Svec(9)= S8(x,y,l);
                        Svec(10)= S9(x,y,l);
                        Svec(11)= S10(x,y,l);
                        Svec(12)= S11(x,y,l);
                        Svec(13)= S12(x,y,l);
                        Svec(14)= S13(x,y,l);
                        Svec(15)= S14(x,y,l);
                        Svec(16)= S15(x,y,l);
                        Svec(17)= S16(x,y,l);
                        Svec(isnan(Svec)) = 0;
                        Svec = Svec/Svec(1);
                        Svec19 = Svec(1:9);
                        Svec49 = Svec([1,10:end]);

                        % Fit the Sub-diffusion model with 2 Delta
                        SUB = @(params, q) ml(-params(1) * q.^2 .* diff_t.^params(2), params(2)) ;
                        [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], qvecAll, Svec, [Dmin betamin], [Dmax betamax], options);
                        D_SUBq(x,y,l) = params(1)/1e6;
                        beta_SUBq(x,y,l)  = params(2);
                        error_SUBq(x,y,l) = sqrt(resnorm/length(Svec)); %RMSE

                        if ~RedDir
                            % Sub-diffusion model Delta = 19 only
                            SUB = @(params, q) ml(-params(1) * q.^2 .* diff_t(idx19).^params(2), params(2)) ;
                            [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], qvecAll(idx19), Svec(idx19), [Dmin betamin], [Dmax betamax], options);
                            D_SUB19(x,y,l) = params(1)/1e6;
                            beta_SUB19(x,y,l)  = params(2);
                            error_SUB19(x,y,l) = sqrt(resnorm/length(Svec(idx19))); %RMSE


                            % Sub-diffusion model Delta = 49 only
                            SUB = @(params, q) ml(-params(1) * q.^2 .* diff_t(idx49).^params(2), params(2)) ;
                            [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], qvecAll(idx49), Svec(idx49), [Dmin betamin], [Dmax betamax], options);
                            D_SUB49(x,y,l) = params(1)/1e6;
                            beta_SUB49(x,y,l)  = params(2);
                            error_SUB49(x,y,l) = sqrt(resnorm/length(Svec(idx49))); %RMSE


                            % DKI model (fit b=0, 50, 350, 800, 1500, 2400 for Delta=19ms)
                            % (fit b = 0 200 950 2300 for Delta=49ms)
                            nb_DKI19 = 6;
                            nb_DKI49 = 4;

                            DKI = @(params, b) exp(-b*params(1) + 1/6*b.^2.*params(1)^2.*params(2));
                            [params, resnorm] = lsqcurvefit(DKI, [D_guess_DKI K_guess], bvec19, Svec19(1:length(bvec19)), [Dmin_DKI Kmin], [Dmax_DKI, Kmax], options);
                            D_DKI19(x,y,l) = params(1);
                            K_DKI19(x,y,l) = params(2);
                            error_DKI19(x,y,l) = sqrt(resnorm/length(Svec(1:nb_DKI19))); %RMSE

                            [params, resnorm] = lsqcurvefit(DKI, [D_guess_DKI K_guess], bvec49, Svec49(1:length(bvec49)), [Dmin_DKI Kmin], [Dmax_DKI, Kmax], options);
                            D_DKI49(x,y,l) = params(1);
                            K_DKI49(x,y,l) = params(2);
                            error_DKI49(x,y,l) = sqrt(resnorm/length(Svec49(1:nb_DKI49))); %RMSE
                        end

                        % Sub-diffusion model 5 b values best
                        SUB = @(params, q) ml(-params(1) * q.^2 .* diff_t(idx5a).^params(2), params(2)) ;
                        [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], qvecAll(idx5a), Svec(idx5a), [Dmin betamin], [Dmax betamax], options);
                        D_SUB5a(x,y,l) = params(1)/1e6;
                        beta_SUB5a(x,y,l)  = params(2);
                        error_SUB5a(x,y,l) = sqrt(resnorm/length(Svec(idx5a))); %RMSE


                        % Sub-diffusion model 5 b values bad
                        SUB = @(params, q) ml(-params(1) * q.^2 .* diff_t(idx5b).^params(2), params(2)) ;
                        [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], qvecAll(idx5b), Svec(idx5b), [Dmin betamin], [Dmax betamax], options);
                        D_SUB5b(x,y,l) = params(1)/1e6;
                        beta_SUB5b(x,y,l)  = params(2);
                        error_SUB5b(x,y,l) = sqrt(resnorm/length(Svec(idx5b))); %RMSE

                        % Sub-diffusion model 4 b values best
                        SUB = @(params, q) ml(-params(1) * q.^2 .* diff_t(idx4a).^params(2), params(2)) ;
                        [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], qvecAll(idx4a), Svec(idx4a), [Dmin betamin], [Dmax betamax], options);
                        D_SUB4a(x,y,l) = params(1)/1e6;
                        beta_SUB4a(x,y,l)  = params(2);
                        error_SUB4a(x,y,l) = sqrt(resnorm/length(Svec(idx4a))); %RMSE


                        % Sub-diffusion model 4 b values bad
                        SUB = @(params, q) ml(-params(1) * q.^2 .* diff_t(idx4b).^params(2), params(2)) ;
                        [params,resnorm] = lsqcurvefit(SUB, [D_guess beta_guess], qvecAll(idx4b), Svec(idx4b), [Dmin betamin], [Dmax betamax], options);
                        D_SUB4b(x,y,l) = params(1)/1e6;
                        beta_SUB4b(x,y,l)  = params(2);
                        error_SUB4b(x,y,l) = sqrt(resnorm/length(Svec(idx4b))); %RMSE


                    end
                end
                disp(['Subject', SUBJECT, 'direction reduction ', num2str(diridx), ' of ', num2str(length(RED_Dir)), ', l = ',num2str(l),', x = ',num2str(x)])
            end
        end

        toc
        % D* and K* computed from Sub-diffusion model
        K_starq = 6*(gamma(1+beta_SUBq)).^2 ./ gamma(1+2*beta_SUBq) - 3;
        D_starq = D_SUBq./gamma(1+beta_SUBq);

        if ~RedDir
            K_star19 = 6*(gamma(1+beta_SUB19)).^2 ./ gamma(1+2*beta_SUB19) - 3;
            D_star19 = D_SUB19./gamma(1+beta_SUB19);

            K_star49 = 6*(gamma(1+beta_SUB49)).^2 ./ gamma(1+2*beta_SUB49) - 3;
            D_star49 = D_SUB49./gamma(1+beta_SUB49);
        end
        K_star5a = 6*(gamma(1+beta_SUB5a)).^2 ./ gamma(1+2*beta_SUB5a) - 3;
        D_star5a = D_SUB5a./gamma(1+beta_SUB5a);

        K_star5b = 6*(gamma(1+beta_SUB5b)).^2 ./ gamma(1+2*beta_SUB5b) - 3;
        D_star5b = D_SUB5b./gamma(1+beta_SUB5b);

        K_star4a = 6*(gamma(1+beta_SUB4a)).^2 ./ gamma(1+2*beta_SUB4a) - 3;
        D_star4a = D_SUB4a./gamma(1+beta_SUB4a);

        K_star4b = 6*(gamma(1+beta_SUB4b)).^2 ./ gamma(1+2*beta_SUB4b) - 3;
        D_star4b = D_SUB4b./gamma(1+beta_SUB4b);

        % Save Results
        if RedDir
            save(['results_qspace_all_data_real_',SUBJECT,'_D_beta_253mTm_dirRed_32_', num2str(numDir32), '_64_', num2str(numDir64),'_P_D_MM.mat'])
        else
            save(['results_qspace_all_data_real_',SUBJECT,'_D_beta_253mTm.mat'])
        end
    end
end
