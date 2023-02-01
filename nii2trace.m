%% This file creates the trace weighted nii file for all data.
%% transfer nii to trace nii data

PATH   = pwd;
subjects = {'sub_001'};%, 'sub_003', 'sub_004', 'sub_005', 'sub_006', 'sub_007', 'sub_001_rescan', 'sub_003_rescan', 'sub_004_rescan', 'sub_005_rescan', 'sub_006_rescan', 'sub_007_rescan'};

for ii = 1:length(subjects)
    SUBJECT    = subjects{ii};
    for DELTA = [19; 49]
        FOLDER = [PATH '\' SUBJECT '\' 'dwi_real\'];
        filename = [FOLDER, SUBJECT,'_dwi_real.nii.gz'];
        data = niftiread(filename);
        bvals = load([FOLDER, SUBJECT, '_dwi_real.bval']);
        bvec = load([FOLDER, SUBJECT, '_dwi_real.bvec']);
        delta = load([FOLDER, SUBJECT, '_dwi_real.delta']);

        bval = unique(bvals);
        num_bvalues = length(bval); % including b0

        
        [N,M,L,b_dirs] = size(data);
        num_bvals = 8; % per delta value
        trace = zeros(N,M,L,num_bvalues);

        data(data<0) = nan;

        for bs = 1:num_bvalues
            trace(:,:,:,bs) = geomean(abs(data(:,:,:,bvals==bval(bs))),4,'omitnan');
        end

        niftiwrite(squeeze(trace(:,:,:,1)), [FOLDER, SUBJECT,'_dwi_real_b0_delta0_image.nii'])

        for bs = 2:num_bvalues
            if unique(delta(bvals==bval(bs)))==19
                niftiwrite(squeeze(trace(:,:,:,bs)), [FOLDER, SUBJECT,'_dwi_real_b', num2str(bval(bs)),'_delta19_image.nii'])
            elseif unique(delta(bvals==bval(bs)))==49
                niftiwrite(squeeze(trace(:,:,:,bs)), [FOLDER, SUBJECT,'_dwi_real_b', num2str(bval(bs)),'_delta49_image.nii'])
            end
        end

    end
end
