%% This file generates all figures and tables of results to the fitting of the connectome data.
% Be sure to run FitConnectome first.

%% Create Brain figures
PATH   = pwd;
fontsize = 10;
diff_time19 = (19 - 8/3)/1000;
diff_time49 = (49 - 8/3)/1000;
image1Sub = 'sub_003_rescan';
image1slice = 71;
image2Sub = 'sub_005';
image2slice = 74;

brainmask = double(niftiread([PATH '\' image1Sub '\' 'dwi_real\' image1Sub '_dwi_real_brainmask.nii.gz']));
mask = squeeze(brainmask(:,:,image1slice));
Height = size(mask,2);
Width = size(mask,1);

xmin = max(find(sum(mask,2)>1,1, 'first')-2,1);
xmax = min(find(sum(mask,2)>1,1, 'last')+2, Width);
ymin = max(find(sum(mask,1)>1,1, 'first')-3,1);
ymax = min(find(sum(mask,1)>1,1, 'last')+2,Height);

Height = ymax-ymin+1;
Width = xmax-xmin+1;
mask = mask(xmin:xmax,ymin:ymax);

xmin2 = xmin + 7;
xmax2 = xmax - 7;
ymin2 = ymin + 6;
ymax2 = ymax - 12;

Height2 = ymax2 - ymin2+1;
Width2 = xmax2 - xmin2+1;

brainmask = double(niftiread([PATH '\' image2Sub '\' 'dwi_real\' image2Sub '_dwi_real_brainmask.nii.gz']));
mask2 = squeeze(brainmask(:,:,image2slice));
mask2 = mask2(xmin2:xmax2, ymin2:ymax2);

resAll_1 = load(['results_qspace_all_data_real_',image1Sub,'_D_beta_253mTm.mat'], 'K_starq', 'K_star19', 'K_star49', 'K_DKI19', 'K_DKI49', 'K_star5a', 'K_star5b', 'beta_SUBq', 'D_SUBq');
res16_1 = load(['results_qspace_all_data_real_',image1Sub,'_D_beta_253mTm_dirRed_32_16_64_16_P_D_MM.mat'], 'K_starq','K_star5a', 'K_star5b');
res4_1 = load(['results_qspace_all_data_real_',image1Sub,'_D_beta_253mTm_dirRed_32_4_64_4P_D_MM.mat'], 'K_starq', 'K_star5a', 'K_star5b');


resAll_2 = load(['results_qspace_all_data_real_',image2Sub,'_D_beta_253mTm.mat'], 'K_starq', 'K_star19', 'K_star49', 'K_DKI19', 'K_DKI49', 'K_star5a', 'K_star5b', 'beta_SUBq', 'D_SUBq');
res16_2 = load(['results_qspace_all_data_real_',image2Sub,'_D_beta_253mTm_dirRed_32_16_64_16_P_D_MM.mat'], 'K_starq', 'K_star5a', 'K_star5b');
res4_2 = load(['results_qspace_all_data_real_',image2Sub,'_D_beta_253mTm_dirRed_32_4_64_4P_D_MM.mat'], 'K_starq', 'K_star5a', 'K_star5b');


figure
set(gcf, 'position', [100,50,524,715])
tileHeight = 4;

t = tiledlayout(5,4*tileHeight, 'tilespacing', 'none', 'padding', 'compact');%,'TileIndexing','columnmajor'

padding = 6;
nexttile([1,tileHeight])
imagesc([squeeze(resAll_1.D_SUBq(xmin:xmax,ymin:ymax,image1slice)).*mask, zeros(Width, 15)])
axis image
clim([0,1e-3])
ylim([0.5,Width+0.5+padding])
xlim([0.5,Height+0.5+7.5])
view(-90,90)
axis off
set(gca, 'ydir', 'normal', 'Clipping','off')
c2 = colorbar('ticklabelinterpreter', 'latex');
c2.Label.String = '$D_\beta$';
c2.Label.Interpreter = 'latex';
c2.Label.FontSize = fontsize;
c2.Layout.Tile = 'East';
text(Height/2, Width+0.5+5, '$\Delta = 19, 49\textrm{ ms}$', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex', 'Rotation',90)
text( Height+0.5+6,Width/2,'$D_\beta$', 'FontSize',fontsize+2, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','color','w')
text( Height+0.5+15,Width,'(A)', 'FontSize',fontsize+2, 'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex')

nexttile([1,tileHeight])
imagesc([squeeze(resAll_1.beta_SUBq(xmin:xmax,ymin:ymax,image1slice)).*mask, zeros(Width, 15)])
axis image
clim([0.5,1])
ylim([0.5-padding,Width+0.5])
xlim([0.5,Height+0.5+7.5])
view(-90,90)
axis off
set(gca, 'ydir', 'normal', 'Clipping','off')
text( Height+0.5+6,Width/2,'$\beta$', 'FontSize',fontsize+2, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','color','w')

nexttile([1,tileHeight])
imagesc([squeeze(resAll_2.D_SUBq(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, zeros(Width2, 12)])
axis image
clim([0,1e-3])
ylim([0.5,Width2+0.5+0.8*padding])
xlim([0.5,Height2+0.5+6])
view(-90,90)
axis off
set(gca, 'ydir', 'normal', 'Clipping','off')
text( Height2+0.5+0.8*6,Width2/2,'$D_\beta$', 'FontSize',fontsize+2, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','color','w', 'FontWeight','bold')
text( Height2+0.5+0.8*15,Width2,'(B)', 'FontSize',fontsize+2, 'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','latex')

nexttile([1,tileHeight])
imagesc([squeeze(resAll_2.beta_SUBq(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, zeros(Width2, 12)])
axis image
clim([0.5,1])
ylim([0.5-0.8*padding,Width2+0.5])
xlim([0.5,Height2+0.5+6])
view(-90,90)
axis off
set(gca, 'ydir', 'normal', 'Clipping','off')
text( Height2+0.5+0.8*6,Width2/2,'$\beta$', 'FontSize',fontsize+2, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','color','w', 'FontWeight','bold')

c3 = colorbar('ticklabelinterpreter', 'latex');
c3.Label.String = '$\beta$';
c3.Label.Interpreter = 'latex';
c3.Layout.Tile = 'East';
c3.Label.FontSize = fontsize;

nexttile(t,[2,tileHeight*2])
imagesc([squeeze(resAll_1.K_starq(xmin:xmax,ymin:ymax,image1slice)).*mask, zeros(Width,9)])
axis image
clim([0,3])
ylim([0.5-0.5*padding, Width+0.5+0.5*padding])
xlim([0.5, Height+0.5+7.5])

view(-90,90)
text(Height/2, Width+0.5+2.5, '$\Delta = 19, 49\textrm{ ms}$', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex', 'Rotation',90)
text(Height+0.5+2.5, Width/2, '$K^*$', 'FontSize',fontsize+2, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','color','w')
axis off
set(gca, 'ydir', 'normal', 'Clipping','off')


nexttile(t,[2,tileHeight*2])
imagesc([squeeze(resAll_2.K_starq(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, zeros(Width2,7)])
axis image
clim([0,3])
ylim([0.5-0.4*padding, Width2+0.5+0.4*padding])
xlim([0.5, Height2+0.5+6])
view(-90,90)
text(Height2+0.5+2, Width2/2, '$K^*$', 'FontSize',fontsize+2, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','color','w')
axis off
set(gca, 'ydir', 'normal', 'Clipping','off')



nexttile(t,[2,tileHeight])
imagesc([squeeze(resAll_1.K_star49(xmin:xmax,ymin:ymax,image1slice)).*mask, squeeze(resAll_1.K_star19(xmin:xmax,ymin:ymax,image1slice)).*mask, zeros(Width,17)])
axis image
clim([0,3])
ylim([0.5,Width+0.5+padding])
xlim([0.5, 2*Height+0.5+15])
axis off
colormap hot
text(Height*2+0.5 + 7, Width/2,'$K^*$', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','latex', 'fontsize', fontsize+2,'color','w')
text( Height/2,Width+5,'$\Delta = 49\textrm{ ms}$', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','rotation',90)
text(Height/2+Height, Width+5,'$\Delta = 19 \textrm{ ms}$', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','rotation',90)
view(-90,90)
set(gca, 'ydir', 'normal', 'Clipping','off')



nexttile(t,[2,tileHeight])
imagesc([squeeze(resAll_1.K_DKI49(xmin:xmax,ymin:ymax,image1slice)).*mask, squeeze(resAll_1.K_DKI19(xmin:xmax,ymin:ymax,image1slice)).*mask, zeros(Width,17)])
view(-90,90)
axis image
clim([0,3])
ylim([0.5-padding,Width+0.5])
xlim([0.5, 2*Height+0.5+15])
axis off
text(Height*2+0.5 + 7, Width/2,'$K_{DKI}$', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','latex', 'fontsize', fontsize+2,'color','w')
set(gca, 'ydir', 'normal', 'Clipping','off')

nexttile(t,[2,tileHeight])
imagesc([squeeze(resAll_2.K_star49(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, squeeze(resAll_2.K_star19(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, zeros(Width2,14)])
axis image
clim([0,3])

ylim([0.5,Width2+0.5+0.8*padding])
xlim([0.5, 2*Height2+0.5+12])
view(-90,90)
axis off
colormap hot
text(Height2*2+0.5 + 5.6, Width2/2,'$K^*$', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','latex', 'fontsize', fontsize+2,'color','w')
set(gca, 'ydir', 'normal', 'Clipping','off')

nexttile(t,[2,tileHeight])
imagesc([squeeze(resAll_2.K_DKI49(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, squeeze(resAll_2.K_DKI19(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, zeros(Width2,14)])
axis image
clim([0,3])

ylim([0.5-0.8*padding,Width2+0.5])
xlim([0.5, 2*Height2+0.5+12])
view(-90,90)
axis off
text(Height2*2+0.5 + 5.6, Width2/2,'$K_{DKI}$', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','latex', 'fontsize', fontsize+2,'color','w')
set(gca, 'ydir', 'normal', 'Clipping','off')
c = colorbar;
c.Layout.Tile = 'East';
c.TickLabelInterpreter = 'latex';
c.FontSize = fontsize;
c.Label.Interpreter = 'latex';
c.Label.String = '$K^*,\, K_{DKI}$';%'Kurtosis - $D_\beta (\times 10^{-3})$ - $\beta$';
c2.FontSize = fontsize;
c3.FontSize = fontsize;
c.Label.FontSize = fontsize;
c2.Label.FontSize = fontsize;
c3.Label.FontSize = fontsize;

exportgraphics(gcf,'BrainImages1.eps')

l = image1slice;
figure
set(gcf, 'position', [100,100,625,550])
tileHeight = 4;
t = tiledlayout(3, 4*tileHeight+1, 'tilespacing', 'none', 'padding', 'loose', 'TileIndexing','columnmajor');
nexttile([3,tileHeight])
imagesc([squeeze(resAll_1.K_star5a(xmin:xmax,ymin:ymax,image1slice)).*mask, squeeze(res16_1.K_star5a(xmin:xmax,ymin:ymax,image1slice)).*mask, squeeze(res4_1.K_star5a(xmin:xmax,ymin:ymax,image1slice)).*mask])
axis image
clim([0,3])
view(-90,90)
axis off
colormap hot
text(Height*3 +0.5 + 5, Width/2,'Optimal $b$ values', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','latex', 'fontsize', fontsize)
text(Height*3 +0.5 + 14, Width,'(A)', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Interpreter','latex', 'fontsize', fontsize+2)
text(Height/2, Width+0.5+5, 'SNR = 20', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex', 'Rotation',90)
text(Height/2+Height, Width+0.5+5, 'SNR = 10', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex', 'Rotation',90)
text(Height/2+2*Height, Width+0.5+5, 'SNR = 5', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex', 'Rotation',90)
set(gca, 'ydir', 'normal', 'Clipping','off')


nexttile([3,tileHeight])
imagesc([squeeze(resAll_1.K_star5b(xmin:xmax,ymin:ymax,image1slice)).*mask, squeeze(res16_1.K_star5b(xmin:xmax,ymin:ymax,image1slice)).*mask, squeeze(res4_1.K_star5b(xmin:xmax,ymin:ymax,image1slice)).*mask])
axis image
clim([0,3])
view(-90,90)
axis off
text(Height*3+0.5+5,Width/2, 'Suboptimal $b$ values', 'HorizontalAlignment','center', 'VerticalAlignment','middle',  'Interpreter','latex', 'fontsize', fontsize)
set(gca, 'ydir', 'normal', 'Clipping','off')

nexttile([3,1])
axis off

nexttile([3,tileHeight])
imagesc([squeeze(resAll_2.K_star5a(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, squeeze(res16_2.K_star5a(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, squeeze(res4_2.K_star5a(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2])
axis image
clim([0,3])
view(-90,90)
axis off
colormap hot
text(Height2*3+0.5 + 4, Width2/2,'Optimal $b$ values', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','latex', 'fontsize', fontsize)
text(Height2*3 +0.5 + 0.8*14, Width2,'(B)', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Interpreter','latex', 'fontsize', fontsize+2)
text(Height2/2, Width2+0.5+4, 'SNR = 20', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex', 'Rotation',90)
text(Height2/2+Height2, Width2+0.5+4, 'SNR = 10', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex', 'Rotation',90)
text(Height2/2+2*Height2, Width2+0.5+4, 'SNR = 5', 'FontSize',fontsize, 'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex', 'Rotation',90)
set(gca, 'ydir', 'normal', 'Clipping','off')

nexttile([3,tileHeight])
imagesc([squeeze(resAll_2.K_star5b(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, squeeze(res16_2.K_star5b(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2, squeeze(res4_2.K_star5b(xmin2:xmax2,ymin2:ymax2,image2slice)).*mask2])
axis image
clim([0,3])
view(-90,90)
axis off
text(Height2*3+5,Width2/2, 'Suboptimal $b$ values', 'HorizontalAlignment','center', 'VerticalAlignment','middle',  'Interpreter','latex', 'fontsize', fontsize)
set(gca, 'ydir', 'normal', 'Clipping','off')
c = colorbar;
c.Layout.Tile = 'South';
c.TickLabelInterpreter = 'latex';
c.FontSize = fontsize;

exportgraphics(gcf,'BrainImages2.eps')

%% Create Table 1
PATH   = pwd;

diff_time19 = (19 - 8/3)/1000;
diff_time49 = (49 - 8/3)/1000;
SUBJECTs    = {'sub_001', 'sub_003', 'sub_004','sub_005','sub_006','sub_007','sub_001_rescan', 'sub_003_rescan','sub_004_rescan','sub_005_rescan','sub_006_rescan','sub_007_rescan'};
regions = {'scGM', 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'cGM', 'Fusiform', 'Lingual', 'WM', 'Cerebral WM', 'Cerebellum WM', 'CC'};
percent = 10;
N = length(SUBJECTs);
M = length(regions);

mubetaAll = zeros(M, N);
sigma2betaAll = zeros(M, N);
nibetaAll = zeros(M, N);
muDbetaAll = zeros(M, N);
sigma2DbetaAll = zeros(M, N);
niDbetaAll = zeros(M, N);

muKstarAll = zeros(M, N);
muKDKI19All = zeros(M, N);
muKDKI49All = zeros(M, N);

sigma2KstarAll = zeros(M, N);
sigma2KDKI19All = zeros(M, N);
sigma2KDKI49All = zeros(M, N);

muKstar5aAll = zeros(M, N);
muKstar5bAll = zeros(M, N);

sigma2Kstar5aAll = zeros(M, N);
sigma2Kstar5bAll = zeros(M, N);

muKstar5a16 = zeros(M, N);
muKstar5b16 = zeros(M, N);

sigma2Kstar5a16 = zeros(M, N);
sigma2Kstar5b16 = zeros(M, N);

muKstar5a4 = zeros(M, N);
muKstar5b4 = zeros(M, N);

sigma2Kstar5a4 = zeros(M, N);
sigma2Kstar5b4 = zeros(M, N);

niKstarAll = zeros(M, N);
niKDKI19All = zeros(M, N);
niKDKI49All = zeros(M, N);

niKstar5aAll = zeros(M, N);
niKstar5bAll = zeros(M, N);

niKstar5a16 = zeros(M, N);
niKstar5b16 = zeros(M, N);

niKstar5a4 = zeros(M, N);
niKstar5b4 = zeros(M, N);

for ii = 1:N
    SUBJECT = SUBJECTs{ii};
    FOLDER = [PATH '\' SUBJECT '\' 'dwi_real\'];
    segmentation = niftiread([FOLDER SUBJECT '_dwi_real_aparc+aseg.nii.gz']);
    brainmask = logical(niftiread([FOLDER SUBJECT '_dwi_real_brainmask.nii.gz']));
    segmentation(~brainmask) = 0;
    resAll = load(['results_qspace_all_data_real_',SUBJECT,'_D_beta_253mTm.mat'], 'K_starq',  'K_DKI19',  'K_DKI49', 'K_star5a', 'K_star5b', 'beta_SUB5a', 'D_SUB5a');
    res16_1 = load(['results_qspace_all_data_real_',SUBJECT,'_D_beta_253mTm_dirRed_32_16_64_16_P_D_MM.mat'], 'K_star5a', 'K_star5b');
    res4_1 = load(['results_qspace_all_data_real_',SUBJECT,'_D_beta_253mTm_dirRed_32_4_64_4P_D_MM.mat'], 'K_star5a', 'K_star5b');
    for jj = 1:M
        REGION = regions{jj};
        map = BrainRegions(segmentation, REGION);
        [mubetaAll(jj,ii), sigma2betaAll(jj,ii), nibetaAll(jj,ii)] = SummaryStats221011(resAll.beta_SUB5a(map),percent);
        [muDbetaAll(jj,ii), sigma2DbetaAll(jj,ii), niDbetaAll(jj,ii)] = SummaryStats221011(resAll.D_SUB5a(map),percent);

        [muKstarAll(jj,ii), sigma2KstarAll(jj,ii), niKstarAll(jj,ii)] = SummaryStats221011(resAll.K_starq(map),percent);
        [muKDKI19All(jj,ii), sigma2KDKI19All(jj,ii), niKDKI19All(jj,ii)] = SummaryStats221011(resAll.K_DKI19(map),percent);
        [muKDKI49All(jj,ii), sigma2KDKI49All(jj,ii), niKDKI49All(jj,ii)] = SummaryStats221011(resAll.K_DKI49(map),percent);

        [muKstar5aAll(jj,ii), sigma2Kstar5aAll(jj,ii), niKstar5aAll(jj,ii)] = SummaryStats221011(resAll.K_star5a(map), percent);
        [muKstar5a16(jj,ii), sigma2Kstar5a16(jj,ii), niKstar5a16(jj,ii)] = SummaryStats221011(res16_1.K_star5a(map), percent);
        [muKstar5a4(jj,ii), sigma2Kstar5a4(jj,ii), niKstar5a4(jj,ii)] = SummaryStats221011(res4_1.K_star5a(map), percent);
        [muKstar5bAll(jj,ii), sigma2Kstar5bAll(jj,ii), niKstar5bAll(jj,ii)] = SummaryStats221011(resAll.K_star5b(map), percent);
        [muKstar5b16(jj,ii), sigma2Kstar5b16(jj,ii), niKstar5b16(jj,ii)] = SummaryStats221011(res16_1.K_star5b(map), percent);
        [muKstar5b4(jj,ii), sigma2Kstar5b4(jj,ii), niKstar5b4(jj,ii)] = SummaryStats221011(res4_1.K_star5b(map), percent);
    end
end

muSigmaCVbeta = [sum(mubetaAll.*nibetaAll,2)./sum(nibetaAll,2), sqrt(sum(sigma2betaAll.*(nibetaAll-1),2)./(sum(nibetaAll,2)-N)), zeros(M,1)];
muSigmaCVbeta(:,3) = muSigmaCVbeta(:,2)./muSigmaCVbeta(:,1)*100;
muSigmaCVDbeta = [sum(muDbetaAll.*niDbetaAll,2)./sum(niDbetaAll,2), sqrt(sum(sigma2DbetaAll.*(niDbetaAll-1),2)./(sum(niDbetaAll,2)-N)), zeros(M,1)];
muSigmaCVDbeta(:,3) = muSigmaCVDbeta(:,2)./muSigmaCVDbeta(:,1)*100;

muSigmaCVKstar = [sum(muKstarAll.*niKstarAll,2)./sum(niKstarAll,2), sqrt(sum(sigma2KstarAll.*(niKstarAll-1),2)./(sum(niKstarAll,2)-N)), zeros(M,1)];
muSigmaCVKstar(:,3) = muSigmaCVKstar(:,2)./muSigmaCVKstar(:,1)*100;


muSigmaCVKDKI19 = [sum(muKDKI19All.*niKDKI19All,2)./sum(niKDKI19All,2), sqrt(sum(sigma2KDKI19All.*(niKDKI19All-1),2)./(sum(niKDKI19All,2)-N)), zeros(M,1)];
muSigmaCVKDKI19(:,3) = muSigmaCVKDKI19(:,2)./muSigmaCVKDKI19(:,1)*100;


muSigmaCVKDKI49 = [sum(muKDKI49All.*niKDKI49All,2)./sum(niKDKI49All,2), sqrt(sum(sigma2KDKI49All.*(niKDKI49All-1),2)./(sum(niKDKI49All,2)-N)), zeros(M,1)];
muSigmaCVKDKI49(:,3) = muSigmaCVKDKI49(:,2)./muSigmaCVKDKI49(:,1)*100;

muSigmaCVKstar5aAll = [sum(muKstar5aAll.*niKstar5aAll,2)./sum(niKstar5aAll,2), sqrt(sum(sigma2Kstar5aAll.*(niKstar5aAll-1),2)./(sum(niKstar5aAll,2)-N)), zeros(M,1)];
muSigmaCVKstar5aAll(:,3) = muSigmaCVKstar5aAll(:,2)./muSigmaCVKstar5aAll(:,1)*100;

muSigmaCVKstar5a16 = [sum(muKstar5a16.*niKstar5a16,2)./sum(niKstar5a16,2), sqrt(sum(sigma2Kstar5a16.*(niKstar5a16-1),2)./(sum(niKstar5a16,2)-N)), zeros(M,1)];
muSigmaCVKstar5a16(:,3) = muSigmaCVKstar5a16(:,2)./muSigmaCVKstar5a16(:,1)*100;


muSigmaCVKstar5a4 = [sum(muKstar5a4.*niKstar5a4,2)./sum(niKstar5a4,2), sqrt(sum(sigma2Kstar5a4.*(niKstar5a4-1),2)./(sum(niKstar5a4,2)-N)), zeros(M,1)];
muSigmaCVKstar5a4(:,3) = muSigmaCVKstar5a4(:,2)./muSigmaCVKstar5a4(:,1)*100;

muSigmaCVKstar5bAll = [sum(muKstar5bAll.*niKstar5bAll,2)./sum(niKstar5bAll,2), sqrt(sum(sigma2Kstar5bAll.*(niKstar5bAll-1),2)./(sum(niKstar5bAll,2)-N)), zeros(M,1)];
muSigmaCVKstar5bAll(:,3) = muSigmaCVKstar5bAll(:,2)./muSigmaCVKstar5bAll(:,1)*100;

muSigmaCVKstar5b16 = [sum(muKstar5b16.*niKstar5b16,2)./sum(niKstar5b16,2), sqrt(sum(sigma2Kstar5b16.*(niKstar5b16-1),2)./(sum(niKstar5b16,2)-N)), zeros(M,1)];
muSigmaCVKstar5b16(:,3) = muSigmaCVKstar5b16(:,2)./muSigmaCVKstar5b16(:,1)*100;


muSigmaCVKstar5b4 = [sum(muKstar5b4.*niKstar5b4,2)./sum(niKstar5b4,2), sqrt(sum(sigma2Kstar5b4.*(niKstar5b4-1),2)./(sum(niKstar5b4,2)-N)), zeros(M,1)];
muSigmaCVKstar5b4(:,3) = muSigmaCVKstar5b4(:,2)./muSigmaCVKstar5b4(:,1)*100;

for ii = 1:M
    REGION = regions{ii};
    if strcmpi(REGION, 'scGM') || strcmpi(REGION, 'cGM') || strcmpi(REGION, 'WM')
        REGION = ['\\textbf{', REGION, '}'];
    end
    fprintf([REGION, ' & $ %.2f \\pm %.2f (%.0f \\%%)$  & $ %.2f \\pm %.2f (%.0f \\%%)$ & $ %.2f \\pm %.2f (%.0f \\%%)$ \\\\ \n'], [muSigmaCVKDKI19(ii,:), muSigmaCVKDKI49(ii,:), muSigmaCVKstar(ii,:)])
end
fprintf(['Optimal $b$ values \n'])
for ii = 1:M
    REGION = regions{ii};
    if strcmpi(REGION, 'scGM') || strcmpi(REGION, 'cGM') || strcmpi(REGION, 'WM')
        REGION = ['\\textbf{', REGION, '}'];
    end
    fprintf([REGION, ' & $ \\textit{ %.2f } \\pm %.2f (%.0f \\%%)$  & $ %.2f \\pm %.2f (%.0f \\%%)$ & $ %.2f \\pm %.2f (%.0f \\%%)$ \\\\ \n'], [muSigmaCVKstar5a4(ii,:), muSigmaCVKstar5a16(ii,:), muSigmaCVKstar5aAll(ii,:)])
end

fprintf(['Suboptimal $b$ values \n'])
for ii = 1:M
    REGION = regions{ii};
    if strcmpi(REGION, 'scGM') || strcmpi(REGION, 'cGM') || strcmpi(REGION, 'WM')
        REGION = ['\\textbf{', REGION, '}'];
    end
    fprintf([REGION, ' & $ %.2f \\pm %.2f (%.0f \\%%)$  & $ %.2f \\pm %.2f (%.0f \\%%)$ & $ %.2f \\pm %.2f (%.0f \\%%)$  \\\\ \n'], [muSigmaCVKstar5b4(ii,:), muSigmaCVKstar5b16(ii,:), muSigmaCVKstar5bAll(ii,:)])
end

%% Table 2
PATH = pwd;
regions = {'scGM', 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'cGM', 'Fusiform', 'Lingual', 'WM', 'Cerebral WM', 'Cerebellum WM', 'CC'};
subjects = [1,3, 4,5,6,7];
p = zeros(length(regions), length(subjects));
p2 = zeros(length(regions), length(subjects));
p3 = zeros(length(regions), length(subjects));
mibar = zeros(length(regions), length(subjects));
intrai = zeros(length(regions), length(subjects));
data1cell = cell(size(regions));
data2cell = cell(size(regions));
segmentation1 = niftiread([PATH, '\sub_001\dwi_real\sub_001_dwi_real_aparc+aseg.nii.gz']);
size1 = size(segmentation1);
regionSize = zeros(size(regions));
Data = cell(size(regions));
Data2 = cell(size(regions));
for ii = 1:length(regions)
REGION = regions{ii};
map = BrainRegions(segmentation1, REGION);
regionSize(ii) = sum(map(:));
Data{ii} = zeros(length(subjects)*2,regionSize(ii));
end

for ii = 1:length(subjects)
    SUBJECT    = ['sub_00', num2str(subjects(ii))];
    FOLDER = [PATH, '\', SUBJECT,'\dwi_real\'];
    FOLDER2 = [PATH '\' SUBJECT '_rescan\' 'dwi_real\'];
    segmentation = niftiread([FOLDER SUBJECT '_dwi_real_aparc+aseg.nii.gz']);
    brainmask = logical(niftiread([FOLDER SUBJECT '_dwi_real_brainmask.nii.gz']));
    segmentation(~brainmask) = 0;

    segmentation2 = niftiread([FOLDER2 SUBJECT '_rescan_dwi_real_aparc+aseg.nii.gz']);
    brainmask2 = logical(niftiread([FOLDER2 SUBJECT '_rescan_dwi_real_brainmask.nii.gz']));
    segmentation2(~brainmask2) = 0;

    resAll = load(['results_qspace_all_data_real_',SUBJECT,'_D_beta_253mTm.mat'], 'K_starq',  'K_DKI19',  'K_DKI49', 'K_star5a');
    res16_1 = load(['results_qspace_all_data_real_',SUBJECT,'_D_beta_253mTm_dirRed_32_16_64_16_P_D_MM.mat'], 'K_star5a', 'K_star4a');
    resAllrescan = load(['results_qspace_all_data_real_',SUBJECT,'_rescan_D_beta_253mTm.mat'], 'K_starq',  'K_DKI19',  'K_DKI49', 'K_star5a');
    res16rescan_1 = load(['results_qspace_all_data_real_',SUBJECT,'_rescan_D_beta_253mTm_dirRed_32_16_64_16_P_D_MM.mat'], 'K_star5a', 'K_star4a');

    load([FOLDER2, SUBJECT, '_rescan_image_registration_Transformation_to_sub_001.mat'])
    tformSimilarity_rescan = tformSimilarity;
    load([FOLDER, SUBJECT, '_image_registration_Transformation_to_sub_001.mat'])

    centerOutput = affineOutputView(size1,tformSimilarity,'BoundsStyle','sameAsInput');
    resAll.K_starqReg = imwarp(resAll.K_starq, tformSimilarity, 'cubic', 'OutputView', centerOutput);
    resAll.K_star5aReg = imwarp(res16_1.K_star5a, tformSimilarity, 'cubic', 'OutputView', centerOutput);
    resAll.K_DKI19Reg = imwarp(resAll.K_DKI19, tformSimilarity, 'cubic', 'OutputView', centerOutput);
    resAll.K_DKI49Reg = imwarp(resAll.K_DKI49, tformSimilarity, 'cubic', 'OutputView', centerOutput);

    centerOutput = affineOutputView(size1,tformSimilarity_rescan,'BoundsStyle','sameAsInput');
    resAllrescan.K_starqReg = imwarp(resAllrescan.K_starq, tformSimilarity_rescan, 'cubic', 'OutputView', centerOutput);
    resAllrescan.K_star5aReg = imwarp(res16rescan_1.K_star5a, tformSimilarity_rescan, 'cubic', 'OutputView', centerOutput);
    resAllrescan.K_DKI19Reg = imwarp(resAllrescan.K_DKI19, tformSimilarity_rescan, 'cubic', 'OutputView', centerOutput);
    resAllrescan.K_DKI49Reg = imwarp(resAllrescan.K_DKI49, tformSimilarity_rescan, 'cubic', 'OutputView', centerOutput);

    for jj = 1:length(regions)
        REGION = regions{jj};
        map1 = BrainRegions(segmentation1, REGION);
        map = BrainRegions(segmentation, REGION);
        map2 = BrainRegions(segmentation2, REGION);
        Data{jj}(ii*2-1,:) = resAll.K_star5aReg(map1);
        Data{jj}(ii*2,:) = resAllrescan.K_star5aReg(map1);
        Data2{jj}(ii*2-1,:) = resAll.K_starqReg(map1);
        Data2{jj}(ii*2,:) = resAllrescan.K_starqReg(map1);
        
    end


end
%%
meanICC = zeros(size(regions));
stdICC = zeros(size(regions));
figure
set(gcf, 'position', [100,100,585,650])
t=tiledlayout('flow', 'tilespacing', 'tight', 'padding','none');
colour = lines(length(subjects));
for ii = 1:length(regions)
    REGION = regions{ii};
data = Data{ii};
data2 = Data2{ii};

Mi = (data(1:2:end,:) + data(2:2:end,:))/2;
Sinter = var(Mi,0,1);
intraData = data;
intraData(1:2:end,:) = intraData(1:2:end,:) - Mi;
intraData(2:2:end,:) = intraData(2:2:end,:) - Mi;
Sintra = mean(intraData.^2, 1);
ICC = Sinter./(Sinter+Sintra);
Mi2 = (data2(1:2:end,:) + data2(2:2:end,:))/2;
Sinter2 = var(Mi2,0,1);
intraData2 = data2;
intraData2(1:2:end,:) = intraData2(1:2:end,:) - Mi2;
intraData2(2:2:end,:) = intraData2(2:2:end,:) - Mi2;
Sintra2 = mean(intraData2.^2, 1);
ICC2 = Sinter2./(Sinter2+Sintra2);
nexttile(t)
histogram(ICC, 0:0.05:1, 'normalization', 'pdf')
title(REGION, 'interpreter', 'latex')
set(gca, 'ticklabelinterpreter', 'latex')
axis([0,1,0,6])
meanICC(ii) = mean(ICC, 'omitnan');
stdICC(ii) = std(ICC, 'omitnan');
meanICC2(ii) = mean(ICC2, 'omitnan');
stdICC2(ii) = std(ICC2, 'omitnan');
text(0.1,5.5, ['$\mu_A = ', num2str(meanICC2(ii),'%.2f'),'$'], 'interpreter', 'latex')
text(0.1,4.5, ['$\sigma_A = ', num2str(stdICC2(ii),'%.2f'),'$'], 'interpreter', 'latex')
text(0.1,3, ['$\mu_O = ', num2str(meanICC(ii),'%.2f'),'$'], 'interpreter', 'latex')
text(0.1,2, ['$\sigma_O = ', num2str(stdICC(ii),'%.2f'),'$'], 'interpreter', 'latex')
end
exportgraphics(gcf, 'ICC_values.eps')
