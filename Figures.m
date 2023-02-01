%% This file creates the figures and table data for the simulated experiments
% Figures 1-4 and Table 1

%% Figure 1
% make sure you have run Simulation_varyDelta_SetParameters with the
% correct parameters first.

fontsize=10;
Variable = {'Dfit', 'beta_fit', 'K_fit'};
SNR = [5,10,20];
Dvec = [3e-4,5e-4];
betavec = [0.75,0.85];
K = 6*(gamma(1+betavec)).^2 ./ gamma(1+2*betavec) - 3;
alpha = 0.25;

figure
set(gcf, 'Position',[100,100,600,600])
t1 = tiledlayout(4,3,"TileSpacing","tight", "Padding","tight");

numDeltaLabels = {'$1$', '$2$', '$3$', '$4$', '$5$', '$2^\prime$'};
numDeltaAxisPos = [1:5,7]; % separate 2'
numDeltaTickLabels = {'$1$','$2$','$3$', '$4$','$5$', '$2^\prime$'};

CVK = zeros(length(SNR)*length(Dvec),length(numDeltaLabels));
CVBeta = zeros(length(SNR)*length(Dvec),length(numDeltaLabels));
CVD = zeros(length(SNR)*length(Dvec),length(numDeltaLabels));

colour = lines(2);

for ii = 1:length(SNR)
    nexttile
    for ij = 1:length(Dvec)
        load(['VaryDeltaSim_SNR_', num2str(SNR(ii)), '_D_', num2str(Dvec(ij)*1e4),'_beta_',num2str(betavec(ij)*100)], 'Dfit','D')
        bw = zeros(1,size(Dfit,2));
        hold on
        F = zeros(100, size(Dfit,2));
        Dvals = zeros(100, size(Dfit,2));
        for ik = 1:size(Dfit,2)
            [~,edges] = histcounts(Dfit(:,ik)); % use histcounts to find optimal bandwidth for kernel density function
            bw(ik) = edges(2)-edges(1);

            [F(:,ik), Dvals(:,ik)]=ksdensity(Dfit(:,ik),'bandwidth',bw(ik));


            CVD(ik,ij+(ii-1)*length(Dvec)) = sqrt(var(Dfit(:,ik)))/mean(Dfit(:,ik))*100; % Compute coefficient of variation percentage
        end

        F = F/max(F(:))*0.45; % Scale the density distributions so they fit without overlapping.
        for ik = 1:size(Dfit,2)
            fill([F(:,ik)+numDeltaAxisPos(ik); flipud(numDeltaAxisPos(ik)-F(:,ik))], [Dvals(:,ik); flipud(Dvals(:,ik))], colour(ij,:), 'FaceAlpha', alpha, 'edgecolor', colour(ij,:)) % Plot violin
            meanD = mean(Dfit(:,ik));
            meanF = interp1(Dvals(:,ik), F(:,ik), meanD);
            plot([numDeltaAxisPos(ik)-meanF, numDeltaAxisPos(ik)+meanF],[meanD, meanD], '-', 'linewidth', 2, 'color',colour(ij,:))
        end

        plot([0,8],[D,D],':', 'color', colour(ij,:), 'linewidth', 2)
        xlabel('No. of diffusion times', 'interpreter', 'latex', 'fontsize', fontsize)
        ylabel([string(['$\textbf{SNR } \mathbf{= ', num2str(SNR(ii)),'}$']); "$D_\beta$"], 'interpreter', 'latex', 'Fontsize', fontsize+2)
        xticks([1:5,7,8])
        xticklabels({'$1$','$2$','$3$', '$4$','$5$', '$2^\prime$'})
        set(gca, 'TickLabelInterpreter','latex', 'fontsize', fontsize)
        ylim([0.1e-3,0.8e-3])
        xlim([0.5,7.5])
    end

    nexttile
    for ij = 1:length(Dvec)
        load(['VaryDeltaSim_SNR_', num2str(SNR(ii)), '_D_', num2str(Dvec(ij)*1e4),'_beta_',num2str(betavec(ij)*100)], 'betaFit','beta')

        F = zeros(100, size(betaFit,2));
        betavals = zeros(100, size(betaFit,2));
        for ik = 1:size(Dfit,2)
            [~,edges] = histcounts(betaFit(:,ik)); % use histcounts to find optimal bandwidth for kernel density function
            bw = edges(2)-edges(1);

            [F(:,ik), betavals(:,ik)]=ksdensity(betaFit(:,ik),'bandwidth',bw);


            CVBeta(ik,ij+(ii-1)*length(Dvec)) = sqrt(var(betaFit(:,ik)))/mean(betaFit(:,ik))*100; % Compute coefficient of variation percentage
        end

        F = F/max(F(:))*0.45; % Scale the density distributions so they fit without overlapping.
        hold on
        for ik = 1:size(Dfit,2)
            fill([F(:,ik)+numDeltaAxisPos(ik); flipud(numDeltaAxisPos(ik)-F(:,ik))], [betavals(:,ik); flipud(betavals(:,ik))], colour(ij,:), 'FaceAlpha', alpha, 'edgecolor', colour(ij,:)) % Plot violin
            meanbeta = mean(betaFit(:,ik));
            meanF = interp1(betavals(:,ik), F(:,ik), meanbeta);
            plot([numDeltaAxisPos(ik)-meanF, numDeltaAxisPos(ik)+meanF],[meanbeta, meanbeta], '-', 'linewidth', 2, 'color',colour(ij,:))
        end

        plot([0,8],[beta,beta],':', 'color', colour(ij,:), 'linewidth', 2)
        xlabel('No. of diffusion times', 'interpreter', 'latex', 'fontsize', fontsize)
        ylabel("$\beta$", 'interpreter', 'latex', 'Fontsize', fontsize+2)
        xticks([1:5,7,8])
        xticklabels({'$1$','$2$','$3$', '$4$','$5$', '$2^\prime$'})
        set(gca, 'TickLabelInterpreter','latex', 'fontsize', fontsize)
        ylim([0.6,1])
        xlim([0.5,7.5])
    end

    nexttile
    for ij = 1:length(Dvec)
        load(['VaryDeltaSim_SNR_', num2str(SNR(ii)), '_D_', num2str(Dvec(ij)*1e4),'_beta_',num2str(betavec(ij)*100)], 'Kfit','K')
        F = zeros(100, size(betaFit,2));
        Kvals = zeros(100, size(betaFit,2));
        for ik = 1:size(Dfit,2)
            [~,edges] = histcounts(Kfit(:,ik)); % use histcounts to find optimal bandwidth for kernel density function
            bw = edges(2)-edges(1);

            [F(:,ik), Kvals(:,ik)]=ksdensity(Kfit(:,ik),'bandwidth',bw);


            CVK(ik,ij+(ii-1)*length(Dvec)) = sqrt(var(Kfit(:,ik)))/mean(Kfit(:,ik))*100; % Compute coefficient of variation percentage
        end

        F = F/max(F(:))*0.45; % Scale the density distributions so they fit without overlapping.
        hold on
        for ik = 1:size(Dfit,2)
            fill([F(:,ik)+numDeltaAxisPos(ik); flipud(numDeltaAxisPos(ik)-F(:,ik))], [Kvals(:,ik); flipud(Kvals(:,ik))], colour(ij,:), 'FaceAlpha', alpha, 'edgecolor', colour(ij,:)) % Plot violin
            meanK = mean(Kfit(:,ik));
            meanF = interp1(Kvals(:,ik), F(:,ik), meanK);
            plot([numDeltaAxisPos(ik)-meanF, numDeltaAxisPos(ik)+meanF],[meanK, meanK], '-', 'linewidth', 2, 'color',colour(ij,:))
        end
        plot([0,8],[K,K],':', 'color', colour(ij,:), 'linewidth', 2)
        xlabel('No. of diffusion times', 'interpreter', 'latex', 'fontsize', fontsize)
        ylabel("$K^*$", 'interpreter', 'latex', 'Fontsize', fontsize+2)
        xticks([1:5,7,8])
        xticklabels({'$1$','$2$','$3$', '$4$','$5$', '$2^\prime$'})
        set(gca, 'TickLabelInterpreter','latex', 'fontsize', fontsize)
        ylim([0.15,1.2])
        xlim([0.5,7.5])
    end
end

nexttile
plot(1:5, CVD(1:5,:), 'linewidth',2)
ylabel('CV of $D_\beta$ (\%)', 'interpreter', 'latex', 'fontsize', fontsize)
xlabel('No. of diffusion times', 'interpreter', 'latex', 'fontsize', fontsize)
% lgd = legend({['$\textrm{SNR} = ', num2str(SNR(1)),'$'],['$\textrm{SNR} = ', num2str(SNR(2)),'$'],['$\textrm{SNR} = ', num2str(SNR(3)),'$']},'Interpreter','latex', 'location', 'best');
set(gca, 'TickLabelInterpreter','latex', 'ColorOrder', lines(2), 'LineStyleOrder', {'-', '--',':'}, 'fontsize', fontsize)
xlim([1,5])


nexttile
plot(1:5, CVBeta(1:5,:),'linewidth',2)
ylabel('CV of $\beta$ (\%)', 'interpreter', 'latex', 'fontsize', fontsize)
xlabel('No. of diffusion times', 'interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'TickLabelInterpreter','latex', 'ColorOrder', lines(2), 'LineStyleOrder', {'-', '--',':'}, 'fontsize', fontsize)
xlim([1,5])

nexttile
h4 = plot(1:5, CVK(1:5,:), 'linewidth',2);
ylabel('CV of $K^*$ (\%)', 'interpreter', 'latex', 'fontsize', fontsize)
xlabel('No. of diffusion times', 'interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'TickLabelInterpreter','latex', 'ColorOrder', lines(2), 'LineStyleOrder', {'-', '--',':'}, 'fontsize', fontsize)
xlim([1,5])
lgd = legend({['WM: $\textrm{SNR} = ', num2str(SNR(1)),'\quad$'],['GM: $\textrm{SNR} = ', num2str(SNR(1)),'\quad$'],['WM: $\textrm{SNR} = ', num2str(SNR(2)),'\quad$'],['GM: $\textrm{SNR} = ', num2str(SNR(2)),'\quad$'],['WM: $\textrm{SNR} = ', num2str(SNR(3)),'$'],['GM: $\textrm{SNR} = ', num2str(SNR(3)),'$']},'Interpreter','latex', 'location', 'best', 'NumColumns',3, 'fontsize', 12);
lgd.Layout.Tile = 'South';

exportgraphics(gcf,'FitComparisonsSNR.eps')


%% Figure 2
% Make sure you run Simulation_VaryDelta_ParameterRanges first. This
% simulation takes a very long time.

h = fspecial('disk', 8);
figure
set(gcf, 'Position', [100,100, 600, 215])
t = tiledlayout(1,3, 'TileSpacing','none', 'padding', 'none');
nexttile
load('DeltaCombinations5_1000.mat', 'deltarange', 'diffRange', 'R_2_K')
imagesc(deltarange, diffRange, R_2_K)
clim([0.8,1])
hold on
contour(deltarange, diffRange, imfilter(R_2_K,h, 'replicate'), [0.85,0.9,0.95,0.99], 'linewidth', 2, 'color','k', 'showtext', 'on')
axis image
xlabel('$\Delta_1$ (ms)', 'interpreter', 'latex', 'fontsize', fontsize)
ylabel('$\Delta_2 - \Delta_1$ (ms)', 'interpreter', 'latex', 'fontsize', fontsize)
title('SNR = 5', 'interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'Ydir', 'normal', 'ticklabelinterpreter', 'latex', 'fontsize', fontsize)


nexttile
load('DeltaCombinations10_1000.mat', 'deltarange', 'diffRange', 'R_2_K')
imagesc(deltarange, diffRange, R_2_K)
clim([0.8,1])
axis image
hold on
contour(deltarange, diffRange, imfilter(R_2_K,h, 'replicate'), [0.85,0.9,0.95,0.99], 'linewidth', 2, 'color','k', 'showtext', 'on')
xlabel('$\Delta_1$ (ms)', 'interpreter', 'latex', 'fontsize', fontsize)
yticks([])
title('SNR = 10', 'interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'Ydir', 'normal', 'ticklabelinterpreter', 'latex', 'fontsize', fontsize)


nexttile
load('DeltaCombinations20_1000.mat', 'deltarange', 'diffRange', 'R_2_K')
imagesc(deltarange, diffRange, R_2_K)
clim([0.8,1])
axis image
hold on
contour(deltarange, diffRange, imfilter(R_2_K,h, 'replicate'), [0.85,0.9,0.95,0.99], 'linewidth', 2, 'color','k', 'showtext', 'on')
xlabel('$\Delta_1$ (ms)', 'interpreter', 'latex', 'fontsize', fontsize)
yticks([])
title('SNR = 20', 'interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'Ydir', 'normal', 'ticklabelinterpreter', 'latex', 'fontsize', fontsize)

c = colorbar;
c.Layout.Tile = 'East';
c.Label.String = '$R^2$';
c.Label.Interpreter = 'latex';
set(c, 'ticklabelinterpreter', 'latex', 'fontsize',fontsize)

exportgraphics(gcf, 'DeltaCombinationsComb.eps')


%% Figure 3
% Make sure you run Simulation_ParameterRanges_Fit first.
SNR = [5, 10, 20, inf];
fontsize = 9;
facealpha = 0.75;

h3 = figure;
set(gcf, 'Position',[100,100,540, 780])
t3 = tiledlayout(4,length(SNR),"TileSpacing","tight", "Padding","loose", "TileIndexing","columnmajor");

min_beta = 0.5;

R_2_String = strings(4, length(SNR));

for ii = 1:length(SNR)
    SNRString = string(['SNR = ', num2str(SNR(ii))]);

load(['BetaDfit_', num2str(SNR(ii)), '_', num2str(min_beta*10), '.mat'])

Dvec = Dvec/1e6;
DFitVec = DFitVec/1e6;
N = length(betaVec);

KFitVec = 6* (gamma(1+betaFitVec)).^2 ./ gamma(1+2*betaFitVec) - 3;
KVec = 6* (gamma(1+betaVec)).^2 ./ gamma(1+2*betaVec) - 3;

R_2_K_1 = 1-sum((KVec - KFitVec(:,1)).^2)/sum( (KVec - mean(KVec)).^2);%((N-1)*var(KFitVec(:,1)));
R_2_K_2 = 1-sum((KVec - KFitVec(:,2)).^2)/sum( (KVec - mean(KVec)).^2);%((N-1)*var(KFitVec(:,2)));
R_2_K_3 = 1-sum((KVec - KFitVec(:,3)).^2)/sum( (KVec - mean(KVec)).^2);%((N-1)*var(KFitVec(:,3)));
R_2_K_4 = 1-sum((KVec - KFitVec(:,4)).^2)/sum( (KVec - mean(KVec)).^2);%((N-1)*var(KFitVec(:,4)));
R_2_String(1,ii) = string(['$R^2 = ', num2str(R_2_K_1, '%.2f'),'$']);
R_2_String(2,ii) = string(['$R^2 = ', num2str(R_2_K_2, '%.2f'),'$']);
R_2_String(3,ii) = string(['$R^2 = ', num2str(R_2_K_3, '%.2f'),'$']);
R_2_String(4,ii) = string(['$R^2 = ', num2str(R_2_K_4, '%.2f'),'$']);

figure(h3)
nexttile
h4(1) = plot(KVec, KFitVec(:,1), '.');
hold on
axis equal
box on
h4(2) = plot(KVec,KFitVec_DKI19(:,1),  '.');
xlim([0,1.7])
ylim([0,3])
plot([0,3],[0,3], 'k-', 'linewidth',2)
xlabel('Simulated $K$','Interpreter','latex', 'fontsize', fontsize)
ylabel('Fitted $K$', 'Interpreter','latex', 'fontsize', fontsize)
title(SNRString, 'fontsize', fontsize+2)
if ii == 4
ax = gca;
    hCopy = copyobj(h4,ax);
    set(hCopy(1), 'XData', NaN', 'YData', NaN);
    set(hCopy(2), 'XData', NaN', 'YData', NaN);

    hCopy(1).MarkerSize = 20;
    hCopy(2).MarkerSize = 20;

l = legend(hCopy,{'Sub-diffusion','DKI'}, 'Numcolumns', 2, 'interpreter', 'latex', 'fontsize', fontsize);
l.Layout.Tile = "north";
end
set(gca, 'fontsize',fontsize, 'TickLabelInterpreter','latex')
nexttile
plot(KVec,KFitVec(:,2),  '.')
axis equal
hold on
box on
plot(KVec, KFitVec_DKI49(:,1), '.')
ylim([0,3])
xlim([0,1.7])
plot([0,3],[0,3], 'k-', 'linewidth',2)
xlabel('Simulated $K$','Interpreter','latex', 'fontsize', fontsize)
ylabel('Fitted $K$', 'Interpreter','latex', 'fontsize', fontsize)
set(gca, 'fontsize',fontsize, 'TickLabelInterpreter','latex')
nexttile
plot(KVec, KFitVec(:,3), '.')
axis equal
box on
hold on
ylim([0,3])
xlim([0,1.7])
plot([0,3],[0,3], 'k-', 'linewidth',2)
xlabel('Simulated $K$','Interpreter','latex', 'fontsize', fontsize)
ylabel('Fitted $K$', 'Interpreter','latex', 'fontsize', fontsize)
set(gca, 'fontsize',fontsize, 'TickLabelInterpreter','latex')
nexttile
plot(KVec, KFitVec(:,4), '.')
axis equal
hold on
box on
ylim([0,3])
xlim([0,1.7])
plot([0,3],[0,3], 'k-', 'linewidth',2)
xlabel('Simulated $K$','Interpreter','latex', 'fontsize', fontsize)
ylabel('Fitted $K$', 'Interpreter','latex', 'fontsize', fontsize)
set(gca, 'fontsize',fontsize, 'TickLabelInterpreter','latex')
end

nexttile(1)
text(-0.9,1.5, "$\mathbf{\Delta = 19}$ ms",'rotation',90, 'interpreter', 'latex', 'VerticalAlignment','middle', 'HorizontalAlignment','center', 'FontWeight','bold', 'fontsize', fontsize+2)
nexttile(2)
text(-0.9,1.5, "$\mathbf{\Delta = 49}$ ms",'rotation',90, 'interpreter', 'latex', 'VerticalAlignment','middle', 'HorizontalAlignment','center', 'FontWeight','bold', 'fontsize', fontsize+2)
nexttile(3)
text(-0.9,1.5, "$\mathbf{\Delta = 19,49}$ ms",'rotation',90, 'interpreter', 'latex', 'VerticalAlignment','middle', 'HorizontalAlignment','center', 'FontWeight','bold', 'fontsize', fontsize+2)
nexttile(4)
text(-0.9,1.5, "$\mathbf{\Delta = 19,34,49}$ ms",'rotation',90, 'interpreter', 'latex', 'VerticalAlignment','middle', 'HorizontalAlignment','center', 'FontWeight','bold', 'fontsize', fontsize+2)

for ii = 1:4*length(SNR)
    nexttile(ii)
    pos = get(gca, 'Position');
    annotation('textbox', [[0.4,0.8].*pos(3:4)+pos(1:2),0.03,0.03],'String', R_2_String(ii), 'FitBoxToText', 'on', 'fontsize', fontsize, 'margin', 2, 'VerticalAlignment','bottom', 'backgroundcolor', 'w','interpreter', 'latex', 'HorizontalAlignment','center', 'EdgeColor','none', 'FaceAlpha',facealpha)

end

exportgraphics(h3, 'SimulatedKComparisons_2.pdf','ContentType','vector')

%% Figure 4
% Run Simulation_OptimalBValues first
SNR = [5,10,20];
numB = [3,4,5];
fontsize = 10;
figure
set(gcf, 'position', [490,310,600,600])
tiledlayout(length(SNR), length(numB),'tilespacing', 'tight', 'padding', 'loose')
above = false(length(SNR), length(numB));
posMult = cell(length(SNR), length(numB));
paru5 = parula(5);
paru4 = parula(4);
paruNew = [paru5(1:2,:); paru4(2,:); paru5(3,:); paru4(3,:); paru5(4:5,:)];

bvalues = [0,50,350,800,1500,2400,3450,4750,6000,200,950,2300,4250,6750,9850,13500,17800];

% Define properties for each of the zoomed in boxes
zoomNum = [4, 7, 8; 5, 11, 10; 6, 10, 10]-1;
lc_x = [0.25,0.25,0.25;0.25,0.25,0.4;0.25,0.55,0.3];
lc_y = [0.5,0.5,0.5;0.4,0.7,0.1;0.5,0.1,0.1];
rw = [0.6,0.6,0.45; 0.5,0.45,0.55; 0.4,0.4,0.65];
rh = [0.4,0.4,0.4;0.5,0.2,0.3;0.4,0.4,0.4];


for ii = 1:length(SNR)
    for ij = 1:length(numB)
        load(['b_red_data_', num2str(numB(ij)), '_', num2str(SNR(ii))])
        type_ii_ij = sum(combsIdx>9,2)/(numB(ij)-1);
        [R_2_K_plot, order_ii_ij] = sort(R_2_K);
        xmax = length(R_2_K);
        xmin = xmax - zoomNum(ii,ij);
        ymax = ceil(R_2_K_plot(xmax)*100)/100;
        ymin = floor(R_2_K_plot(xmin)*100)/100;
        nexttile
        scatter(1:xmax, R_2_K_plot, 20, type_ii_ij(order_ii_ij), 'filled')
        xlim([0,xmax])
        ylim([0,1])
        ylabel('$R^2$', 'interpreter', 'latex', 'fontsize', fontsize)
        colormap(paruNew)
        set(gca, 'fontsize',fontsize,'ticklabelinterpreter','latex')
        if ii == 1
            title([num2str(numB(ij)-1), ' non-zero b-values'], 'fontsize', fontsize, 'interpreter', 'latex')
        end
        if ij == 1
            text(-0.5*xmax, 0.5, ['\textbf{SNR = ', num2str(SNR(ii)),'}'], 'interpreter', 'latex', 'rotation',90, 'fontsize', fontsize, 'verticalalignment', 'middle', 'horizontalalignment', 'center')
        end

        % Add overlay
        hold on
        rectangle('Position',[lc_x(ii,ij)*xmax, lc_y(ii,ij), rw(ii,ij)*xmax, rh(ii,ij)])
        text(lc_x(ii,ij)*xmax, lc_y(ii,ij)-0.01, num2str(xmin), 'horizontalalignment', 'center', 'verticalalignment', 'top','interpreter', 'latex', 'fontsize', fontsize)
        text((lc_x(ii,ij)+rw(ii,ij))*xmax, lc_y(ii,ij)-0.01, num2str(xmax), 'horizontalalignment', 'center', 'verticalalignment', 'top','interpreter', 'latex', 'fontsize', fontsize)
        text((lc_x(ii,ij)-0.01)*xmax, lc_y(ii,ij), num2str(ymin), 'horizontalalignment', 'right', 'verticalalignment', 'middle','interpreter', 'latex', 'fontsize', fontsize)
        text((lc_x(ii,ij)-0.01)*xmax, lc_y(ii,ij)+rh(ii,ij), num2str(ymax), 'horizontalalignment', 'right', 'verticalalignment', 'middle','interpreter', 'latex', 'fontsize', fontsize)
        scatter(lc_x(ii,ij)*xmax+(0:zoomNum(ii,ij))*xmax*rw(ii,ij)/zoomNum(ii,ij), (R_2_K_plot(xmin:xmax)-ymin)/(ymax-ymin)*rh(ii,ij)+lc_y(ii,ij), 20,type_ii_ij(order_ii_ij(xmin:xmax)), 'filled')
        if ymin<lc_y(ii,ij) 
            plot([xmin,lc_x(ii,ij)*xmax], [ymin, lc_y(ii,ij)], 'k-')
        else
            plot([xmax,(lc_x(ii,ij)+rw(ii,ij))*xmax], [ymin, lc_y(ii,ij)], 'k-')
        end
        if ymax<lc_y(ii,ij)+rh(ii,ij)
            plot([xmax, (lc_x(ii,ij)+rw(ii,ij))*xmax], [ymax, lc_y(ii,ij)+rh(ii,ij)], 'k-')
        else
            plot([xmin, (lc_x(ii,ij))*xmax], [ymax, lc_y(ii,ij)+rh(ii,ij)], 'k-')
        end

        % Display the best 5 combinations of b-values
        disp([num2str(numB(ij)-1), ' b-values '])
        bident = char(repmat([37,100,44,32], 1, numB(ij)-1));
        fprintf(['b= ',bident,' Rsquared = %.2f \n'], [bvalues(combsIdx(order_ii_ij(end-4:end),:)),R_2_K_plot(end-4:end)]')

        
    end
end
c = colorbar('ticklabelinterpreter', 'latex', 'ticks', [1/14, 3/14, 5/14, 7/14, 9/14, 11/14, 13/14], 'TickLabels', {'1:0', '3:1', '2:1', '1:1', '1:2', '1:3', '0:1'}, 'fontsize', fontsize);
c.Label.String = 'Proportion of $b$ values from $\Delta =19$ to $\Delta = 49$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = fontsize;
c.Layout.Tile = 'East';


%%
