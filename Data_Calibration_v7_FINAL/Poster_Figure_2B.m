%%=========================================================================
clear; clc; close all; 
format long g;

%% Obtaining Fluorescence Data in Concentration Space
cd 'Concentration_Data_Avg'

list = dir;

fileNames = [];
for i = 3:length(list) %ignores . and ..
    fileNames = [fileNames,convertCharsToStrings(list(i).name)];
end

sort(fileNames);

filesToPlot = [];

for i = 1:(length(fileNames)) 
    currFileName = fileNames(i);

    if (contains(currFileName,'_LR_') == 0 && ...
            contains(currFileName,'_New_') == 0)
        filesToPlot = [filesToPlot,currFileName];
    end

end

%% ========================================================================

% Plot control

figure(1)
set(gcf,'Units','inches','Position',[0.5 0.5 10 9.85]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 10 9.85],'color','white');

currData = load(filesToPlot(2)); %Ctrl blood
subplot(3,2,1)
%plot ctrl blood on IxR60
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-^g',...
        'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g','CapSize',4); hold on;

subplot(3,2,2)
%plot ctrl blood on IxR4H
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-^g',...
        'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g','CapSize',4); hold on;

currData = load(filesToPlot(3)); %Ctrl plasma
subplot(3,2,3)
%plot ctrl plasma on IxR60
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-^g',...
        'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g','CapSize',4); hold on;

subplot(3,2,4)
%plot ctrl plasma on IxR4H
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-^g',...
        'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g','CapSize',4); hold on;

currData = load(filesToPlot(1)); %Ctrl plasma
subplot(3,2,5)
%plot ctrl bile on IxR60
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-^g',...
        'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g','CapSize',4); hold on;

subplot(3,2,6)
%plot ctrl bile on IxR4H
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-^g',...
        'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g','CapSize',4); hold on;
%% ========================================================================

% Plot IRI

for i = 1:(length(filesToPlot)/3 - 1)
    currBile = load(filesToPlot(3*i+1));
    currBlood = load(filesToPlot(3*i+2));
    currPlasma = load(filesToPlot(3*i+3));
    
    if i > 2
        colorStr = "-or";
        colorChar = 'r';
        % I60Rxx
    else
        colorStr = "-sb";
        colorChar = 'b';
        % I30Rxx
    end

    if mod(i,2) == 0
        %IxxR60
        subplot(3,2,1)
%         errorbar(currBlood(:,1),currBlood(:,2)*10^2,currBlood(:,4)*10^2,colorStr,...
%         'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',colorChar,...
%         'MarkerEdgeColor',colorChar,'CapSize',4); hold on;
        axis([0 60 0 1.5]);
        set(gca,'XTick',(0:10:60),'YTick',(0:0.3:1.5));
        set(gca,'LineWidth',1.5,'FontSize',18); box off;
        ytickformat('%.1f')
        xlabel('Time (min)'); ylabel(sprintf('Concentration/10^2\n(mg/mL)'));
        title(sprintf('Whole Blood'));
        legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','NorthEast',...
        'FontSize',15);
        legend('boxoff')

        subplot(3,2,3)
%         errorbar(currPlasma(:,1),currPlasma(:,2)*10^2,currPlasma(:,4)*10^2,colorStr,...
%         'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',colorChar,...
%         'MarkerEdgeColor',colorChar,'CapSize',4); hold on;
        axis([0 60 0 2.5]);
        set(gca,'XTick',(0:10:60),'YTick',(0:0.5:2.5));
        set(gca,'LineWidth',1.5,'FontSize',18); box off;
        ytickformat('%.1f')
        xlabel('Time (min)'); ylabel(sprintf('Concentration/10^2\n(mg/mL)'));
        title(sprintf('Plasma'));

        subplot(3,2,5)
%         errorbar(currBile(:,1),currBile(:,2)*10^2,currBile(:,4)*10^2,colorStr,...
%         'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',colorChar,...
%         'MarkerEdgeColor',colorChar,'CapSize',4); hold on;
        axis([0 60 0 5]);
        set(gca,'XTick',(0:10:60),'YTick',(0:1:5));
        set(gca,'LineWidth',1.5,'FontSize',18); box off;
        ytickformat('%.1f')
        xlabel('Time (min)'); ylabel(sprintf('Concentration/10^2\n(mg/mL)'));
        title(sprintf('Bile'));
    
    else
        %IxxR4H
        subplot(3,2,2)
%         errorbar(currBlood(:,1),currBlood(:,2)*10^2,currBlood(:,4)*10^2,colorStr,...
%         'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',colorChar,...
%         'MarkerEdgeColor',colorChar,'CapSize',4); hold on;
        axis([0 60 0 1.5]);
        set(gca,'XTick',(0:10:60),'YTick',(0:0.3:1.5));
        set(gca,'LineWidth',1.5,'FontSize',18); box off;
        h = gca;
        h.YAxis.Visible = 'off';
        xlabel('Time (min)'); %ylabel(sprintf('Concentration\n(mg/mL)'));
        title(sprintf('Whole Blood'));

        subplot(3,2,4)
%         errorbar(currPlasma(:,1),currPlasma(:,2)*10^2,currPlasma(:,4)*10^2,colorStr,...
%         'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',colorChar,...
%         'MarkerEdgeColor',colorChar,'CapSize',4); hold on;
        axis([0 60 0 2.5]);
        set(gca,'XTick',(0:10:60),'YTick',(0:0.5:2.5));
        set(gca,'LineWidth',1.5,'FontSize',18); box off;
        h = gca;
        h.YAxis.Visible = 'off';
        xlabel('Time (min)'); %ylabel(sprintf('Concentration\n(mg/mL)'));
        title(sprintf('Plasma'));

        subplot(3,2,6)
%         errorbar(currBile(:,1),currBile(:,2)*10^2,currBile(:,4)*10^2,colorStr,...
%         'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',colorChar,...
%         'MarkerEdgeColor',colorChar,'CapSize',4); hold on;
        axis([0 60 0 5.0]);
        set(gca,'XTick',(0:10:60),'YTick',(0:1.0:5.0));
        set(gca,'LineWidth',1.5,'FontSize',18); box off;
        h = gca;
        h.YAxis.Visible = 'off';
        xlabel('Time (min)'); %ylabel(sprintf('Concentration\n(mg/mL)'));
        title(sprintf('Bile'));
    end

end

%% ========================================================================

