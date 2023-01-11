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

figure(1); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

currData = load(filesToPlot(2)); %Ctrl blood
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-or',...
        'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor','r',...
        'MarkerEdgeColor','r','CapSize',6); hold on;


currData = load(filesToPlot(3)); %Ctrl plasma
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-sb',...
        'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor','b',...
        'MarkerEdgeColor','b','CapSize',6); hold on;


currData = load(filesToPlot(1)); %Ctrl bile
errorbar(currData(:,1),currData(:,2)*10^2,currData(:,4)*10^2,'-^g',...
        'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g','CapSize',6); hold on;

axis([0 60 0 4]);
set(gca,'XTick',(0:10:60),'YTick',(0:1:4));
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
ytickformat('%.1f')
xlabel('Time (min)'); ylabel(sprintf('Concentration x 10^2 (mg/mL)'));
legend('Blood data','Plasma data','Bile data','Location','North',...
'FontSize',16);
legend('boxoff')