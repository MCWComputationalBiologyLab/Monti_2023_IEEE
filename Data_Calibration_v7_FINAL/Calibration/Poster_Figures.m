%%=========================================================================
clear; clc; close all; 

%% Load Full Calibration Dataset
cd Calibration_Data_Avg;
BloodData = load('Blood_CC_Full_Avg_SD_SEM.txt');
PlasmaData = load('Plasma_CC_Full_Avg_SD_SEM.txt');
BileData = load('Bile_CC_Full_Avg_SD_SEM.txt');
cd ..;

pars_Bl = [-145.8,7.769E5];
pars_Pl = [3.992E4,0.02224];
pars_Bi = [2.636E4,0.04628];
cSpan = (0:0.001:0.05);

M_Bl = linearModel(cSpan,pars_Bl);
M_Pl = nonLinearModel(cSpan,pars_Pl);
M_Bi = nonLinearModel(cSpan,pars_Bi);

figure(1); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

errorbar(BloodData(:,1),BloodData(:,2)/10^4,BloodData(:,4)/10^4,'or','MarkerSize',8,'CapSize',6,...
    'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r','LineWidth',1.5); hold on; 
errorbar(PlasmaData(:,1),PlasmaData(:,2)/10^4,PlasmaData(:,4)/10^4,'sb','MarkerSize',8,'CapSize',6,...
    'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b','LineWidth',1.5); hold on;
errorbar(BileData(:,1),BileData(:,2)/10^4,BileData(:,4)/10^4,'^g','MarkerSize',8,'CapSize',6,...
    'MarkerFaceColor','g','MarkerEdgeColor','g','Color','g','LineWidth',1.5'); hold on;

plot(cSpan,M_Bl/10^4,'-r','LineWidth',2.5); hold on;
plot(cSpan,M_Pl/10^4,'-b','LineWidth',2.5); hold on;
plot(cSpan,M_Bi/10^4,'-g','LineWidth',2.5); hold on;

axis([0 0.05 0 5]); box off; grid off
set(gca,'XTick',(0:0.01:0.05),'YTick',(0:1:5)); box off
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman');
ytickformat('%.1f')
xlabel('Concentration (mg/mL)'); ylabel(sprintf('Fluorescence x 10^4 (A.U.)'));
legend('Blood calibration data','Plasma calibration data','Bile calibration data','Blood calibration curve','Plasma calibration curve',...
    'Bile calibration curve','Location','North','FontSize',16);
legend('boxoff')

%% ========================================================================
function F = linearModel(C,pars)
    B = pars(1);
    m = pars(2);
    F = B+m*C;
end

function F = nonLinearModel(C,pars)
    Vmax = pars(1);
    Km = pars(2);
    F = Vmax*C./(Km+C);
end