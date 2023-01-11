%%=========================================================================
clear; clc; close all; 

%% Load Full Calibration Dataset
cd Calibration_Data_Indiv;
BloodData = load('Blood_CC_Full_Indiv.txt');
PlasmaData = load('Plasma_CC_Full_Indiv.txt');
BileData = load('Bile_CC_Full_Indiv.txt');
cd ..;

%% Parse Full Calibration Dataset
Bl_C_Data = BloodData(:,1);
Bl_F_Data = BloodData(:,2:end);
Bl_F_Avg = mean(Bl_F_Data,2,'omitnan');
Bl_F_SD = std(Bl_F_Data,0,2,'omitnan');

for i = 1:length(Bl_F_Data(:,1))
    Bl_F_SEM (i,:) = Bl_F_SD(i)/sqrt(length(Bl_F_Data(~isnan(Bl_F_Data(i,:)))));
end

Pl_C_Data = PlasmaData(:,1);
Pl_F_Data = PlasmaData(:,2:end);
Pl_F_Avg = mean(Pl_F_Data,2,'omitnan');
Pl_F_SD = std(Pl_F_Data,0,2,'omitnan');

for i = 1:length(Pl_F_Data(:,1))
    Pl_F_SEM (i,:) = Pl_F_SD(i)/sqrt(length(Pl_F_Data(~isnan(Pl_F_Data(i,:)))));
end

Bi_C_Data = BileData(:,1);
Bi_F_Data = BileData(:,2:end);
Bi_F_Avg = mean(Bi_F_Data,2,'omitnan');
Bi_F_SD = std(Bi_F_Data,0,2,'omitnan');

for i = 1:length(Bi_F_Data(:,1))
    Bi_F_SEM (i,:) = Bi_F_SD(i)/sqrt(length(Bi_F_Data(~isnan(Bi_F_Data(i,:)))));
end

% Save mean, SD, SEM data
if(exist('Calibration_Data_Avg','file') == 0)
    mkdir('Calibration_Data_Avg') 
end

cd Calibration_Data_Avg;

Bl_Avg_Full = [Bl_C_Data,Bl_F_Avg,Bl_F_SD,Bl_F_SEM];
Pl_Avg_Full = [Pl_C_Data,Pl_F_Avg,Pl_F_SD,Pl_F_SEM];
Bi_Avg_Full = [Bi_C_Data,Bi_F_Avg,Bi_F_SD,Bi_F_SEM];

% writematrix(Bl_Avg_Full,"Blood_CC_Full_Avg_SD_SEM");
% writematrix(Pl_Avg_Full,"Plasma_CC_Full_Avg_SD_SEM");
% writematrix(Bi_Avg_Full,"Bile_CC_Full_Avg_SD_SEM");

cd ..;

%% Load Calibration Dataset for Experimental Data Ranges
cd Calibration_Data_Indiv;
BloodData = load('Blood_CC_Exp Range_Indiv.txt');
PlasmaData = load('Plasma_CC_Exp Range_Indiv.txt');
BileData = load('Bile_CC_Exp Range_Indiv.txt');
cd ..;

%% Parse Experimental Range Calibration Dataset
Bl_C_Data_ER = BloodData(:,1);
Bl_F_Data_ER = BloodData(:,2:end);
Bl_F_Avg_ER = mean(Bl_F_Data_ER,2,'omitnan');
Bl_F_SD_ER = std(Bl_F_Data_ER,0,2,'omitnan');

for i = 1:length(Bl_F_Data_ER(:,1))
    Bl_F_SEM_ER (i,:) = Bl_F_SD_ER(i)/sqrt(length(Bl_F_Data_ER(~isnan(Bl_F_Data_ER(i,:)))));
end

Pl_C_Data_ER = PlasmaData(:,1);
Pl_F_Data_ER = PlasmaData(:,2:end);
Pl_F_Avg_ER = mean(Pl_F_Data_ER,2,'omitnan');
Pl_F_SD_ER = std(Pl_F_Data_ER,0,2,'omitnan');

for i = 1:length(Pl_F_Data_ER(:,1))
    Pl_F_SEM_ER (i,:) = Pl_F_SD_ER(i)/sqrt(length(Pl_F_Data_ER(~isnan(Pl_F_Data_ER(i,:)))));
end

Bi_C_Data_ER = BileData(:,1);
Bi_F_Data_ER = BileData(:,2:end);
Bi_F_Avg_ER = mean(Bi_F_Data_ER,2,'omitnan');
Bi_F_SD_ER = std(Bi_F_Data_ER,0,2,'omitnan');

for i = 1:length(Bi_F_Data_ER(:,1))
    Bi_F_SEM_ER (i,:) = Bi_F_SD_ER(i)/sqrt(length(Bi_F_Data_ER(~isnan(Bi_F_Data_ER(i,:)))));
end

cd Calibration_Data_Avg;

Bl_Avg_ExpRng = [Bl_C_Data_ER,Bl_F_Avg_ER,Bl_F_SD_ER,Bl_F_SEM_ER];
Pl_Avg_ExpRng = [Pl_C_Data_ER,Pl_F_Avg_ER,Pl_F_SD_ER,Pl_F_SEM_ER];
Bi_Avg_ExpRng = [Bi_C_Data_ER,Bi_F_Avg_ER,Bi_F_SD_ER,Bi_F_SEM_ER];

% writematrix(Bl_Avg_ExpRng,"Blood_CC_Exp Range_Avg_SD_SEM");
% writematrix(Pl_Avg_ExpRng,"Plasma_CC_Exp Range_Avg_SD_SEM");
% writematrix(Bi_Avg_ExpRng,"Bile_CC_Exp Range_Avg_SD_SEM");

cd ..;

%If this is the first time running the script, you will need to stop the
%script here and run the MATLAB curvefitting tool to get the nonlinear and
%linear calibration parameters. All of the variables will be loaded into
%the current workspace and you just need to select the appropriate one from
%the dropdown menu");

%% Nonlinear Calibration Curves for Full Calibration Dataset
% Nonlinear Model Parameters Using Full Calibration Data Ranges
pars_Bl = [-145.8,7.769E5];
pars_Pl = [3.992E4,0.02224];
pars_Bi = [2.636E4,0.04628];
cSpan = (0:0.001:0.05);

M_Bl = linearModel(cSpan,pars_Bl);
M_Pl = nonLinearModel(cSpan,pars_Pl);
M_Bi = nonLinearModel(cSpan,pars_Bi);

% Linear Model Parameters Using Experimental Data Range
pars_Bl = [-106.8,7.841E5];
pars_Pl = [4800,5.85E5];
pars_Bi = [1778,2.623E5];

M_Bl_L = linearModel(cSpan,pars_Bl);
M_Pl_L = linearModel(cSpan,pars_Pl);
M_Bi_L = linearModel(cSpan,pars_Bi);

%% Plotting Nonlinear Calibration Curves for Full Calibration Data Range
figure(1);
set(gcf,'Units','inches','Position',[0.5 0.5 18 3.75]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 18 3.75],'color','white');

subplot(1,3,1)
errorbar(Bl_C_Data,Bl_F_Avg/10^4,Bl_F_SEM/10^4,'o','MarkerSize',6,'CapSize',4,...
    'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k'); hold on;
plot(cSpan,M_Bl/10^4,'-k','LineWidth',2.0); hold on;
plot(cSpan,M_Bl_L/10^4,'--k','LineWidth',2.0); hold on;
set(gca,'LineWidth',1.25,'FontSize',16); box off;
axis([0 0.05 0 5]);
set(gca,'XTick',(0:0.01:0.05),'YTick',(0:1:5));
set(gca,'LineWidth',1.5,'FontSize',18); box off;
ytickformat('%.1f')
xlabel('Concentration (mg/mL)'); ylabel(sprintf('Fluorescence x 10^4\n (A.U.)'));
title('Blood Calibration');
h=legend("Calibration Data","Nonlinear Regression","Linear Regression","Location","NorthWest");
set(h, "fontSize",15); legend("boxoff")

subplot(1,3,2)
errorbar(Pl_C_Data,Pl_F_Avg/10^4,Pl_F_SEM/10^4,'o','MarkerSize',6,'CapSize',4,...
    'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k'); hold on;
plot(cSpan,M_Pl/10^4,'-k','LineWidth',2.0); hold on;
plot(cSpan,M_Pl_L/10^4,'--k','LineWidth',2.0); hold on;
set(gca,'LineWidth',1.25,'FontSize',16); box off;
axis([0 0.05 0 3.5]);
set(gca,'XTick',(0:0.01:0.05),'YTick',(0:0.7:3.5));
ytickformat('%.1f')
set(gca,'LineWidth',1.5,'FontSize',18); box off;
xlabel('Concentration (mg/mL)'); ylabel(sprintf('Fluorescence x 10^4\n (A.U.)'));
title('Plasma Calibration');

subplot(1,3,3)
errorbar(Bi_C_Data,Bi_F_Avg/10^4,Bi_F_SEM/10^4,'o','MarkerSize',6,'CapSize',4,...
    'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k'); hold on;
plot(cSpan,M_Bi/10^4,'-k','LineWidth',2.0); hold on;
plot(cSpan,M_Bi_L/10^4,'--k','LineWidth',2.0); hold on;
set(gca,'LineWidth',1.25,'FontSize',16); box off;
axis([0 0.05 0 1.5]);
set(gca,'XTick',(0:0.01:0.05),'YTick',(0:0.3:1.5));
ytickformat('%.1f')
set(gca,'LineWidth',1.5,'FontSize',18); box off;
xlabel('Concentration (mg/mL)'); ylabel(sprintf('Fluorescence x 10^4\n (A.U.)'));
title('Bile Calibration');

%% Linear Calibration Curves for Calibration Data in Range of Experimental Data
% Nonlinear Model Parameters Using Full Calibration Data Ranges
pars_Bl = [-145.8,7.769E5];
pars_Pl = [3.992E4,0.02224];
pars_Bi = [2.636E4,0.04628];
cSpan = (0:0.001:0.05);

M_Bl = linearModel(cSpan,pars_Bl);
M_Pl = nonLinearModel(cSpan,pars_Pl);
M_Bi = nonLinearModel(cSpan,pars_Bi);

% Linear Model Parameters Using Experimental Data Range
pars_Bl = [-106.8,7.841E5];
pars_Pl = [4800,5.85E5];
pars_Bi = [1778,2.623E5];

M_Bl_L = linearModel(cSpan,pars_Bl);
M_Pl_L = linearModel(cSpan,pars_Pl);
M_Bi_L = linearModel(cSpan,pars_Bi);

%% Plotting Linear Calibration Curves for Calibration Data in Range of Experimental Data
% subplot(3,2,2)
% errorbar(Bl_C_Data_ER,Bl_F_Avg_ER/10^4,Bl_F_SEM_ER/10^4,'o','MarkerSize',6,'CapSize',4,...
%     'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k'); hold on;
% plot(cSpan,M_Bl/10^4,'-k','LineWidth',2.0); hold on;
% plot(cSpan,M_Bl_L/10^4,'--k','LineWidth',2.0); hold on;
% axis([0 0.025 0 2.5]);
% set(gca,'XTick',(0:0.005:0.025),'YTick',(0:0.5:2.5));
% set(gca,'LineWidth',1.25,'FontSize',16); box off;
% xlabel('Concentration (mg/mL)'); ylabel('Fluorescence (A.U.)/10^4');
% title(sprintf('Blood Calibration\n(In Range of Experimental Data)\n'));
% 
% subplot(3,2,4)
% errorbar(Pl_C_Data_ER,Pl_F_Avg_ER/10^4,Pl_F_SEM_ER/10^4,'o','MarkerSize',6,'CapSize',4,...
%     'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k'); hold on;
% plot(cSpan,M_Pl/10^4,'-k','LineWidth',2.0); hold on;
% plot(cSpan,M_Pl_L/10^4,'--k','LineWidth',2.0); hold on;
% axis([0 0.04 0 3]);
% set(gca,'XTick',(0:0.008:0.04),'YTick',(0:0.6:3));
% set(gca,'LineWidth',1.25,'FontSize',16); box off;
% xlabel('Concentration (mg/mL)'); ylabel('Fluorescence (A.U.)/10^4');
% title(sprintf('Plasma Calibration\n(In Range of Experimental Data)\n'));
% 
% subplot(3,2,6)
% errorbar(Bi_C_Data_ER,Bi_F_Avg_ER/10^4,Bi_F_SEM_ER/10^4,'o','MarkerSize',6,'CapSize',4,...
%     'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k'); hold on;
% plot(cSpan,M_Bi/10^4,'-k','LineWidth',2.0); hold on;
% plot(cSpan,M_Bi_L/10^4,'--k','LineWidth',2.0); hold on;
% axis([0 0.05 0 1.5]);
% set(gca,'XTick',(0:0.01:0.05),'YTick',(0:0.3:1.5));
% set(gca,'LineWidth',1.25,'FontSize',16); box off;
% xlabel('Concentration (mg/mL)'); ylabel('Fluorescence (A.U.)/10^4');
% title(sprintf('Bile Calibration\n(In Range of Experimental Data)\n'));

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