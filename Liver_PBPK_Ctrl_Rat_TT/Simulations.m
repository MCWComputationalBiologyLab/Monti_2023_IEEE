%%========================================================================
%close all; clear all; clc;

%%========================================================================
%% Load calibrated experimental data
C1 = load('Control_Pub_Blood_C_Data_NLR_Avg_SD_SEM.txt');
Blood_SEM = C1(:,4);
Blood_SD = C1(:,3);
BloodData = C1(:,2);
tdata = C1(:,1);

C2 = load('Control_Pub_Plasma_C_Data_NLR_Avg_SD_SEM.txt');
Plasma_SEM = C2(:,4);
Plasma_SD = C2(:,3);
PlasmaData = C2(:,2);

C3 = load('Control_Pub_Bile_C_Data_NLR_Avg_SD_SEM.txt');
Bile_SEM = C3(:,4);
Bile_SD = C3(:,3);
BileData = C3(:,2);

%%========================================================================
%% Simulate model
load('Ctrl_mpar.mat');
%mpar(2) = mpar(2); %change Vmax2

% Physiological parameters
p = Parameters(mpar);
Vrbc = p.Vrbc;
Vpl = p.Vpl;
lam_max = p.lamda_max;
lam_Km = p.lamda_Km;
lam_nH = p.lamda_nH;

% Simulation step size (dt) and maximum simulation time (tmax):
n = 18;
dt = 0.01; 
tmax = n*60;
t = (0:dt:tmax)';

% Calculate C in blood (input), C in RBCs, and C in plasma
Cin = zeros(length(t),1);
for i = 1:length(t)
    Cin(i) = InputGenerator(t(i));
end

lamda = 0.5 + lam_max*Cin.^lam_nH./(lam_Km^lam_nH+Cin.^lam_nH); % partition coefficient
Crbc_In = ((Vrbc+Vpl)./(Vrbc+lamda.*Vpl)).*Cin; %concentration in rbc based on concentration in blood
Cpl_In = ((Vrbc+Vpl)./(Vrbc./lamda+Vpl)).*Cin; %concentration in plasma based on concentration in blood

%Calculate plasma concentration of SFG
One_miunus_alpha = InputGenerator_SFC(t); % ratio of SF_F to SF_C
denom = 1./(One_miunus_alpha)-1; %for reading clarity
Cin_C_pl = Cpl_In./denom;
%Calculate blood concentration of SFG
Cin_C = (Vrbc./lamda+Vpl).*Cin_C_pl/(Vrbc+Vpl);
Cin_C_rbc = ((Vrbc+Vpl)./(Vrbc+lamda.*Vpl)).*Cin_C;

% Options for ODE15s solver
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:3), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Solve ODEs using ode15s
Y0 = [0,0,0,0,0,0]; % initial concentrations in mg/mL
ODE_FH = @(t,y)(LiverModel(t,y,p));
[T,Y] = ode15s(ODE_FH,t,Y0,options);
C1M = Y(:,1);   %conc of dye in blood (venous)
C2M = Y(:,2);   %conc of conjugated dye in blood (venous)
C3M = Y(:,3);   %conc of dye in hepatocytes
C4M = Y(:,4);   %conc of dye conjugate in hepatocytes
C5M = Y(:,5);   %conc of dye in bile
C6M = Y(:,6);   %conc of dye conjugate in bile

lamda = 0.5 + lam_max*C1M.^lam_nH./(lam_Km^lam_nH+C1M.^lam_nH); % partition coefficient
CplM = ((Vrbc+Vpl)./(Vrbc./lamda+Vpl)).*C1M; %concentration in plasma based on concentration in blood
CrbcM = ((Vrbc+Vpl)./(Vrbc+lamda*Vpl)).*C1M; %concentration in rbc based on concentration in blood
CplM_C = ((Vrbc+Vpl)./(Vrbc./lamda+Vpl)).*C2M; %conjugated concentration in plasma based on concentration in blood
CrbcM_C = ((Vrbc+Vpl)./(Vrbc+lamda*Vpl)).*C2M; %conjugated concentration in rbc based on concentration in blood

%%========================================================================
%% Plot model simulation results along with experimental data
figure(1); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(t,Cin*100,'-r','LineWidth',2.5); hold on;
plot(t,Cpl_In*100,'-b','LineWidth',2.5); hold on;
plot(t,Crbc_In*100,'-c','LineWidth',2.5); hold on;

plot(t,Cin_C*100,':r','LineWidth',2.5); hold on;
plot(t,Cin_C_pl*100,':b','LineWidth',2.5); hold on;
plot(t,Cin_C_rbc*100,':c','LineWidth',2.5); hold on;

errorbar(tdata,BloodData*100,Blood_SEM*100,'o','MarkerSize',6,'MarkerFaceColor','r',...
    'MarkerEdgeColor','r','LineWidth',2,'Color','r','CapSize',6); hold on;
errorbar(tdata,PlasmaData*100,Plasma_SEM*100,'o','MarkerSize',6,'MarkerFaceColor','b',...
    'MarkerEdgeColor','b','LineWidth',2,'Color','b','CapSize',6); hold on;

% Formatting
axis([0 n*60 0 4.5]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',18); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:0.9:4.5));
xlabel("Time (min)"); ylabel("Concentration x 10^2 (mg/mL)")
h=legend("Blood Inflow Model (SF)","Plasma Inflow Model (SF)","RBC Inflow Model (SF)",...
    "Blood Inflow Model (SFG)","Plasma Inflow Model (SFG)","RBC Inflow Model (SFG)",...
    "Blood Inflow Data","Plasma Inflow Data","Location","NorthWest");
set(h, "fontSize",15); legend("boxoff")

%% Plot model and data
figure(2); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(T, C1M*100, '-r', 'LineWidth', 2.5); hold on;
% SF Outflow in Blood
plot(T, C2M*100, ':r', 'LineWidth', 2.5); hold on;
% SF_C Outflow Blood
plot(T, CplM*100, '-b', 'LineWidth', 2.5); hold on;
% SF Outflow in Plasma
plot(T, CplM_C*100, ':b', 'LineWidth', 2.5); hold on;
% SF_C Outflow Plasma
plot(T, CrbcM*100, '-c', 'LineWidth', 2.5); hold on;
% SF Outflow in RBC
plot(T, CrbcM_C*100, ':c', 'LineWidth', 2.5); hold on;
% SF_C Outflow Plasma
plot(T, C3M*100, '-m', 'LineWidth', 2.5); hold on;
% SF in Hepatocytes
plot(T, C4M*100, ':m', 'LineWidth', 2.5); hold on;
% SF_C in Hepatocytes
plot(T, C5M*100, '-g', 'LineWidth', 2.5); hold on;
%SF in Bile
plot(T, C6M*100, ':g', 'LineWidth', 2.5); hold on;
%SF_C in Bile

errorbar(tdata,BileData*100,Bile_SEM*100,'o','MarkerSize',6,'MarkerFaceColor',...
    'g','MarkerEdgeColor','g','LineWidth',2,'Color','g','CapSize',6); hold on;

% Formatting
axis([0 n*60 0 6.2]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',18); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:1.2:6.2));
ytickformat('%.1f');
xlabel("Time (min)"); ylabel("Concentration x 10^2 (mg/mL)")
h=legend("Blood Outflow Model (SF)","Blood Outflow Model (SFG)","Plasma Outflow Model (SF)",...
    "Plasma Outflow Model (SFG)","RBC Outflow Model (SF)","RBC Outflow Model (SFG)",...
    "Hepatocyte Model (SF)","Hepatocyte Model (SFG)","Bile Model (SF)","Bile Model (SFG)",...
    "Bile Data (SF)","Location","NorthWest",'NumColumns',2);
set(h, "fontSize",14); legend("boxoff")

%% Plot model and data
figure(3); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(t,Cin*100,'-r','LineWidth',2.5); hold on;
plot(t,Cpl_In*100,'-b','LineWidth',2.5); hold on;
plot(t,Crbc_In*100,'-c','LineWidth',2.5); hold on;

plot(t,Cin_C*100,'--r','LineWidth',2.5); hold on;
plot(t,Cin_C_pl*100,'--b','LineWidth',2.5); hold on;
plot(t,Cin_C_rbc*100,'--c','LineWidth',2.5); hold on;

plot(T, C1M*100, ':r', 'LineWidth', 2.5); hold on;
plot(T, CplM*100, ':b', 'LineWidth', 2.5); hold on;
plot(T, CrbcM*100, ':c', 'LineWidth', 2.5); hold on;
% SF Outflow in Blood
plot(T, C2M*100, '-.r', 'LineWidth', 2.5); hold on;
% SF_C Outflow Blood
plot(T, CplM_C*100, '-.b', 'LineWidth', 2.5); hold on;
% SF_C Outflow Plasma
plot(T, CrbcM_C*100, '-.c', 'LineWidth', 2.5); hold on;
% SF_C Outflow RBC
plot(T, C3M*100, '-m', 'LineWidth', 2.5); hold on;
% SF in Hepatocytes
plot(T, C4M*100, '-.m', 'LineWidth', 2.5); hold on;
% SF_C in Hepatocytes
plot(T, C5M*100, '-g', 'LineWidth', 2.5); hold on;
%SF in Bile
plot(T, C6M*100, '-.g', 'LineWidth', 2.5); hold on;
%SF_C in Bile

errorbar(tdata,BileData*100,Bile_SEM*100,'o','MarkerSize',6,'MarkerFaceColor',...
    'g','MarkerEdgeColor','g','LineWidth',2,'Color','g','CapSize',6); hold on;

% Formatting
axis([0 n*60 0 6.2]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',18,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:1.2:6.2));
ytickformat('%.2f');
xlabel("Time (min)"); ylabel("Concentration x 10^2 (mg/mL)")
h=legend("Blood Inflow Model (SF)","Plasma Inflow Model (SF)","RBC Inflow Model (SF)",...
    "Blood Inflow Model (SFG)","Plasma Inflow Model (SFG)","RBC Inflow Model (SFG)",...
    "Blood Outflow Model (SF)","Plasma Outflow Model (SF)","RBC Outflow Model (SF)",...
    "Blood Outflow Model (SFG)","Plasma Outflow Model (SFG)","RBC Outflow Model (SFG)",...
    "Hepatocyte Model (SF)","Hepatocyte Model (SFG)","Bile Model (SF)","Bile Model (SFG)",...
    "Bile Data (SF)","Location","NorthWest",'NumColumns',2);
set(h, "fontSize",12,'FontName','Times New Roman'); legend("boxoff")

%% Manuscript Figure - 5C
figure(4); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(T, C1M*100, '-r', 'LineWidth', 2.5); hold on;
% SF Outflow in Blood
plot(T, C2M*100, ':r', 'LineWidth', 2.5); hold on;
% SF_C Outflow Blood
plot(T, CplM*100, '-b', 'LineWidth', 2.5); hold on;
% SF Outflow in Plasma
plot(T, CplM_C*100, ':b', 'LineWidth', 2.5); hold on;
% SF_C Outflow Plasma
plot(T, CplM_C*100, ':w', 'LineWidth', 0.1); hold on;
% SF_C Outflow Plasma (Fake)

plot(T, C3M*100, '-m', 'LineWidth', 2.5); hold on;
% SF in Hepatocytes
plot(T, C4M*100, ':m', 'LineWidth', 2.5); hold on;
% SF_C in Hepatocytes
plot(T, C5M*100, '-g', 'LineWidth', 2.5); hold on;
%SF in Bile
plot(T, C6M*100, ':g', 'LineWidth', 2.5); hold on;
%SF_C in Bile

errorbar(tdata,BileData*100,Bile_SEM*100,'o','MarkerSize',6,'MarkerFaceColor',...
    'g','MarkerEdgeColor','g','LineWidth',2,'Color','g','CapSize',6); hold on;

% Formatting
axis([0 n*60 0 6.2]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:1.2:6.2));
ytickformat('%.1f');
xlabel("Time (min)"); ylabel("Concentration x 10^2 (mg/mL)")
% h=legend("Blood outflow model (SF)","Blood outflow model (SFG)","Plasma outflow model (SF)",...
%     "Plasma outflow model (SFG)"," ","Hepatocyte model (SF)","Hepatocyte model (SFG)",...
%     "Bile model (SF)","Bile model (SFG)","Bile data (SF)","Location","NorthWest",'NumColumns',2);
% set(h, "fontSize",15,'FontName','Times New Roman'); legend("boxoff")

%% Manuscript Figure - 5A
figure(5); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(T, C1M*100, '-r', 'LineWidth', 2.5); hold on;
% SF Outflow in Blood
plot(T, C2M*100, ':r', 'LineWidth', 2.5); hold on;
% SF_C Outflow Blood
plot(T, CplM*100, '-b', 'LineWidth', 2.5); hold on;
% SF Outflow in Plasma
plot(T, CplM_C*100, ':b', 'LineWidth', 2.5); hold on;
% SF_C Outflow Plasma

% Formatting
axis([0 60 0 2.45]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:10:60),'YTick',(0:0.6:2.4));
ytickformat('%.1f');
xlabel("Time (min)"); ylabel("Concentration x 10^2 (mg/mL)");

h=legend("Blood outflow model (SF)","Blood outflow model (SFG)","Plasma outflow model (SF)",...
    "Plasma outflow model (SFG)",'NumColumns',1);
set(h, "fontSize",15,'FontName','Times New Roman','Location','SouthEast'); 
legend("boxoff")

%% Manuscript Figure - 5B
figure(6); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(T, C3M*100, '-m', 'LineWidth', 2.5); hold on;
% SF in Hepatocytes
plot(T, C4M*100, ':m', 'LineWidth', 2.5); hold on;
% SF_C in Hepatocytes
plot(NaN,NaN,':w','LineWidth',0.5); hold on;
%Fake plot
plot(T, C5M*100, '-g', 'LineWidth', 2.5); hold on;
%SF in Bile
plot(T, C6M*100, ':g', 'LineWidth', 2.5); hold on;
%SF_C in Bile

errorbar(tdata,BileData*100,Bile_SEM*100,'o','MarkerSize',6,'MarkerFaceColor',...
    'g','MarkerEdgeColor','g','LineWidth',2,'Color','g','CapSize',6); hold on;

% Formatting
axis([0 60 0 6.2]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:10:60),'YTick',(0:1.2:6.2));
ytickformat('%.1f');
xlabel("Time (min)"); ylabel("Concentration x 10^2 (mg/mL)");

h=legend("Hepatocyte model (SF)","Hepatocyte model (SFG)","",...
    "Bile model (SF)","Bile model (SFG)","Bile data (SF)","Location","NorthWest",'NumColumns',1);
set(h, "fontSize",15,'FontName','Times New Roman'); legend("boxoff")