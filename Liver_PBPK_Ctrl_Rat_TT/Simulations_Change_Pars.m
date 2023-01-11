%%========================================================================
close all; clear all; clc;

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
%% Simulate model for decreasing Vmax1
% Add in legend
load('Ctrl_mpar.mat');

% Physiological parameters
p = Parameters(mpar);
Vrbc = p.Vrbc;
Vpl = p.Vpl;
lam_max = p.lamda_max;
lam_Km = p.lamda_Km;
lam_nH = p.lamda_nH;

% Simulation step size (dt) and maximum simulation time (tmax):
n = 1;
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

figure(1); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');
plot(T, C1M*100, '-r', 'LineWidth', 2.5); hold on;
% SF Outflow in Blood
plot(T, C2M*100, ':r', 'LineWidth', 2.5); hold on;
% SF_C Outflow Blood
plot(T, CplM*100, '-b', 'LineWidth', 2.5); hold on;
% SF Outflow in Plasma
plot(T, CplM_C*100, ':b', 'LineWidth', 2.5); hold on;
% SF_C Outflow Plasma
plot(T, CplM_C*100, ':w', 'LineWidth', 2.5); hold on;
% Formatting 
plot(T, C3M*100, '-m', 'LineWidth', 2.5); hold on;
% SF in Hepatocytes
plot(T, C4M*100, ':m', 'LineWidth', 2.5); hold on;
% SF_C in Hepatocytes
a = plot(T, C5M*100, '-g', 'LineWidth', 2.5); hold on;
%SF in Bile
b = plot(T, C6M*100, ':g', 'LineWidth', 2.5); hold on;
%SF_C in Bile

errorbar(tdata,BileData*100,Bile_SEM*100,'o','MarkerSize',6,'MarkerFaceColor',...
    'g','MarkerEdgeColor','k','LineWidth',2,'Color','k','CapSize',6); hold on;

% Formatting
axis([0 n*60 0 6.2]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:1.2:6.0));
ytickformat('%.1f');
xlabel('Time (min)'); ylabel('Concentration x 10^2 (mg/mL)');

for j = 1:21
load('Ctrl_mpar.mat');
mpar(1) = mpar(1)*(31-j)/20; %change Vmax2

% Physiological parameters
p = Parameters(mpar);
Vrbc = p.Vrbc;
Vpl = p.Vpl;
lam_max = p.lamda_max;
lam_Km = p.lamda_Km;
lam_nH = p.lamda_nH;

% Simulation step size (dt) and maximum simulation time (tmax):
n = 1;
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
%% Manuscript Figure
figure(1); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(T, C1M*100, '-r', 'LineWidth', 2.5); hold on;
% SF Outflow in Blood
plot(T, C2M*100, ':r', 'LineWidth', 2.5); hold on;
% SF_C Outflow Blood
plot(T, CplM*100, '-b', 'LineWidth', 2.5); hold on;
% SF Outflow in Plasma
plot(T, CplM_C*100, ':b', 'LineWidth', 2.5); hold on;
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
    'g','MarkerEdgeColor','k','LineWidth',2,'Color','k','CapSize',6); hold on;

% Formatting
axis([0 n*60 0 6.2]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:1.2:6.0));
ytickformat('%.1f');
xlabel('Time (min)'); ylabel('Concentration x 10^2 (mg/mL)');
end

% Make sure black is on top
load('Ctrl_mpar.mat');

% Physiological parameters
p = Parameters(mpar);
Vrbc = p.Vrbc;
Vpl = p.Vpl;
lam_max = p.lamda_max;
lam_Km = p.lamda_Km;
lam_nH = p.lamda_nH;

% Simulation step size (dt) and maximum simulation time (tmax):
n = 1;
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

figure(1); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');
plot(T, C5M*100, '-k', 'LineWidth', 2.5); hold on;
%SF in Bile
plot(T, C6M*100, ':k', 'LineWidth', 2.5); hold on;
%SF_C in Bile

h=legend("Blood outflow model (SF)","Blood outflow model (SFG)","Plasma outflow model (SF)",...
    "Plasma outflow model (SFG)"," ","Hepatocyte model (SF)","Hepatocyte model (SFG)",...
    "Bile model (SF)","Bile model (SFG)","Bile data (SF)","Location","North",'NumColumns',2);
set(h, "fontSize",15,'FontName','Times New Roman'); legend("boxoff")

%%========================================================================

%%========================================================================
%% Simulate model for decreasing Vmax2
for j = 1:21
load('Ctrl_mpar.mat');
mpar(2) = mpar(2)*(31-j)/20; %change Vmax2

% Physiological parameters
p = Parameters(mpar);
Vrbc = p.Vrbc;
Vpl = p.Vpl;
lam_max = p.lamda_max;
lam_Km = p.lamda_Km;
lam_nH = p.lamda_nH;

% Simulation step size (dt) and maximum simulation time (tmax):
n = 1;
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
%% Manuscript Figure
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
plot(T, C3M*100, '-m', 'LineWidth', 2.5); hold on;
% SF in Hepatocytes
plot(T, C4M*100, ':m', 'LineWidth', 2.5); hold on;
% SF_C in Hepatocytes
plot(T, C5M*100, '-g', 'LineWidth', 2.5); hold on;
%SF in Bile
plot(T, C6M*100, ':g', 'LineWidth', 2.5); hold on;
%SF_C in Bile

if j == 11
    a = plot(T, C5M*100, '-k', 'LineWidth', 2.5); hold on;
    b = plot(T, C6M*100, ':k', 'LineWidth', 2.5); hold on;
end

errorbar(tdata,BileData*100,Bile_SEM*100,'o','MarkerSize',6,'MarkerFaceColor',...
    'g','MarkerEdgeColor','k','LineWidth',2,'Color','k','CapSize',6); hold on;

% Formatting
axis([0 n*60 0 6.2]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:1.2:6.0));
ytickformat('%.1f');
xlabel('Time (min)'); ylabel('Concentration x 10^2 (mg/mL)');
end

%Make sure black on top
uistack(a,'top');
uistack(b,'top');
%%========================================================================

%%========================================================================
%% Simulate model for decreasing Vmax3
for j = 1:21
load('Ctrl_mpar.mat');
mpar(3) = mpar(3)*(31-j)/20; %change Vmax3

% Physiological parameters
p = Parameters(mpar);
Vrbc = p.Vrbc;
Vpl = p.Vpl;
lam_max = p.lamda_max;
lam_Km = p.lamda_Km;
lam_nH = p.lamda_nH;

% Simulation step size (dt) and maximum simulation time (tmax):
n = 1;
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
%% Manuscript Figure
figure(3); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(T, C1M*100, '-r', 'LineWidth', 2.5); hold on;
% SF Outflow in Blood
plot(T, C2M*100, ':r', 'LineWidth', 2.5); hold on;
% SF_C Outflow Blood
plot(T, CplM*100, '-b', 'LineWidth', 2.5); hold on;
% SF Outflow in Plasma
plot(T, CplM_C*100, ':b', 'LineWidth', 2.5); hold on;
% SF_C Outflow Plasma
plot(T, C3M*100, '-m', 'LineWidth', 2.5); hold on;
% SF in Hepatocytes
plot(T, C4M*100, ':m', 'LineWidth', 2.5); hold on;
% SF_C in Hepatocytes
plot(T, C5M*100, '-g', 'LineWidth', 2.5); hold on;
%SF in Bile
plot(T, C6M*100, ':g', 'LineWidth', 2.5); hold on;
%SF_C in Bile

if j == 11
    plot(T, C5M*100, '-k', 'LineWidth', 2.5); hold on;
    plot(T, C6M*100, ':k', 'LineWidth', 2.5); hold on;
end

errorbar(tdata,BileData*100,Bile_SEM*100,'o','MarkerSize',6,'MarkerFaceColor',...
    'g','MarkerEdgeColor','k','LineWidth',2,'Color','k','CapSize',6); hold on;

% Formatting
axis([0 n*60 0 9.3]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:2.25:9));
ytickformat('%.1f');
xlabel('Time (min)'); ylabel('Concentration x 10^2 (mg/mL)');
end

%Make sure black on top
uistack(a,'top');
uistack(b,'top');
%%========================================================================