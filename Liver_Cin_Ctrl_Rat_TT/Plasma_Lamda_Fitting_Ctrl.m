%%========================================================================
% Fitting lamda function to parse out plasma dye concentration
%
% Christopher Monti 10/18/22
%%========================================================================

clear; clc; close all; clc;
format short e
rng(1);

%%========================================================================
%% Load and parse out data
C1 = load('Control_Pub_Blood_C_Data_NLR_Avg_SD_SEM.txt');
C1_SEM = C1(:,4);
C1_SD = C1(:,3);  %C1_SEM = C1_SD;
C1data = C1(:,2);
tdata = C1(:,1);

Cpl = load('Control_Pub_Plasma_C_Data_NLR_Avg_SD_SEM.txt');
C2_SEM = Cpl(:,4);
C2_SD = Cpl(:,3);  %C2_SEM = C2_SD;
C2data = Cpl(:,2);

%%========================================================================
%% Fitting model to data
options = optimset('fmincon');
%options = optimset(options,'algorithm','sqp');
options = optimset(options,'TolFun',1e-6,'TolX',1e-6);
options = optimset(options,'Display','iter');
%options = optimset(options,'Maxiter',100,'MaxFunEvals',1000);
%options = optimset(options,'UseParallel','always');
%options = []; % We want to use only default optimizer options

numIter = 1;
mpars = []; fvals = [];

randLb = 0.5; randUb = 1.5;

p0 = [9,5e-3];

for i = 1:numIter
    temp = p0;
    randNum = (randLb + (randUb-randLb)*rand(1,2));
    %Generates random number b/w 0.5 and 1.5 effectively allowing p0 to
    %vary between 0.5p0 and 1.5p0 (a maximal increase or decrease of 50%) 

    p0 = p0.*randNum;
    lb = p0./10; ub = p0.*10;
    p0(1) = 9; lb(1) = p0(1); ub(1) = p0(1); % fix p0
    [mpar,fval] = fmincon(@Error_Lamda,p0,[],[],[],[],lb,ub,[],...
        options,tdata,C2data);
    
    mpars(i,:) = mpar;
    fvals(i) = fval;
    p0 = temp;
end

[minFval,minIdx] = min(fvals);
lamda_mpar = mpars(minIdx,:);
minFval
lamda_mpar

save('Plasma_Lamda_mpar_Ctrl.mat','lamda_mpar');

%%========================================================================
%% Simulate the model with estimated parameters and plot results
% Calculate C in blood (input), C in RBCs, and C in plasma
load('Plasma_Lamda_mpar_Ctrl.mat')

% Simulation step size (dt) and maximum simulation time (tmax):
n = 1;
dt = 0.01; 
tmax = n*60;
t = (0:dt:tmax)';

% Model simulation
Cin = Blood_Input(t);
[Cpl,Crbc,lamda] = myfunc(lamda_mpar,t);

% Plot model simulations
figure(1); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');

plot(t,Cin*100,'-r','LineWidth',2.5); hold on;
plot(t,Cpl*100,'-b','LineWidth',2.5); hold on;
plot(t,Crbc*100,'-c','LineWidth',2.5); hold on;
errorbar(tdata,C1data*100,C1_SEM*100,'o','MarkerSize',6,'MarkerFaceColor','r',...
    'MarkerEdgeColor','r','LineWidth',2,'Color','r','CapSize',6); hold on;
errorbar(tdata,C2data*100,C2_SEM*100,'o','MarkerSize',6,'MarkerFaceColor','b',...
    'MarkerEdgeColor','b','LineWidth',2,'Color','b','CapSize',6); hold on;

axis([0 n*60 0 4.5]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',n*(0:10:60),'YTick',(0:0.9:4.5));
xlabel("Time (min)"); ylabel(sprintf("Concentration x 10^2 (mg/mL)"))
% h=legend("Arterial Blood Model","Arterial Plasma Model","Arterial RBC Model",...
%     "Arterial Blood Data","Arterial Plasma Data","Location","NorthEast");
h=legend("Blood inflow model","Plasma inflow model","RBC inflow model",...
    "Blood inflow data","Plasma inflow data","Location","North");
set(h, "fontSize",16,'FontName','Times New Roman'); legend("boxoff")

figure(2); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');
plot(t,lamda,'-k','linewidth',2.5); hold on
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
xlabel('Time (min)'); ylabel('Partition Coefficient (Unitless)');
xlim(n*[0,60]); ylim([0,10])

%%========================================================================
%% Error function for optimization
function Err = Error_Lamda(mpar,tdata,C2data)

C2model = myfunc(mpar,tdata);
Err = sum(((C2data - C2model)/max(C2data)).^2);

end

%%========================================================================
%% Compute plasma and RBC concentrations from blood concentration (Cin)
function varargout = myfunc(lamda_mpar,tdata)

% Necessary data
Cin = Blood_Input(tdata);
p = Parameters();
Vrbc = p.Vrbc;
Vpl = p.Vpl;

% Lamda parameters
lam_max = lamda_mpar(1);
lam_Km = lamda_mpar(2);
lam_nH = 3; % Hill coefficient
lamda = 0.5+lam_max*Cin.^lam_nH./(lam_Km^lam_nH+Cin.^lam_nH);

% Calculating plasma and RBC concentrations of the dye
Cpl = ((Vrbc+Vpl)./(Vrbc./lamda+Vpl)).*Cin; %concentration in plasma based on concentration in blood
Crbc = ((Vrbc+Vpl)./(Vrbc+lamda.*Vpl)).*Cin; %concentration in rbc based on concentration in blood

if (nargout == 1)
    varargout = {Cpl};
else
    varargout = {Cpl,Crbc,lamda};
end

end

%%========================================================================
%% Double exponential with ramp function up to six parameters
function Cin = Blood_Input(tdata)

load('Input_TE_mpar_Ctrl.mat');
%TE_mpar = [2.5E-2,1.36E0,6.26E-3,1.23E-1,2.75E-3,3.77E-3];

A1 = TE_mpar(1);
k1 = TE_mpar(2);
A2 = TE_mpar(3);
k2 = TE_mpar(4);
A3 = TE_mpar(5);
k3 = TE_mpar(6);

p = Parameters();

t_peak = p.TBV/p.CO; % Transit time through circulation
C_peak = p.Dose*p.Weight/p.TBV; % Peak dye concentration at transit time

for i = 1:length(tdata)
    if tdata(i) <= 0
        Cin(i) = 0;
    elseif tdata(i) > 0 && tdata(i) <= t_peak
        Slope = C_peak/t_peak; % Slope of linear ramp
        Cin(i) = tdata(i)*Slope; % Linear ramp to C_peak
    else
        Cin(i) = A1*exp(-k1*tdata(i))+A2*exp(-k2*tdata(i))+...
            A3*exp(-k3*tdata(i));
    end
end
Cin = Cin';

end

%%========================================================================
%% Function defining parameters for liver input function (blood and plasma)
function p = Parameters()

p.CO = 85.05;     % Cardiac output (mL/min)
p.TBV = 20.70;    % Total blood volume (mL)
p.Weight = 0.277; % Body weight (kg)
p.Dose = 2;       % Amount of dye added (mg/kg)

p.Vtot = 11.17;         %volume of total liver of rat, units: ml
p.fV1 = 0.18;           %fractional volume of liver that is capillary
p.fVrbc = 0.40;         %fractional volume of rbc in capillary
p.fVpl = 0.60;          %fractional volume of plasma in capillary
p.V1 = p.Vtot*p.fV1;    %volume of capilary; units: ml
p.Vrbc = p.fVrbc*p.V1;  %volume of red blood cells; units: ml
p.Vpl = p.fVpl*p.V1;    %volume of plasma compartment; units: ml

end
