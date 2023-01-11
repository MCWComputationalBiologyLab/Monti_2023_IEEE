%%========================================================================
% Using fitted Cin for FREE SF (SF_F) to determine Cin for SF CONJUGATED to
% glcuronic acid (SF_C)
% 
% Rationale: SF_F is measured in blood (Cin) via fluorescence. At the 
% wavelength used, SF_C has very little fluorescence signal. Therefore,
% given the ratios of SF_F to SF_C (see Dr. Kim email), we can deduce the
% amount of SF_C at any time given SF_F from Cin.
%
% Christopher Monti 12/2/22
%%========================================================================

clear; clc; close all; clc;
format short e
rng(1);

%%========================================================================
%% Load and parse out data for Cin (for comparison purposes)
C1 = load('Control_Pub_Blood_C_Data_NLR_Avg_SD_SEM.txt');
C1_SEM = C1(:,4);
C1_SD = C1(:,3); %C1_SEM = C1_SD;
C1data = C1(:,2);
tdata = C1(:,1);

C1_pl = load('Control_Pub_Plasma_C_Data_NLR_Avg_SD_SEM.txt');
C1_pl_SEM = C1_pl(:,4);
C1_pl_SD = C1_pl(:,3);
C1_pl_data = C1_pl(:,2);
C1_pl_tdata = C1_pl(:,1);

%% Set data for ratio between SF_F and SF_C
SF_C_Ratio_Data = [0,15,30,60,120,180,240,300,360,420,480;0,0.13,0.33,...
    0.68,0.80,0.87,0.92,0.93,0.94,0.94,0.94]; %from Dr. Kim email
SF_C_Ratio_Data(1,:) = SF_C_Ratio_Data(1,:)/4; %Account for transit time

%%========================================================================
%% Fitting model to data
options = optimset('fmincon');
%options = optimset(options,'algorithm','sqp');
options = optimset(options,'TolFun',1e-6,'TolX',1e-6);
%options = optimset(options,'Display','iter');
%options = optimset(options,'Maxiter',100,'MaxFunEvals',1000);
%options = optimset(options,'UseParallel','always');
%options = []; % We want to use only default optimizer options

numIter = 1;
mpars = []; fvals = [];

randLb = 0.5; randUb = 1.5;

p0 = [1,20,2];

for i = 1:numIter
    temp = p0;
    randNum = (randLb + (randUb-randLb)*rand(1,3));
    %Generates random number b/w 0.5 and 1.5 effectively allowing p0 to
    %vary between 0.5p0 and 1.5p0 (a maximal increase or decrease of 50%) 

    p0 = p0.*randNum;
    lb = p0./100; ub = p0.*100;
    [mpar,fval] = fmincon(@Error,p0,[],[],[],[],lb,ub,[],...
        options,SF_C_Ratio_Data(1,:),SF_C_Ratio_Data(2,:));
    
    mpars(i,:) = mpar;
    fvals(i) = fval;
    p0 = temp;
end

[minFval,minIdx] = min(fvals);
SF_C_mpar = mpars(minIdx,:);
minFval
SF_C_mpar

save('Input_SF_C_mpar_Ctrl','SF_C_mpar');

%% Visualize SF_C_Ratio fit to data
n = 2;
tspan = (0:0.5:n*60);
model = SF_C_Ratio(SF_C_mpar,tspan);

figure(1)
set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');
plot(SF_C_Ratio_Data(1,:),SF_C_Ratio_Data(2,:),'o','MarkerSize',9,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','Visible','on'); hold on;
plot(SF_C_Ratio_Data(1,:),1-SF_C_Ratio_Data(2,:),'s','MarkerSize',9,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','Visible','on'); hold on;
plot(tspan,model,'-.k','LineWidth',2.5,'Visible','on'); hold on;
plot(tspan,1-model,'-k','LineWidth',2.5,'Visible','on'); hold on;
axis([0 n*60 0 1]); box on; grid off
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',(0:n*10:n*60),'YTick',(0:0.2:1));
xlabel("Time (min)"); ylabel("Fraction (Unitless)")
h = legend('Plasma SFG','Plasma SF','','','FontName','Times New Roman',...
    'FontSize',16,'NumColumns',1,'Box','off','Location','East');

%%========================================================================
%% Simulate the model with estimated parameters
load('Input_TE_mpar_Ctrl.mat');
%TE_mpar = [2.5E-2,1.36E0,6.26E-3,1.23E-1,2.75E-3,3.77E-3];

n = 4;
dt = 0.01; tend = n*60;
t = (0:dt:tend)';
Cin = myfunc(TE_mpar,t); %SF_F in blood

%Calculate plasma
Vrbc = 8.0424e-1;
Vpl = 1.2064e0;
lam_max = 9;
lam_Km = 5e-3;
lam_nH = 3;
lamda = 0.5 + lam_max*Cin.^lam_nH./(lam_Km^lam_nH+Cin.^lam_nH); % partition coefficient
Crbc_In = ((Vrbc+Vpl)./(Vrbc+lamda.*Vpl)).*Cin; %concentration in rbc based on concentration in blood
Cpl_In = ((Vrbc+Vpl)./(Vrbc./lamda+Vpl)).*Cin; %concentration in plasma based on concentration in blood

%Ratio relates plasma SF to SFG
One_miunus_alpha = SF_C_Ratio(SF_C_mpar,t); % ratio of SF_F to SF_C
denom = 1./(One_miunus_alpha)-1; %for reading clarity
SF_C_pl = Cpl_In./denom;
SF_C_bl = (Vrbc./lamda+Vpl).*SF_C_pl/(Vrbc+Vpl);

%Convert data for plasma from SF to SFG
C1_pl_SFG = [];
C1_bl_SFG = [];

for i = 1:length(C1_pl_data)
    One_miunus_alpha = SF_C_Ratio(SF_C_mpar,tdata(i)); % ratio of SF_F to SF_C
    denom = 1./(One_miunus_alpha)-1; %for reading clarity
    plTemp = C1_pl_data(i)./denom;
    lamda = 0.5 + lam_max*myfunc(TE_mpar,tdata(i))^lam_nH/(lam_Km^lam_nH+myfunc(TE_mpar,tdata(i))^lam_nH); % partition coefficient
    blTemp = (Vrbc/lamda+Vpl)*plTemp/(Vrbc+Vpl);
    C1_pl_SFG = [C1_pl_SFG,plTemp];
    C1_bl_SFG = [C1_bl_SFG,blTemp];
end

%Plot
figure(2); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');
plot(t,Cin*100,'-r','LineWidth',2.5); hold on;
errorbar(tdata,C1data*100,C1_SEM*100,'o','MarkerSize',6,'MarkerFaceColor','r',...
    'MarkerEdgeColor','r','LineWidth',2,'Color','r','CapSize',6); hold on;
% plot(tdata,C1_bl_SFG*100,'s','MarkerSize',6,'MarkerFaceColor','r',...
%     'MarkerEdgeColor','r','LineWidth',2,'Color','r'); hold on;
plot(t,SF_C_bl*100,':r','LineWidth',3.0);

plot(t,Cpl_In*100,'-b','LineWidth',2.5); hold on;
errorbar(C1_pl_tdata,C1_pl_data*100,C1_pl_SEM*100,'o','MarkerSize',6,'MarkerFaceColor','b',...
    'MarkerEdgeColor','b','LineWidth',2,'Color','b','CapSize',6); hold on;
% plot(C1_pl_tdata,C1_pl_SFG*100,'s','MarkerSize',6,'MarkerFaceColor','b',...
%     'MarkerEdgeColor','b','LineWidth',2,'Color','b'); hold on;
plot(t,SF_C_pl*100,':b','LineWidth',3.0);
axis([0 n*60 0 2.8]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',20,'FontName','Times New Roman'); box off;
set(gca,'XTick',n*(0:10:60),'YTick',(0:0.7:2.8));
ytickformat('%.1f');
xlabel("Time (min)"); ylabel("Concentration x 10^2 (mg/mL)")
% h = legend('Blood model (SF)','Blood data (SF)','Blood data (SFG)',...
%     'Blood model (SFG)','Plasma model (SF)','Plasma data (SF)',...
%     'Plasma data (SFG)','Plasma model (SFG)','FontName','Arial',...
%     'FontSize',14,'NumColumns',2,'Box','off');

h = legend('Blood inflow model (SF)', 'Blood data (SF)',...
    'Blood inflow model (SFG)','Plasma inflow model (SF)','Plasma data (SF)','Plasma inflow model (SFG)','FontName','Times New Roman',...
    'FontSize',16,'NumColumns',1,'Box','off','Location','North');

% Check my work
% figure(3); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
% set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');
% total = SF_C_pl+Cpl_In;
% ratio = SF_C_pl./total;
% plot(t,ratio,'-k','LineWidth',2.5); hold on;
% plot(SF_C_Ratio_Data(1,:),SF_C_Ratio_Data(2,:),'o','MarkerSize',6,...
%     'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
% title("Check");

%%========================================================================
%% Error function for optimization
function Err = Error(mpar,tdata,data)

model = SF_C_Ratio(mpar,tdata);
Err = sum(((data - model)/max(data)).^2);

end

%%========================================================================
%% Double exponential with ramp function up to six parameters
function Cin = myfunc(TE_mpar,tdata) 
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
%% Function defining parameters for liver input function
function p = Parameters()

p.CO = 85.05;     % Cardiac output (mL/min)
p.TBV = 20.70;    % Total blood volume (mL)
p.Weight = 0.277; % Body weight (kg)
p.Dose = 2;       % Amount of dye added (mg/kg)

end

%%========================================================================
%% Function for SF_C Ratio
function SF_C_R = SF_C_Ratio(mpar,tdata) 
Vmax = mpar(1);
Km = mpar(2);
n = mpar(3);

SF_C_R = Vmax*tdata.^n./(Km^n+tdata.^n);

end