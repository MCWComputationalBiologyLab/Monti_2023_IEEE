%%========================================================================
% Fitting input function to whole blood data involving variations on a 
% triple exponential (TE) model. 
% 
% Rationale: Experimentally, blood is sampled from the carotid artery and 
% hence Cin, in Model 1, is being sampled rather than C1 (venous blood
% leaving the liver). Thus, we should fit Cin to whole blood data rather
% than C1. Co is the actual drug input into the jugular vein, which then,
% in Model 1, passes directly to the liver. 
%
% Christopher Monti 9/22/22
%%========================================================================

clear; clc; close all; clc;
format short e
rng(1);

%%========================================================================
%% Load and parse out data
C1 = load('Control_Pub_Blood_C_Data_NLR_Avg_SD_SEM.txt');
C1_SEM = C1(:,4);
C1_SD = C1(:,3); %C1_SEM = C1_SD;
C1data = C1(:,2);
tdata = C1(:,1);

%%========================================================================
%% Fitting model to data
options = optimset('fmincon');
%options = optimset(options,'algorithm','sqp');
options = optimset(options,'TolFun',1e-6,'TolX',1e-6);
%options = optimset(options,'Display','iter');
%options = optimset(options,'Maxiter',100,'MaxFunEvals',1000);
%options = optimset(options,'UseParallel','always');
%options = []; % We want to use only default optimizer options

numIter = 10;
mpars = []; fvals = [];

randLb = 0.5; randUb = 1.5;

p0 = [2.5E-2,1E0,6E-3,1E-1,3E-3,4E-3];

for i = 1:numIter
    temp = p0;
    randNum = (randLb + (randUb-randLb)*rand(1,6));
    %Generates random number b/w 0.5 and 1.5 effectively allowing p0 to
    %vary between 0.5p0 and 1.5p0 (a maximal increase or decrease of 50%) 

    p0 = p0.*randNum;
    lb = p0./100; ub = p0.*100;
    p0(1) = 2.5E-2; lb(1) = p0(1); ub(1) = p0(1); % fix p0
    [mpar,fval] = fmincon(@Error,p0,[],[],[],[],lb,ub,@mycon,...
        options,tdata,C1data);
    
    mpars(i,:) = mpar;
    fvals(i) = fval;
    p0 = temp;
end

[minFval,minIdx] = min(fvals);
TE_mpar = mpars(minIdx,:);
minFval
TE_mpar

save('Input_TE_mpar_Ctrl','TE_mpar');

N = 8; P = 5; K = P+1; E = minFval;
AIC = N*log(E/N) + 2*K + (2*K*(K+1))/(N-K-1)

%%========================================================================
%% Simulate the model with estimated parameters
load('Input_TE_mpar_Ctrl.mat');
%TE_mpar = [2.5E-2,1.36E0,6.26E-3,1.23E-1,2.75E-3,3.77E-3];

n = 20;
dt = 0.01; tend = n*60;
t = (0:dt:tend)';
Cin = myfunc(TE_mpar,t);
AUC = dt*trapz(Cin)

figure(1); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');
plot(t,Cin*100,'-r','LineWidth',2.5); hold on;
errorbar(tdata,C1data*100,C1_SEM*100,'s','MarkerSize',6,'MarkerFaceColor','r',...
    'MarkerEdgeColor','r','LineWidth',2,'Color','r','CapSize',6); hold on;
axis([0 n*60 0 3]); box on; grid off
%set(gca,'LineWidth',1.5,'FontSize',18,'FontWeight','bold'); box off;
set(gca,'LineWidth',1.5,'FontSize',18); box off;
set(gca,'XTick',n*(0:10:60),'YTick',(0:0.6:3));
xlabel("Time (min)"); ylabel("100*SF Concentration (mg/mL)")
title('Triple Exponential Decay','FontWeight','normal')

%%========================================================================
%% Error function for optimization
function Err = Error(mpar,tdata,C1data)

C1model = myfunc(mpar,tdata);
Err = sum(((C1data - C1model)/max(C1data)).^2);

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
%% Nonlinear constraint relating parameters for opimization
function [NL_ineq,NL_eq] = mycon(x,~,~,~,~)
A1 = x(1);
k1 = x(2);
A2 = x(3);
k2 = x(4);
A3 = x(5);
k3 = x(6);

p = Parameters();

t_peak = p.TBV/p.CO; 
C_peak = p.Dose*p.Weight/p.TBV;

Cin_peak = A1*exp(-k1*t_peak)+A2*exp(-k2*t_peak)+A3*exp(-k3*t_peak);

NL_ineq = 0;
NL_eq = [C_peak-Cin_peak];

end

%%========================================================================
%% Function defining parameters for liver input function
function p = Parameters()

p.CO = 85.05;     % Cardiac output (mL/min)
p.TBV = 20.70;    % Total blood volume (mL)
p.Weight = 0.277; % Body weight (kg)
p.Dose = 2;       % Amount of dye added (mg/kg)

end
