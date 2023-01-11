%%========================================================================
close all; clear all; clc;
rng(1); % For reproducibility

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
%% fmincon optimization
options = optimset('fmincon');
%options = optimset(options,'algorithm','sqp');
options = optimset(options,'TolFun',1e-8,'TolX',1e-8);
options = optimset(options,'Display','iter');
%options = optimset(options,'Maxiter',100,'MaxFunEvals',1000);
%options = optimset(options,'UseParallel','always');
%options = []; % We want to use only default optimizer options

numIters = 1000;
mpars = []; fvals = [];

randLb = 0.5; randUb = 1.5;

p0 = [2.4E-2,4.3E-1,2.7E-3];
numPars = length(p0);

parfor i = 1:numIters
    p = p0;
    randNum = (randLb + (randUb-randLb)*rand(1,numPars));
    %Generates random number b/w 0.5 and 1.5 effectively allowing p0 to
    %vary between 0.5p0 and 1.5p0 (a maximal increase or decrease of 50%) 

    p = p.*randNum;
    lbp = p/100; ubp = p*100;
%     p0(1) = 2.4E-2; lbp(1) = p0(1)/1; ubp(1) = p0(1)*1;
%     p0(2) = 4.3E-1; lbp(2) = p0(2)/1; ubp(2) = p0(2)*1;
%     p0(3) = 2.7E-3; lbp(3) = p0(3)/1; ubp(3) = p0(3)*1;
    [mpar,fval] = fmincon(@Error,p,[],[],[],[],lbp,ubp,[],...
        options,tdata,BileData);
    
    mpars(i,:) = mpar;
    fvals(i) = fval;
end

[minFval,minIdx] = min(fvals);
mpar = mpars(minIdx,:);
minFval
mpar

save('Ctrl_mpar.mat','mpar');
save('Ctrl_resid.mat','minFval');
%% =======================================================================
figure(5); set(gcf,'Units','inches','Position',[0.5 0.5 7 5]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 7 5],'color','white');
for i = 1:numPars
    subplot(1,3,i);
    histogram(mpars(:,i)); hold on
end

%%========================================================================
%% Simulate model and plot results
Simulations();