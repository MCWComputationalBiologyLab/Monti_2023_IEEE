%%=========================================================================
clear; clc; close all; 
format compact;

%% Obtaining Fluorescence Data in Concentration Space
cd 'Fluorescence_Data_Indiv'

list = dir;

fileNames = [];
for i = 3:length(list) %ignores . and ..
    fileNames = [fileNames,convertCharsToStrings(list(i).name)];
end

sort(fileNames);

for i = 1:(length(fileNames)/3) 
    BileData = load(fileNames(3*i-2));
    BileData = BileData(:,2:end);
    
    BloodData = load(fileNames(3*i-1));
    tData = BloodData(:,1);
    BloodData = BloodData(:,2:end);
    cond = split(fileNames(3*i-1),"_Blood");
    cond = cond(1);
    
    PlasmaData = load(fileNames(3*i));
    PlasmaData = PlasmaData(:,2:end);
    
    % Best fit models
    pars_Bl = [-145.8,7.769E5];
    pars_Pl = [3.992E4,0.02224];
    pars_Bi = [2.636E4,0.04628];
    
    Bl_C_BF = undoLinear(BloodData,pars_Bl);
    Pl_C_BF = undoNonLinear(PlasmaData,pars_Pl);
    Bi_C_BF = undoNonLinear(BileData,pars_Bi);

    % Obtain average and SD, then arrange data
    Bl_C_BF_Avg = mean(Bl_C_BF,2,'omitnan');
    Pl_C_BF_Avg = mean(Pl_C_BF,2,'omitnan');
    Bi_C_BF_Avg = mean(Bi_C_BF,2,'omitnan');
    
    Bl_C_BF_SD = std(Bl_C_BF,0,2,'omitnan');
    Pl_C_BF_SD = std(Pl_C_BF,0,2,'omitnan');
    Bi_C_BF_SD = std(Bi_C_BF,0,2,'omitnan');

    Bl_C_BF_SEM = [];
    for j = 1:length(BloodData(:,1))
        Bl_C_BF_SEM(j,:) = Bl_C_BF_SD(j)/sqrt(length(BloodData(~isnan(BloodData(j,:)))));
    end

    Pl_C_BF_SEM = [];
    for j = 1:length(PlasmaData(:,1))
        Pl_C_BF_SEM(j,:) = Pl_C_BF_SD(j)/sqrt(length(PlasmaData(~isnan(PlasmaData(j,:)))));
    end

    Bi_C_BF_SEM = [];
    for j = 1:length(BileData(:,1))
        Bi_C_BF_SEM(j,:) = Bi_C_BF_SD(j)/sqrt(length(BileData(~isnan(BileData(j,:)))));
    end
    
    Blood_Cal_BF = [tData, Bl_C_BF_Avg, Bl_C_BF_SD, Bl_C_BF_SEM];
    Plasma_Cal_BF = [tData, Pl_C_BF_Avg, Pl_C_BF_SD, Pl_C_BF_SEM];
    Bile_Cal_BF = [tData, Bi_C_BF_Avg, Bi_C_BF_SD, Bi_C_BF_SEM];
    
    Bl_C_BF = [tData,Bl_C_BF];
    Pl_C_BF = [tData,Pl_C_BF];
    Bi_C_BF = [tData,Bi_C_BF];

    % Linear Model Parameters
    pars_Bl = [-106.8,7.841E5];
    pars_Pl = [4800,5.85E5];
    pars_Bi = [1778,2.623E5];
    
    Bl_C_LF = undoLinear(BloodData,pars_Bl);
    Pl_C_LF = undoLinear(PlasmaData,pars_Pl);
    Bi_C_LF = undoLinear(BileData,pars_Bi);
    
    % Remove negatives
    Bl_C_LF(Bl_C_LF < 0) = 0;
    Pl_C_LF(Pl_C_LF < 0) = 0;
    Bi_C_LF(Bi_C_LF < 0) = 0;
    
    % Obtain average and SD, then arrange data
    Bl_C_LF_Avg = mean(Bl_C_LF,2,'omitnan');
    Pl_C_LF_Avg = mean(Pl_C_LF,2,'omitnan');
    Bi_C_LF_Avg = mean(Bi_C_LF,2,'omitnan');
    
    Bl_C_LF_SD = std(Bl_C_LF,0,2,'omitnan');
    Pl_C_LF_SD = std(Pl_C_LF,0,2,'omitnan');
    Bi_C_LF_SD = std(Bi_C_LF,0,2,'omitnan');

    Bl_C_LF_SEM = [];
    for j = 1:length(BloodData(:,1))
        Bl_C_LF_SEM (j,:) = Bl_C_LF_SD(j)/sqrt(length(BloodData(~isnan(BloodData(j,:)))));
    end

    Pl_C_LF_SEM = [];
    for j = 1:length(PlasmaData(:,1))
        Pl_C_LF_SEM (j,:) = Pl_C_LF_SD(j)/sqrt(length(PlasmaData(~isnan(PlasmaData(j,:)))));
    end

    Bi_C_LF_SEM = [];
    for j = 1:length(BileData(:,1))
        Bi_C_LF_SEM (j,:) = Bi_C_LF_SD(j)/sqrt(length(BileData(~isnan(BileData(j,:)))));
    end
    
    Blood_Cal_LF = [tData, Bl_C_LF_Avg, Bl_C_LF_SD, Bl_C_LF_SEM];
    Plasma_Cal_LF = [tData, Pl_C_LF_Avg, Pl_C_LF_SD, Pl_C_LF_SEM];
    Bile_Cal_LF = [tData, Bi_C_LF_Avg, Bi_C_LF_SD, Bi_C_LF_SEM];

    Bl_C_LF = [tData,Bl_C_LF];
    Pl_C_LF = [tData,Pl_C_LF];
    Bi_C_LF = [tData,Bi_C_LF];
    
    %% Save Data
    cd ..;
    
    
    % Save mean data
    if(exist('Concentration_Data_Avg','file') == 0)
        mkdir('Concentration_Data_Avg') 
    end

    cd 'Concentration_Data_Avg';
    
    cd ..;

    % Save individual data
    if(exist('Concentration_Data_Indiv','file') == 0)
        mkdir('Concentration_Data_Indiv') 
    end

    cd 'Concentration_Data_Indiv';
    
    figure(i);
    set(gcf,'Units','inches','Position',[0.5 0.5 10 9.85]);
    set(gcf,'Units','inches','PaperPosition',[0.5 0.5 10 9.85],'color','white');
    temp = split(cond,"_");
    plotTitle = strcat(temp(1)," ",temp(2));
    s = sgtitle(plotTitle,'FontSize',18,'FontWeight','Bold');

    subplot(3,2,2)
    errorbar(tData,mean(BloodData,2,'omitnan'),...
        std(BloodData,0,2,'omitnan')/sqrt(length(BloodData(1,:))),'-or',...
        'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','r',...
        'MarkerEdgeColor','r','CapSize',4); hold on;
    errorbar(tData,mean(PlasmaData,2,'omitnan'),...
        std(PlasmaData,0,2,'omitnan')/sqrt(length(PlasmaData(1,:))),...
        '-sb','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','b',...
        'MarkerEdgeColor','b','CapSize',4); hold on;
    errorbar(tData,mean(BileData,2,'omitnan'),...
        std(BileData,0,2,'omitnan')/sqrt(length(BileData(1,:))),'-^g',...
        'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g','CapSize',4); hold on;
    axis([0 60 0 25000]);
    set(gca,'XTick',(0:10:60),'YTick',(0:5000:25000));
    set(gca,'LineWidth',1.5,'FontSize',16); box off;
    xlabel('Time (min)'); ylabel(sprintf('Fluorescence\n(A.U.)'));
    legend('Whole Blood','Plasma','Bile','Location','NorthEast',...
        'FontSize',12); 
    legend('boxoff')
    title(sprintf('Fluorescence Mean +/- SEM\n'));

    subplot(3,2,1)
    plot(tData,BloodData,'-or','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','r','MarkerEdgeColor','r'); hold on;
    plot(tData,PlasmaData,'-sb','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','b','MarkerEdgeColor','b'); hold on;
    plot(tData,BileData,'-^g','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','g','MarkerEdgeColor','g'); hold on;
    axis([0 60 0 25000]);
    set(gca,'XTick',(0:10:60),'YTick',(0:5000:25000));
    set(gca,'LineWidth',1.5,'FontSize',16); box off;
    xlabel('Time (min)'); ylabel(sprintf('Fluorescence\n(A.U.)'));
    title(sprintf('Fluorescence Individual Data\n'));
    
    subplot(3,2,6)
    errorbar(tData,Bl_C_BF_Avg,Bl_C_BF_SEM,'-or','LineWidth',2.0,...
        'MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r',...
        'Color','r','CapSize',4); hold on;
    errorbar(tData,Pl_C_BF_Avg,Pl_C_BF_SEM,'-sb','LineWidth',2.0,...
        'MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','b',...
        'Color','b','CapSize',4); hold on
    errorbar(tData,Bi_C_BF_Avg,Bi_C_BF_SEM,'-^g','LineWidth',2.0,...
        'MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','g',...
        'Color','g','CapSize',4); hold on
    axis([0 60 0 0.07]);
    set(gca,'XTick',(0:10:60),'YTick',(0:0.014:0.07));
    set(gca,'LineWidth',1.5,'FontSize',16); box off;
    xlabel('Time (min)'); ylabel(sprintf('Concentration\n(mg/mL)'));
    title(sprintf('Nonlinear Regression Mean +/- SEM\n'));

    subplot(3,2,5)
    plot(tData,Bl_C_BF(:,2:end),'-or','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r'); hold on;
    plot(tData,Pl_C_BF(:,2:end),'-sb','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b'); hold on
    plot(tData,Bi_C_BF(:,2:end),'-^g','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','g','MarkerEdgeColor','g','Color','g'); hold on
    axis([0 60 0 0.07]);
    set(gca,'XTick',(0:10:60),'YTick',(0:0.014:0.07));
    set(gca,'LineWidth',1.5,'FontSize',16); box off;
    xlabel('Time (min)'); ylabel(sprintf('Concentration\n(mg/mL)'));
    title(sprintf('Nonlinear Regression Individual Data\n'));
    
    subplot(3,2,4)
    errorbar(tData,Bl_C_LF_Avg,Bl_C_LF_SEM,'-or','LineWidth',2.0,...
        'MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r',...
        'Color','r','CapSize',4); hold on;
    errorbar(tData,Pl_C_LF_Avg,Pl_C_LF_SEM,'-sb','LineWidth',2.0,...
        'MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','b',...
        'Color','b','CapSize',4); hold on
    errorbar(tData,Bi_C_LF_Avg,Bi_C_LF_SEM,'-^g','LineWidth',2.0,...
        'MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','g',...
        'Color','g','CapSize',4); hold on
    axis([0 60 0 0.07]);
    set(gca,'XTick',(0:10:60),'YTick',(0:0.014:0.07));
    set(gca,'LineWidth',1.5,'FontSize',16); box off;
    xlabel('Time (min)'); ylabel(sprintf('Concentration\n(mg/mL)'));
    title(sprintf('Linear Regression Mean +/- SEM\n'));
  

    subplot(3,2,3)
    plot(tData,Bl_C_LF(:,2:end),'-or','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r'); hold on;
    plot(tData,Pl_C_LF(:,2:end),'-sb','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','b','MarkerEdgeColor','b','Color','b'); hold on
    plot(tData,Bi_C_LF(:,2:end),'-^g','LineWidth',2.0,'MarkerSize',6,...
        'MarkerFaceColor','g','MarkerEdgeColor','g','Color','g'); hold on
    axis([0 60 0 0.07]);
    set(gca,'XTick',(0:10:60),'YTick',(0:0.014:0.07));
    set(gca,'LineWidth',1.5,'FontSize',16); box off;
    xlabel('Time (min)'); ylabel(sprintf('Concentration\n(mg/mL)'));
    title(sprintf('Linear Regression Individual Data\n'));
    
    % Setup for next loop and save mean fluorescence data
    cd ..;

    if(exist('Fluorescence_Data_Avg','file') == 0)
        mkdir('Fluorescence_Data_Avg') 
    end
    
    cd 'Fluorescence_Data_Avg';

    Bl_F = [tData,mean(BloodData,2,'omitnan'),std(BloodData,0,2,...
        'omitnan'),std(BloodData,0,2,'omitnan')/sqrt(length(...
        BloodData(1,:)))];
    Pl_F = [tData,mean(PlasmaData,2,'omitnan'),std(PlasmaData,0,2,...
        'omitnan'),std(PlasmaData,0,2,'omitnan')/sqrt(length(...
        PlasmaData(1,:)))];
    Bi_F = [tData,mean(BileData,2,'omitnan'),std(BileData,0,2,...
        'omitnan'),std(BileData,0,2,'omitnan')/sqrt(length(...
        BileData(1,:)))];

    % Ends up in Fluorescence_Data_Avg
    cd ..;
    cd Fluorescence_Data_Indiv;
end

%% ========================================================================
function C = undoLinear(F,pars)
    B = pars(1);
    m = pars(2);
    C = (F-B)/m;
end

function C =  undoNonLinear(F,pars)
    Vmax = pars(1);
    Km = pars(2);
    C = Km*F./(Vmax-F);
end