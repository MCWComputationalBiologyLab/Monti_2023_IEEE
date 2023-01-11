clear; clc; close all; 

% Load in data
cd 'Concentration_Data_Indiv';

list = dir;

fileNames = [];
for i = 3:length(list) %ignores . and .. and mean fluorescence folder
    fileNames = [fileNames,convertCharsToStrings(list(i).name)];
end

cd ..;

if(exist('Ratio_Data','file') == 0)
    mkdir('Ratio_Data') 
end

%% Obtain ratios in Figure 4
for i = 1:(length(fileNames)/3)
    % Get into folder containing Concentration_Data
    cd 'Concentration_Data_Indiv';

    % Loads data and obtains which experimental condition (e.g., control,
    % I30R4H, etc.)
    BileData = load(fileNames(3*i-2));
    BileData = BileData(:,2:end);
    
    BloodData = load(fileNames(3*i-1));
    tData = BloodData(:,1);
    BloodData = BloodData(:,2:end);
    cond = split(fileNames(3*i-1),"_Blood");
    cond = cond(1);
    
    PlasmaData = load(fileNames(3*i));
    PlasmaData = PlasmaData(:,2:end);

    % Take ratios shown in Figure 4. Results in indicated ratios for
    % individual biological replicates in an experimental group. 
    Bi_Bl = BileData./BloodData;
    Bi_Pl = BileData./PlasmaData;

    % Back out of Concentration_Data folder into super folder
    cd ..;
    % Get into folder where things should be saved
    cd 'Ratio_Data';
    
    Bi_Bl = [tData, Bi_Bl];
    Bi_Pl = [tData, Bi_Pl];

    writematrix(Bi_Bl,strcat(cond,"_Bile_to_Blood_Ratio_Indiv.txt"));
    writematrix(Bi_Pl,strcat(cond,"_Bile_to_Plasma_Ratio_Indiv.txt"));

    % Back out of Concentration_Data folder into super folder
    cd ..;
end

cd 'Ratio_Data';

list = dir;

fileNames = [];
for i = 3:length(list) %ignores . and .. and mean fluorescence folder
    fileNames = [fileNames,convertCharsToStrings(list(i).name)];
end

sort(fileNames);

Bi_to_Bl_Files_60 = [];
Bi_to_Bl_Files_4H = [];
Bi_to_Pl_Files_60 = [];
Bi_to_Pl_Files_4H = [];
for i = 1:(length(fileNames))
    if contains(fileNames(i),"Control") == 1
        if contains(fileNames(i),"Bl")
            Bi_to_Bl_Files_60 = [Bi_to_Bl_Files_60 fileNames(i)];
            Bi_to_Bl_Files_4H = [Bi_to_Bl_Files_4H fileNames(i)];
        else
            Bi_to_Pl_Files_60 = [Bi_to_Pl_Files_60 fileNames(i)];
            Bi_to_Pl_Files_4H = [Bi_to_Pl_Files_4H fileNames(i)];
        end
        continue;
    end

    if contains(fileNames(i),"Bl") == 1
        if contains(fileNames(i),"R60") == 1
            Bi_to_Bl_Files_60 = [Bi_to_Bl_Files_60 fileNames(i)];
        else
            Bi_to_Bl_Files_4H = [Bi_to_Bl_Files_4H fileNames(i)];
        end
    else
        if contains(fileNames(i),"R60") == 1
            Bi_to_Pl_Files_60 = [Bi_to_Pl_Files_60 fileNames(i)];
        else
            Bi_to_Pl_Files_4H = [Bi_to_Pl_Files_4H fileNames(i)];
        end
    end
end

figure(1);
set(gcf,'Units','inches','Position',[0.5 0.5 10 8]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 10 8],'color','white');
shapes = ["-^","-s","-o"];
colors = ["g","b","r"];

for i = 1:length(Bi_to_Bl_Files_60)
    currData = load(Bi_to_Bl_Files_60(i));
    tData = currData(:,1);
    Data = currData(:,2:end);
    
    subplot(2,2,1);
    errorbar(tData(2:end),mean(Data(2:end,:),2,'omitnan'),...
    std(Data(2:end,:),0,2,'omitnan')/sqrt(length(Data(1,:))),shapes(i),'Color',...
    colors(i),'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
    colors(i),'MarkerEdgeColor',colors(i),'CapSize',4); hold on;
end

axis([0 60 0 15])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:3:15)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Concentration Ratio\n(A.U.)'));
legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','South',...
    'FontSize',10); legend boxoff;
title(sprintf("Reperfusion 1h\nBile/Whole Blood"));

for i = 1:length(Bi_to_Bl_Files_4H)
    currData = load(Bi_to_Bl_Files_4H(i));
    tData = currData(:,1);
    Data = currData(:,2:end);
    
    subplot(2,2,2);
    errorbar(tData(2:end),mean(Data(2:end,:),2,'omitnan'),...
    std(Data(2:end,:),0,2,'omitnan')/sqrt(length(Data(1,:))),shapes(i),'Color',...
    colors(i),'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
    colors(i),'MarkerEdgeColor',colors(i),'CapSize',4); hold on;
end

axis([0 60 0 15])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:3:15)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Concentration Ratio\n(A.U.)'));
title(sprintf("Reperfusion 4h\nBile/Whole Blood"));

for i = 1:length(Bi_to_Pl_Files_60)
    currData = load(Bi_to_Pl_Files_60(i));
    tData = currData(:,1);
    Data = currData(:,2:end);
    
    subplot(2,2,3);
    errorbar(tData(2:end),mean(Data(2:end,:),2,'omitnan'),...
    std(Data(2:end,:),0,2,'omitnan')/sqrt(length(Data(1,:))),shapes(i),'Color',...
    colors(i),'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
    colors(i),'MarkerEdgeColor',colors(i),'CapSize',4); hold on;
end

axis([0 60 0 15])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:3:15)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Concentration Ratio\n(A.U.)'));
legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','South',...
    'FontSize',10); legend boxoff;
title(sprintf("Bile/Plasma"));

for i = 1:length(Bi_to_Pl_Files_4H)
    currData = load(Bi_to_Pl_Files_4H(i));
    tData = currData(:,1);
    Data = currData(:,2:end);
    
    subplot(2,2,4);
    errorbar(tData(2:end),mean(Data(2:end,:),2,'omitnan'),...
    std(Data(2:end,:),0,2,'omitnan')/sqrt(length(Data(1,:))),shapes(i),'Color',...
    colors(i),'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
    colors(i),'MarkerEdgeColor',colors(i),'CapSize',4); hold on;
end

axis([0 60 0 15])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:3:15)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Concentration Ratio\n(A.U.)'));
title(sprintf("Bile/Plasma"));

cd ..;

cd Concentration_Data_Avg;

list = dir;
fileNames = [];
for i = 3:length(list) %ignores . and .. and mean fluorescence folder
    fileNames = [fileNames,convertCharsToStrings(list(i).name)];
end

figure (2);
set(gcf,'Units','inches','Position',[0.5 0.5 10 8]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 10 8],'color','white');
shapes = ["-^","-s","-o"];
colors = ["g","b","r"];

for i = 1:(length(fileNames)/3)
    BiData = load(fileNames(3*i-2));
    BlData = load(fileNames(3*i-1));
    PlData = load(fileNames(3*i));
    
    Bi2Bl = BiData(2:end,2)./BlData(2:end,2);
    Bi2Pl = BiData(2:end,2)./PlData(2:end,2);
    
    Bi2BlErr = ((BiData(2:end,4)./BiData(2:end,2)).^2 + (BlData(2:end,4)./BlData(2:end,2)).^2).^(1/2);
    Bi2PlErr = ((BiData(2:end,4)./BiData(2:end,2)).^2 + (PlData(2:end,4)./PlData(2:end,2)).^2).^(1/2);
    
    % Plotting
    
    if i == 1
        %Ctrl
        subplot(2,2,1)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-^','Color','g','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'g','MarkerEdgeColor','g','CapSize',8); hold on;
        subplot(2,2,2)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-^','Color','g','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'g','MarkerEdgeColor','g','CapSize',8); hold on;
        subplot(2,2,3)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-^','Color','g','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'g','MarkerEdgeColor','g','CapSize',8); hold on;
        subplot(2,2,4)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-^','Color','g','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'g','MarkerEdgeColor','g','CapSize',8); hold on;
    elseif i == 2
        %I30R4H
        subplot(2,2,2)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-s','Color','b','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'b','MarkerEdgeColor','b','CapSize',8); hold on;
        subplot(2,2,4)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-s','Color','b','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'b','MarkerEdgeColor','b','CapSize',8); hold on;  
    elseif i == 3
        %I30R60
        subplot(2,2,1)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-s','Color','b','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'b','MarkerEdgeColor','b','CapSize',8); hold on;
        subplot(2,2,3)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-s','Color','b','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'b','MarkerEdgeColor','b','CapSize',8); hold on;
    elseif i == 4
        %I60R4H
        subplot(2,2,2)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-o','Color','r','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'r','MarkerEdgeColor','r','CapSize',8); hold on;
        subplot(2,2,4)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-o','Color','r','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'r','MarkerEdgeColor','r','CapSize',8); hold on;
    else
        %I60R60
        subplot(2,2,1)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-o','Color','r','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'r','MarkerEdgeColor','r','CapSize',8); hold on;
        subplot(2,2,3)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-o','Color','r','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'r','MarkerEdgeColor','r','CapSize',8); hold on;
    end
end
subplot(2,2,1)
axis([0 60 0 15])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:3:15)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Concentration Ratio\n(A.U.)'));
legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','South',...
    'FontSize',10); legend boxoff;
title(sprintf("Reperfusion 1h\nBile/Whole Blood"));

subplot(2,2,2)
axis([0 60 0 15])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:3:15)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Concentration Ratio\n(A.U.)'));
title(sprintf("Reperfusion 4h\nBile/Whole Blood"));

subplot(2,2,3)
axis([0 60 0 15])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:3:15)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Concentration Ratio\n(A.U.)'));
legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','South',...
    'FontSize',10); legend boxoff;
title(sprintf("Bile/Plasma"));

subplot(2,2,4)
axis([0 60 0 15])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:3:15)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Concentration Ratio\n(A.U.)'));
title(sprintf("Bile/Plasma"));

cd ..;

% Load in data
cd 'Fluorescence_Data_Indiv';

list = dir;

fileNames = [];
for i = 3:length(list) %ignores . and .. and mean fluorescence folder
    fileNames = [fileNames,convertCharsToStrings(list(i).name)];
end

cd ..;

if(exist('Fl_Ratio_Data','file') == 0)
    mkdir('Fl_Ratio_Data') ;
end

%% Obtain ratios in Figure 4
for i = 1:(length(fileNames)/3)
    % Get into folder containing Concentration_Data
    cd Fluorescence_Data_Indiv;

    % Loads data and obtains which experimental condition (e.g., control,
    % I30R4H, etc.)
    BileData = load(fileNames(3*i-2));
    BileData = BileData(:,2:end);
    
    BloodData = load(fileNames(3*i-1));
    tData = BloodData(:,1);
    BloodData = BloodData(:,2:end);
    cond = split(fileNames(3*i-1),"_Blood");
    cond = cond(1);
    
    PlasmaData = load(fileNames(3*i));
    PlasmaData = PlasmaData(:,2:end);

    % Take ratios shown in Figure 4. Results in indicated ratios for
    % individual biological replicates in an experimental group. 
    Bi_Bl = BileData./BloodData;
    Bi_Pl = BileData./PlasmaData;

    % Back out of Concentration_Data folder into super folder
    cd ..;
    % Get into folder where things should be saved
    cd 'Fl_Ratio_Data';
    
    Bi_Bl = [tData, Bi_Bl];
    Bi_Pl = [tData, Bi_Pl];

    writematrix(Bi_Bl,strcat(cond,"_Bile_to_Blood_Ratio_Indiv.txt"));
    writematrix(Bi_Pl,strcat(cond,"_Bile_to_Plasma_Ratio_Indiv.txt"));

    % Back out of Concentration_Data folder into super folder
    cd ..;
end

cd 'Fl_Ratio_Data';

list = dir;

fileNames = [];
for i = 3:length(list) %ignores . and .. and mean fluorescence folder
    fileNames = [fileNames,convertCharsToStrings(list(i).name)];
end

sort(fileNames);

Bi_to_Bl_Files_60 = [];
Bi_to_Bl_Files_4H = [];
Bi_to_Pl_Files_60 = [];
Bi_to_Pl_Files_4H = [];
for i = 1:(length(fileNames))
    if contains(fileNames(i),"Control") == 1
        if contains(fileNames(i),"Bl")
            Bi_to_Bl_Files_60 = [Bi_to_Bl_Files_60 fileNames(i)];
            Bi_to_Bl_Files_4H = [Bi_to_Bl_Files_4H fileNames(i)];
        else
            Bi_to_Pl_Files_60 = [Bi_to_Pl_Files_60 fileNames(i)];
            Bi_to_Pl_Files_4H = [Bi_to_Pl_Files_4H fileNames(i)];
        end
        continue;
    end

    if contains(fileNames(i),"Bl") == 1
        if contains(fileNames(i),"R60") == 1
            Bi_to_Bl_Files_60 = [Bi_to_Bl_Files_60 fileNames(i)];
        else
            Bi_to_Bl_Files_4H = [Bi_to_Bl_Files_4H fileNames(i)];
        end
    else
        if contains(fileNames(i),"R60") == 1
            Bi_to_Pl_Files_60 = [Bi_to_Pl_Files_60 fileNames(i)];
        else
            Bi_to_Pl_Files_4H = [Bi_to_Pl_Files_4H fileNames(i)];
        end
    end
end

figure(3);
set(gcf,'Units','inches','Position',[0.5 0.5 10 8]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 10 8],'color','white');
shapes = ["-^","-s","-o"];
colors = ["g","b","r"];

for i = 1:length(Bi_to_Bl_Files_60)
    currData = load(Bi_to_Bl_Files_60(i));
    tData = currData(:,1);
    Data = currData(:,2:end);
    
    subplot(2,2,1);
    errorbar(tData(2:end),mean(Data(2:end,:),2,'omitnan'),...
    std(Data(2:end,:),0,2,'omitnan')/sqrt(length(Data(1,:))),shapes(i),'Color',...
    colors(i),'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
    colors(i),'MarkerEdgeColor',colors(i),'CapSize',4); hold on;
end

axis([0 60 0 10])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:2:10)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Fluorescence Ratio\n(A.U.)'));
legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','South',...
    'FontSize',10); legend boxoff;
title(sprintf("Reperfusion 1h\nBile/Whole Blood"));

for i = 1:length(Bi_to_Bl_Files_4H)
    currData = load(Bi_to_Bl_Files_4H(i));
    tData = currData(:,1);
    Data = currData(:,2:end);
    
    subplot(2,2,2);
    errorbar(tData(2:end),mean(Data(2:end,:),2,'omitnan'),...
    std(Data(2:end,:),0,2,'omitnan')/sqrt(length(Data(1,:))),shapes(i),'Color',...
    colors(i),'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
    colors(i),'MarkerEdgeColor',colors(i),'CapSize',4); hold on;
end

axis([0 60 0 10])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:2:10)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Fluorescence Ratio\n(A.U.)'));
title(sprintf("Reperfusion 4h\nBile/Whole Blood"));

for i = 1:length(Bi_to_Pl_Files_60)
    currData = load(Bi_to_Pl_Files_60(i));
    tData = currData(:,1);
    Data = currData(:,2:end);
    
    subplot(2,2,3);
    errorbar(tData(2:end),mean(Data(2:end,:),2,'omitnan'),...
    std(Data(2:end,:),0,2,'omitnan')/sqrt(length(Data(1,:))),shapes(i),'Color',...
    colors(i),'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
    colors(i),'MarkerEdgeColor',colors(i),'CapSize',4); hold on;
end

axis([0 60 0 4])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:1:4)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Fluorescence Ratio\n(A.U.)'));
legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','North',...
    'FontSize',10); legend boxoff;
title(sprintf("Bile/Plasma"));

for i = 1:length(Bi_to_Pl_Files_4H)
    currData = load(Bi_to_Pl_Files_4H(i));
    tData = currData(:,1);
    Data = currData(:,2:end);
    
    subplot(2,2,4);
    errorbar(tData(2:end),mean(Data(2:end,:),2,'omitnan'),...
    std(Data(2:end,:),0,2,'omitnan')/sqrt(length(Data(1,:))),shapes(i),'Color',...
    colors(i),'LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
    colors(i),'MarkerEdgeColor',colors(i),'CapSize',4); hold on;
end

axis([0 60 0 4])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:1:4));
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Fluorescence Ratio\n(A.U.)'));
title(sprintf("Bile/Plasma"));

cd ..;

cd Fluorescence_Data_Avg;

list = dir;
fileNames = [];
for i = 3:length(list) %ignores . and .. and mean fluorescence folder
    fileNames = [fileNames,convertCharsToStrings(list(i).name)];
end

figure(4);
set(gcf,'Units','inches','Position',[0.5 0.5 10 8]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 10 8],'color','white');
shapes = ["-^","-s","-o"];
colors = ["g","b","r"];

for i = 1:(length(fileNames)/3)
    BiData = load(fileNames(3*i-2));
    BlData = load(fileNames(3*i-1));
    PlData = load(fileNames(3*i));
    
    Bi2Bl = BiData(2:end,2)./BlData(2:end,2);
    Bi2Pl = BiData(2:end,2)./PlData(2:end,2);
    
    Bi2BlErr = ((BiData(2:end,4)./BiData(2:end,2)).^2 + (BlData(2:end,4)./BlData(2:end,2)).^2).^(1/2);
    Bi2PlErr = ((BiData(2:end,4)./BiData(2:end,2)).^2 + (PlData(2:end,4)./PlData(2:end,2)).^2).^(1/2);
    
    % Plotting
    
    if i == 1
        %Ctrl
        subplot(2,2,1)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-^','Color','g','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'g','MarkerEdgeColor','g','CapSize',8); hold on;
        subplot(2,2,2)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-^','Color','g','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'g','MarkerEdgeColor','g','CapSize',8); hold on;
        subplot(2,2,3)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-^','Color','g','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'g','MarkerEdgeColor','g','CapSize',8); hold on;
        subplot(2,2,4)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-^','Color','g','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'g','MarkerEdgeColor','g','CapSize',8); hold on;
    elseif i == 2
        %I30R4H
        subplot(2,2,2)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-s','Color','b','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'b','MarkerEdgeColor','b','CapSize',8); hold on;
        subplot(2,2,4)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-s','Color','b','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'b','MarkerEdgeColor','b','CapSize',8); hold on;  
    elseif i == 3
        %I30R60
        subplot(2,2,1)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-s','Color','b','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'b','MarkerEdgeColor','b','CapSize',8); hold on;
        subplot(2,2,3)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-s','Color','b','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'b','MarkerEdgeColor','b','CapSize',8); hold on;
    elseif i == 4
        %I60R4H
        subplot(2,2,2)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-o','Color','r','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'r','MarkerEdgeColor','r','CapSize',8); hold on;
        subplot(2,2,4)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-o','Color','r','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'r','MarkerEdgeColor','r','CapSize',8); hold on;
    else
        %I60R60
        subplot(2,2,1)
        errorbar(BiData(2:end,1),Bi2Bl,Bi2BlErr,'-o','Color','r','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'r','MarkerEdgeColor','r','CapSize',8); hold on;
        subplot(2,2,3)
        errorbar(BiData(2:end,1),Bi2Pl,Bi2PlErr,'-o','Color','r','LineWidth',2.0,'MarkerSize',6,'MarkerFaceColor',...
        'r','MarkerEdgeColor','r','CapSize',8); hold on;
    end
end
subplot(2,2,1)
axis([0 60 0 10])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:2:10));
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Fluorescence Ratio\n(A.U.)'));
legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','South',...
    'FontSize',10); legend boxoff;
title(sprintf("Reperfusion 1h\nBile/Whole Blood"));

subplot(2,2,2)
axis([0 60 0 10])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:2:10));
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Fluorescence Ratio\n(A.U.)'));
title(sprintf("Reperfusion 4h\nBile/Whole Blood"));

subplot(2,2,3)
axis([0 60 0 4])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:1:4)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Fluorescence Ratio\n(A.U.)'));
legend('WIT = 0 min','WIT = 30 min','WIT = 60 min','Location','North',...
    'FontSize',10); legend boxoff;
title(sprintf("Bile/Plasma"));

subplot(2,2,4)
axis([0 60 0 4])
set(gca,'XTick',(0:10:60)); set(gca,'YTick',(0:1:4)); 
set(gca,'LineWidth',1.5,'FontSize',16); box off;
xlabel('Time (min)'); ylabel(sprintf('Fluorescence Ratio\n(A.U.)'));
title(sprintf("Bile/Plasma"));

cd ..;