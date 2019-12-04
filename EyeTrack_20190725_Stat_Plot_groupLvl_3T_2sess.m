%% Plot for all participants altogether

%% work log

%   25_07_2019  created the script

%% Main task

%% preparation

clear;close all;
% warning('off')
addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20180613'))

path_gen    = '/Users/yeojin/Desktop/';

path_par    = [ path_gen 'E_data/' ];
path_dat    = [ path_par 'EA_raw/EAA_pupil/pilot_3T/'];
path_behav  = [ path_par 'EA_raw/EAC_behav/pilot_3T_2sess/tmp/' ];
path_clean  = [ path_par 'EB_cleaned/EBB_pupil/pilot_3T/' ];
path_fig    = [ path_gen 'C_writings/CB_figures/pilot_3T_2sess/MainTask/zscores/' ];
path_stim   = [ path_par 'EA_raw/EAA_pupil/stiminfo_all_20181004.mat' ];

% create file list and VP ID list
VPs = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221];
tag = [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 0 2; 1 0; 0 0; 1 2; 1 2; 0 0]; % ID 2113, Day 1 data excluded for bad quality

expdat = [];
for ids = 1:length(VPs)
    for d = 1:2
        if tag(ids,d) == 0
            expdat{ids,d} = {};
            rewcond(ids,d) = 0;
        else
            expdat{ids,d}    = load([path_behav num2str(VPs(ids)) '_' num2str(d) '.mat']);
            rewcond(ids,d) = eval(['str2num(expdat{ids,d}.dat.day' num2str(d) '.RewardCategory(9));']);
        end
    end
end

%% now run

dat = [];
for i2 = 1:length(VPs)
    for d = 1:2
        if tag(i2,d) ~= 0
            dat{i2,d} = load([ path_clean 're' num2str(d) num2str(VPs(i2)) '_preproc_dat.mat']);
            dat_raw{i2,d} = load([path_clean 're' num2str(d) num2str(VPs(i2)) '_eye_rv200_stim.mat'],'eye_rv');
            accuracies{i2,d} = dat_raw{i2,d}.eye_rv.accuracy;
            rew{i2,d} = cell2mat(dat_raw{i2,d}.eye_rv.trlord(:,2))==rewcond(i2,d);
            
        else
            dat{i2,d} = {NaN};
            dat_raw{i2,d} = {NaN};
            accuracies{i2,d} = {NaN};
            rew{i2,d} = {NaN};
        end
    end
end


% everything collapsed
sR=[];sN=[];fR=[];fN=[];
for i3 = 1:size(dat,1) % number of IDs
    for d = 1:2
        if tag(i3,d) == 0
        else
            if isempty(sR)
                sR = dat{i3,d}.plotdat_preproc.stim.rew;
                sN = dat{i3,d}.plotdat_preproc.stim.neu;
                fR = dat{i3,d}.plotdat_preproc.fb.rew;
                fN = dat{i3,d}.plotdat_preproc.fb.neu;
                sR_corr = dat{i3,d}.plotdat_preproc.stim.rew(accuracies{i3,d}(rew{i3,d}==1)==1,:);
                sN_corr = dat{i3,d}.plotdat_preproc.stim.neu(accuracies{i3,d}(rew{i3,d}~=1)==1,:);
                
                %         rew_h  = dat{i3,1}.plotdat_preproc.stim.hiMem_rew;
                %         pun_h = dat{i3,1}.plotdat_preproc.stim.hiMem_pun;
                % %         fR_h  = dat{i3,1}.plotdat_preproc.fb.R_mH;
                % %         fnR_h = dat{i3,1}.plotdat_preproc.fb.nR_mH;
                %
                %         rew_l  = dat{i3,1}.plotdat_preproc.stim.loMem_rew;
                %         pun_l = dat{i3,1}.plotdat_preproc.stim.loMem_pun;
                % %         fR_l  = dat{i3,1}.plotdat_preproc.fb.R_mL;
                % %         fnR_l = dat{i3,1}.plotdat_preproc.fb.nR_mL;
                %
                %         highmem = dat{i3,1}.plotdat_preproc.stim.hiMem;
                %         lowmem  = dat{i3,1}.plotdat_preproc.stim.loMem;
                %
                %         smH{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mH,dat{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mL,dat{i3,1}.plotdat_preproc.stim.R_mL);
                
            else
                sR  = vertcat(sR,dat{i3,d}.plotdat_preproc.stim.rew);
                sN = vertcat(sN,dat{i3,d}.plotdat_preproc.stim.neu);
                fR  = vertcat(fR,dat{i3,d}.plotdat_preproc.fb.rew);
                fN = vertcat(fN,dat{i3,d}.plotdat_preproc.fb.neu);
                sR_corr = vertcat(sR_corr,dat{i3,d}.plotdat_preproc.stim.rew(accuracies{i3,d}(rew{i3,d}==1)==1,:));
                sN_corr = vertcat(sN_corr,dat{i3,d}.plotdat_preproc.stim.neu(accuracies{i3,d}(rew{i3,d}~=1)==1,:));
                
                %         rew_h  = vertcat(rew_h,dat{i3,1}.plotdat_preproc.stim.hiMem_rew);
                %         pun_h = vertcat(pun_h,dat{i3,1}.plotdat_preproc.stim.hiMem_pun);
                % %         fR_h  = vertcat(fR_h,dat{i3,1}.plotdat_preproc.fb.R_mH);
                % %         fnR_h = vertcat(fnR_h,dat{i3,1}.plotdat_preproc.fb.nR_mH);
                %
                %         rew_l  = vertcat(rew_l,dat{i3,1}.plotdat_preproc.stim.loMem_rew);
                %         pun_l = vertcat(pun_l,dat{i3,1}.plotdat_preproc.stim.loMem_pun);
                % %         fR_l  = vertcat(fR_l,dat{i3,1}.plotdat_preproc.fb.R_mL);
                % %         fnR_l = vertcat(fnR_l,dat{i3,1}.plotdat_preproc.fb.nR_mL);
                %
                %         highmem = vertcat(highmem,dat{i3,1}.plotdat_preproc.stim.hiMem);
                %         lowmem  = vertcat(lowmem,dat{i3,1}.plotdat_preproc.stim.loMem);
                %
                %         smH{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mH,dat{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mL,dat{i3,1}.plotdat_preproc.stim.R_mL);
                
                
            end
        end
    end
    
end

sR_all = sR; sN_all = sN; fR_all = fR; fN_all = fN;
sR_corr_all = sR_corr; sN_corr_all = sN_corr; 

% only session1 collapsed
sR=[];sN=[];fR=[];fN=[];
for i3 = 1:size(dat,1) % number of IDs
    for d = 1
        if tag(i3,d) == 0
        else
            if isempty(sR)
                sR = dat{i3,d}.plotdat_preproc.stim.rew;
                sN = dat{i3,d}.plotdat_preproc.stim.neu;
                fR = dat{i3,d}.plotdat_preproc.fb.rew;
                fN = dat{i3,d}.plotdat_preproc.fb.neu;
                sR_corr = dat{i3,d}.plotdat_preproc.stim.rew(accuracies{i3,d}(rew{i3,d}==1)==1,:);
                sN_corr = dat{i3,d}.plotdat_preproc.stim.neu(accuracies{i3,d}(rew{i3,d}~=1)==1,:);
                
                %         rew_h  = dat{i3,1}.plotdat_preproc.stim.hiMem_rew;
                %         pun_h = dat{i3,1}.plotdat_preproc.stim.hiMem_pun;
                % %         fR_h  = dat{i3,1}.plotdat_preproc.fb.R_mH;
                % %         fnR_h = dat{i3,1}.plotdat_preproc.fb.nR_mH;
                %
                %         rew_l  = dat{i3,1}.plotdat_preproc.stim.loMem_rew;
                %         pun_l = dat{i3,1}.plotdat_preproc.stim.loMem_pun;
                % %         fR_l  = dat{i3,1}.plotdat_preproc.fb.R_mL;
                % %         fnR_l = dat{i3,1}.plotdat_preproc.fb.nR_mL;
                %
                %         highmem = dat{i3,1}.plotdat_preproc.stim.hiMem;
                %         lowmem  = dat{i3,1}.plotdat_preproc.stim.loMem;
                %
                %         smH{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mH,dat{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mL,dat{i3,1}.plotdat_preproc.stim.R_mL);
                
            else
                sR  = vertcat(sR,dat{i3,d}.plotdat_preproc.stim.rew);
                sN = vertcat(sN,dat{i3,d}.plotdat_preproc.stim.neu);
                fR  = vertcat(fR,dat{i3,d}.plotdat_preproc.fb.rew);
                fN = vertcat(fN,dat{i3,d}.plotdat_preproc.fb.neu);
                sR_corr = vertcat(sR_corr,dat{i3,d}.plotdat_preproc.stim.rew(accuracies{i3,d}(rew{i3,d}==1)==1,:));
                sN_corr = vertcat(sN_corr,dat{i3,d}.plotdat_preproc.stim.neu(accuracies{i3,d}(rew{i3,d}~=1)==1,:));
                
                %         rew_h  = vertcat(rew_h,dat{i3,1}.plotdat_preproc.stim.hiMem_rew);
                %         pun_h = vertcat(pun_h,dat{i3,1}.plotdat_preproc.stim.hiMem_pun);
                % %         fR_h  = vertcat(fR_h,dat{i3,1}.plotdat_preproc.fb.R_mH);
                % %         fnR_h = vertcat(fnR_h,dat{i3,1}.plotdat_preproc.fb.nR_mH);
                %
                %         rew_l  = vertcat(rew_l,dat{i3,1}.plotdat_preproc.stim.loMem_rew);
                %         pun_l = vertcat(pun_l,dat{i3,1}.plotdat_preproc.stim.loMem_pun);
                % %         fR_l  = vertcat(fR_l,dat{i3,1}.plotdat_preproc.fb.R_mL);
                % %         fnR_l = vertcat(fnR_l,dat{i3,1}.plotdat_preproc.fb.nR_mL);
                %
                %         highmem = vertcat(highmem,dat{i3,1}.plotdat_preproc.stim.hiMem);
                %         lowmem  = vertcat(lowmem,dat{i3,1}.plotdat_preproc.stim.loMem);
                %
                %         smH{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mH,dat{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mL,dat{i3,1}.plotdat_preproc.stim.R_mL);
                
                
            end
        end
    end
    
end

sR_1 = sR; sN_1 = sN; fR_1 = fR; fN_1 = fN;

% session 2 collapsed
sR=[];sN=[];fR=[];fN=[];
for i3 = 1:size(dat,1) % number of IDs
    for d = 2
        if tag(i3,d) == 0
        else
            if isempty(sR)
                sR = dat{i3,d}.plotdat_preproc.stim.rew;
                sN = dat{i3,d}.plotdat_preproc.stim.neu;
                fR = dat{i3,d}.plotdat_preproc.fb.rew;
                fN = dat{i3,d}.plotdat_preproc.fb.neu;
                sR_corr = vertcat(sR_corr,dat{i3,d}.plotdat_preproc.stim.rew(accuracies{i3,d}(rew{i3,d}==1)==1,:));
                sN_corr = vertcat(sN_corr,dat{i3,d}.plotdat_preproc.stim.neu(accuracies{i3,d}(rew{i3,d}~=1)==1,:));
                
                %         rew_h  = dat{i3,1}.plotdat_preproc.stim.hiMem_rew;
                %         pun_h = dat{i3,1}.plotdat_preproc.stim.hiMem_pun;
                % %         fR_h  = dat{i3,1}.plotdat_preproc.fb.R_mH;
                % %         fnR_h = dat{i3,1}.plotdat_preproc.fb.nR_mH;
                %
                %         rew_l  = dat{i3,1}.plotdat_preproc.stim.loMem_rew;
                %         pun_l = dat{i3,1}.plotdat_preproc.stim.loMem_pun;
                % %         fR_l  = dat{i3,1}.plotdat_preproc.fb.R_mL;
                % %         fnR_l = dat{i3,1}.plotdat_preproc.fb.nR_mL;
                %
                %         highmem = dat{i3,1}.plotdat_preproc.stim.hiMem;
                %         lowmem  = dat{i3,1}.plotdat_preproc.stim.loMem;
                %
                %         smH{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mH,dat{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mL,dat{i3,1}.plotdat_preproc.stim.R_mL);
                
            else
                sR  = vertcat(sR,dat{i3,d}.plotdat_preproc.stim.rew);
                sN = vertcat(sN,dat{i3,d}.plotdat_preproc.stim.neu);
                fR  = vertcat(fR,dat{i3,d}.plotdat_preproc.fb.rew);
                fN = vertcat(fN,dat{i3,d}.plotdat_preproc.fb.neu);
                sR_corr = vertcat(sR_corr,dat{i3,d}.plotdat_preproc.stim.rew(accuracies{i3,d}(rew{i3,d}==1)==1,:));
                sN_corr = vertcat(sN_corr,dat{i3,d}.plotdat_preproc.stim.neu(accuracies{i3,d}(rew{i3,d}~=1)==1,:));
                
                %         rew_h  = vertcat(rew_h,dat{i3,1}.plotdat_preproc.stim.hiMem_rew);
                %         pun_h = vertcat(pun_h,dat{i3,1}.plotdat_preproc.stim.hiMem_pun);
                % %         fR_h  = vertcat(fR_h,dat{i3,1}.plotdat_preproc.fb.R_mH);
                % %         fnR_h = vertcat(fnR_h,dat{i3,1}.plotdat_preproc.fb.nR_mH);
                %
                %         rew_l  = vertcat(rew_l,dat{i3,1}.plotdat_preproc.stim.loMem_rew);
                %         pun_l = vertcat(pun_l,dat{i3,1}.plotdat_preproc.stim.loMem_pun);
                % %         fR_l  = vertcat(fR_l,dat{i3,1}.plotdat_preproc.fb.R_mL);
                % %         fnR_l = vertcat(fnR_l,dat{i3,1}.plotdat_preproc.fb.nR_mL);
                %
                %         highmem = vertcat(highmem,dat{i3,1}.plotdat_preproc.stim.hiMem);
                %         lowmem  = vertcat(lowmem,dat{i3,1}.plotdat_preproc.stim.loMem);
                %
                %         smH{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mH,dat{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(dat{i3,1}.plotdat_preproc.stim.nR_mL,dat{i3,1}.plotdat_preproc.stim.R_mL);
                
                
            end
        end
    end
    
end

sR_2 = sR; sN_2 = sN; fR_2 = fR; fN_2 = fN;


cd(path_fig) % saving path for the figures

% fancier: shadedErrorBar

brightred = [1 , 0 , 0]; % red is first condition
brightblue = [0 , 0 , 1];
darkred = [0.5 , 0 , 0]; 
darkblue = [0 , 0 , 0.5];

% plot with expectancy: Rewards
% stim
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_2,{@nanmean,errorbar},'lineprops',{'-','LineWidth',10,'color',brightred,'markerfacecolor',brightred}); hold on; % rare reward vs frequent reward
shadedErrorBar((1:2701),sR_1,{@nanmean,errorbar},'lineprops',{'--','LineWidth',10,'color',brightred,'markerfacecolor',brightred}); hold on;
H = ttest2(sR_2,sR_1,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('stimulus','Fontsize',20,'Fontweight','bold')
clear H sigbar

% fb
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),fR_2,{@nanmean,errorbar},'lineprops',{'-','LineWidth',10,'color',brightred,'markerfacecolor',brightred}); hold on; % rare reward vs frequent reward
shadedErrorBar((1:2701),fR_1,{@nanmean,errorbar},'lineprops',{'--','LineWidth',10,'color',brightred,'markerfacecolor',brightred}); hold on;
H = ttest2(fR_2,fR_1,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('fb','Fontsize',20,'Fontweight','bold')
clear H sigbar

suptitle(['Rewards: rare(red) vs. frequent(blue)'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';

saveas(gcf,['GroupLvl_3T_2sess_maintask_Rewards_rare(solid)_vs_frequent(dashed).fig'])
print(['GroupLvl_3T_2sess_maintask_Rewards_rare(solid)_vs_frequent(dashed)'],'-dpdf','-r0')
clear fig H sigbar



% plot with expectancy: Neutrals
% stim
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sN_1,{@nanmean,errorbar},'lineprops',{'-','LineWidth',10,'color',brightblue,'markerfacecolor',brightblue}); hold on; % rare reward vs frequent reward
shadedErrorBar((1:2701),sN_2,{@nanmean,errorbar},'lineprops',{'--','LineWidth',10,'color',brightblue,'markerfacecolor',brightblue}); hold on;
H = ttest2(sN_1,sN_2,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('stimulus','Fontsize',20,'Fontweight','bold')
clear H sigbar

% fb
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),fN_1,{@nanmean,errorbar},'lineprops',{'-','LineWidth',10,'color',brightblue,'markerfacecolor',brightblue}); hold on; % rare reward vs frequent reward
shadedErrorBar((1:2701),fN_2,{@nanmean,errorbar},'lineprops',{'--','LineWidth',10,'color',brightblue,'markerfacecolor',brightblue}); hold on;
H = ttest2(fN_1,fN_2,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('fb','Fontsize',20,'Fontweight','bold')
clear H sigbar

suptitle(['Neutrals: rare(solid) vs. frequent(dashed)'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';

saveas(gcf,['GroupLvl_3T_2sess_maintask_Neutrals_rare(solid)_vs_frequent(dashed).fig'])
print(['GroupLvl_3T_2sess_maintask_Neutrals_rare(solid)_vs_frequent(dashed)'],'-dpdf','-r0')
clear fig H sigbar



% plot with expectancy, only stimulus segment
maxy = 0.3; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
plot(onsetbar,y,'g-','linewidth',2)
shadedErrorBar((1:2701),sN_1,{@nanmean,errorbar},'lineprops',{'-','LineWidth',7,'color',brightblue,'markerfacecolor',brightblue}); hold on; % rare N
shadedErrorBar((1:2701),sN_2,{@nanmean,errorbar},'lineprops',{'--','LineWidth',7,'color',brightblue,'markerfacecolor',brightblue}); hold on; % frequent N
shadedErrorBar((1:2701),sR_2,{@nanmean,errorbar},'lineprops',{'-','LineWidth',7,'color',brightred,'markerfacecolor',brightred}); hold on; % rare R
shadedErrorBar((1:2701),sR_1,{@nanmean,errorbar},'lineprops',{'--','LineWidth',7,'color',brightred,'markerfacecolor',brightred}); hold on; % frequent R
H = ttest2(sN_1,sN_2,'Alpha',0.05); sigbar = H*-0.2; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'c-','LineWidth',7); clear H sigbar % sigbar for rare vs freq neutral
H = ttest2(sR_2,sR_1,'Alpha',0.05); sigbar = H*-0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'m-','LineWidth',7); clear H sigbar % sigbar for rare vs freq reward
H = ttest2(sR_2,sN_1,'Alpha',0.05); sigbar = H*0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-','LineWidth',7); clear H sigbar % sigbar for rare reward vs neutral
H = ttest2(sR_1,sN_2,'Alpha',0.05); sigbar = H*0.2; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-.','LineWidth',7); clear H sigbar % sigbar for freq reward vs neutral

xlim([1 2701]); ylim([miny maxy]);

suptitle(['rare(solid) vs. frequent(dashed) & reward(red) vs. neutral(blue)'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';

saveas(gcf,['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_rare(solid)_vs_freq(dashed).fig'])
print(['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_rare(solid)_vs_freq(dashed)'],'-dpdf','-r0')
clear fig H sigbar




% all sessions
% stim
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_all,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_all,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_all,sN_all,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('stimulus','Fontsize',20,'Fontweight','bold')
clear H sigbar

% feedback
subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),fR_all,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),fN_all,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(fR_all,fN_all,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('feedback','Fontsize',20,'Fontweight','bold')

suptitle(['Reward(red) vs. neutral(blue): Session 1 and 2'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';

saveas(gcf,['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_all_sessions.fig'])
print(['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_all_sessions'],'-dpdf','-r0')
clear fig H sigbar



% all sessions, correct trials
% stim
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure; plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_corr_all,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_corr_all,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_corr_all,sN_corr_all,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('stimulus','Fontsize',20,'Fontweight','bold')
clear H sigbar

% save
fig = gcf;
fig.PaperPositionMode = 'auto';

saveas(gcf,['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_all_sessions_corrects.fig'])
print(['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_all_sessions_corrects'],'-dpdf','-r0')
clear fig H sigbar


% each sessions

for d = 1:2

% stim
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

eval(['shadedErrorBar((1:2701),sR_' num2str(d) ',{@nanmean,errorbar},''lineprops'',{''r-o'',''markerfacecolor'',''r''})']); hold on;
eval(['shadedErrorBar((1:2701),sN_' num2str(d) ',{@nanmean,errorbar},''lineprops'',{''b-o'',''markerfacecolor'',''b''})']); hold on;
eval(['H = ttest2(sR_' num2str(d) ',sN_' num2str(d) ');'])
sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('stimulus','Fontsize',20,'Fontweight','bold')
clear H sigbar

% feedback
subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

eval(['shadedErrorBar((1:2701),fR_' num2str(d) ',{@nanmean,errorbar},''lineprops'',{''r-o'',''markerfacecolor'',''r''})']); hold on;
eval(['shadedErrorBar((1:2701),fN_' num2str(d) ',{@nanmean,errorbar},''lineprops'',{''b-o'',''markerfacecolor'',''b''})']); hold on;
eval(['H = ttest2(fR_' num2str(d) ',fN_' num2str(d) ');'])
 sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('feedback','Fontsize',20,'Fontweight','bold')

suptitle(['Reward(red) vs. neutral(blue): Session ' num2str(d)])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';

saveas(gcf,['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_session' num2str(d) '.fig'])
print(['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_session' num2str(d)],'-dpdf','-r0')
clear fig H sigbar


end



% session 1 & 2
% stim_1
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(2,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_1,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_1,sN_1,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('session1 stim','Fontsize',20,'Fontweight','bold')
clear H sigbar

% feedback_1
subplot(2,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),fR_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),fN_1,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(fR_1,fN_1,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('session1 feedback','Fontsize',20,'Fontweight','bold')

% stim_2
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

subplot(2,2,3);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_2,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_2,sN_2,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('session2 stim','Fontsize',20,'Fontweight','bold')
clear H sigbar

% feedback_2
subplot(2,2,4);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),fR_2,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),fN_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(fR_2,fN_2,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('session2 feedback','Fontsize',20,'Fontweight','bold')


suptitle(['Reward(red) vs. neutral(blue): Session 1 and 2'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
saveas(gcf,'GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_all_sessions-cross.fig')
print(['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_all_sessions-cross'],'-dpdf','-r0')
clear fig H sigbar



% session 1 vs 2
% stim_reward
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(2,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sR_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_1,sR_2,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('Stim & Reward','Fontsize',17,'Fontweight','bold')
clear H sigbar

% fb_reward
subplot(2,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),fR_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),fR_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(fR_1,fR_2,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('Feedback & Reward','Fontsize',17,'Fontweight','bold')

% stim_neutral
maxy = 1; miny = -1;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

subplot(2,2,3);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sN_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sN_1,sN_2,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('Stim & Neutral','Fontsize',17,'Fontweight','bold')
clear H sigbar

% feedback_2
subplot(2,2,4);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),fN_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),fN_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(fN_1,fN_2,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('Feedback & Neutral','Fontsize',17,'Fontweight','bold')

suptitle(['Session 1(red) vs. Session 2(blue)'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
saveas(gcf,'GroupLvl_3T_2sess_maintask_Session1(red)_vs_Session2(blue)_all_sessions_compared.fig')
print(['GroupLvl_3T_2sess_maintask_Session1(red)_vs_Session2(blue)_all_sessions_compared'],'-dpdf','-r0')
clear fig H sigbar




%
% % reward: high mem vs. low mem
% figure;
% % subplot(1,2,1);
% plot(onsetbar,y,'g-','linewidth',2)
%
% shadedErrorBar((1:2701),rew_h,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),rew_l,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(rew_h,rew_l,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('Reward(himem-red, lomem-blue)','Fontsize',20,'Fontweight','bold')
% clear H sigbar
% % feedback seg: reward vs. non-reward
% % subplot(1,2,2);
% % plot(onsetbar,y,'g-','linewidth',2)
% %
% % shadedErrorBar((1:2701),fR_h,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% % shadedErrorBar((1:2701),fR_l,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% % H = ttest2(fR_h,fR_l,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% % xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% % xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% % ylabel('zscores','Fontsize',20,'Fontweight','bold');
% % title('feedback','Fontsize',20,'Fontweight','bold')
% % suptitle('Rewarded: high Mem vs. low Mem')
%
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_3T_maintask_Reward_highMem_vs_lowMem'],'-dpdf','-r0')
% clear fig H sigbar
%
%
% % no-reward: high mem vs. low mem
% figure;
% % subplot(1,2,1);
% plot(onsetbar,y,'g-','linewidth',2)
%
% shadedErrorBar((1:2701),pun_h,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),pun_l,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(pun_h,pun_l,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('punishment(himem-red, lomem-blue)','Fontsize',20,'Fontweight','bold')
% clear H sigbar
% % feedback seg: reward vs. non-reward
% % subplot(1,2,2);
% % plot(onsetbar,y,'g-','linewidth',2)
% %
% % shadedErrorBar((1:2701),fnR_h,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% % shadedErrorBar((1:2701),fnR_l,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% % H = ttest2(fnR_h,fnR_l,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% % xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% % xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% % ylabel('zscores','Fontsize',20,'Fontweight','bold');
% % title('feedback','Fontsize',20,'Fontweight','bold')
% % suptitle('Non-Rewarded: high Mem vs. low Mem')
%
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_3T_maintask_Punishment_highMem_vs_lowMem'],'-dpdf','-r0')
% clear fig H sigbar
%
%
% % memorability
% figure;
% % subplot(1,2,1);
% plot(onsetbar,y,'g-','linewidth',2)
%
% shadedErrorBar((1:2701),highmem,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),lowmem,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(highmem,lowmem,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('himem-red, lomem-blue','Fontsize',20,'Fontweight','bold')
% clear H sigbar
%
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_3T_maintask_highMem_vs_lowMem'],'-dpdf','-r0')
% clear fig H sigbar

% 
% 
% % single participants in one plot!
% % stim, reward vs. punishment
% figure; hold on;
% for i4 = 1:length(VPs)
%     maxy = 3; miny = -3;
%     onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
%     errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
%     
%     subplot(ceil(sqrt(length(VPs))),ceil(sqrt(length(VPs))),i4);
%     plot(onsetbar,y,'g-','linewidth',2)
%     
%     shadedErrorBar((1:2701),dat{i4,1}.plotdat_preproc.stim.rew,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:2701),dat{i4,1}.plotdat_preproc.stim.pun,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(dat{i4,1}.plotdat_preproc.stim.rew,dat{i4,1}.plotdat_preproc.stim.pun,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',8,'Fontweight','bold');
%     ylabel('zscores','Fontsize',8,'Fontweight','bold');
%     title(num2str(VPs(i4)))
%     clear H sigbar
%     
% end
% suptitle('stimulus seg: Reward(red) vs. Punishment(blue)')
% 
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_3T_maintask_Indiv_stim_Reward(red)_vs_Punishment(blue)'],'-dpdf','-r0')
% clear fig H sigbar

%
% % high memorability vs. low memorability
% close all
% figure; hold on;
% for i4 = 1:length(VPs)
%     maxy = 3; miny = -3;
%     onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
%     errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
%
%     subplot(ceil(sqrt(length(VPs))),ceil(sqrt(length(VPs))),i4);
%     plot(onsetbar,y,'g-','linewidth',2)
%
%     shadedErrorBar((1:2701),smH{i4,1},{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:2701),smL{i4,1},{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(smH{i4,1},smL{i4,1},'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',8,'Fontweight','bold');
%     ylabel('zscores','Fontsize',8,'Fontweight','bold');
%     title(num2str(VPs(i4)))
%     clear H sigbar
%
% end
% suptitle('HighMem(red) vs. LowMem(blue)')
%
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_maintask_pilot' num2str(pilotNr) '_Indiv_highMem(red)_vs_lowMem(blue)'],'-dpdf','-r0')
% clear fig H sigbar
%
%
%
%


% feedback
% close all;
% figure; hold on;
% for i4 = 1:length(VPs)
%     maxy = 3; miny = -3;
%     onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
%     errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
%
%     subplot(ceil(sqrt(length(VPs))),ceil(sqrt(length(VPs))),i4);
%     plot(onsetbar,y,'g-','linewidth',2)
%
%     shadedErrorBar((1:2701),dat{i4,1}.plotdat_preproc.fb.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:2701),dat{i4,1}.plotdat_preproc.fb.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(dat{i4,1}.plotdat_preproc.fb.R,dat{i4,1}.plotdat_preproc.fb.nR,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',8,'Fontweight','bold');
%     ylabel('zscores','Fontsize',8,'Fontweight','bold');
%     title(num2str(VPs(i4)))
%
%     clear H sigbar
% end
% suptitle('feedback seg: Reward(red) vs. Non-Reward(blue)')

% save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_Indiv_fb_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
% clear fig H sigbar


%% Memory test

%% immediate & delayed

%% preparation

clear;clc;close all;
% warning('off')
addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20180613'))

path_gen    = '/Users/yeojin/Desktop/';

path_par    = [ path_gen 'E_data/' ];
path_dat    = [ path_par 'EA_raw/EAA_pupil/pilot_3T/'];
path_behav  = [ path_par 'EA_raw/EAC_behav/pilot_3T/' ];
path_clean  = [ path_par 'EB_cleaned/EBB_pupil/pilot_3T/' ];
path_fig    = [ path_gen 'C_writings/CB_figures/pilot_3T_2sess/MainTask/zscores/' ];
path_stim   = [ path_par 'EA_raw/EAA_pupil/stiminfo_all_20181004.mat' ];

% create file list and VP ID list
VPs = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221]; % 2113(D1), 2207(D2), 2116(D1) 2117(D2) -> bad data / 2201, 2111 -> dropout
tag = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2];
d1m = [1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

%% now run

dat = [];
for i2 = 1:length(VPs)
    for d = 1:2
        if tag(i2,d) ~= 0
            dat{i2,d} = load([ path_clean 'rr' num2str(d) num2str(VPs(i2)) '_preproc_dat.mat']);
            dat_raw_i{i2,d} = load([path_clean 'ri' num2str(d) num2str(VPs(i2)) '_eye_rv200_stim.mat'],'eye_rv');
            accuracies_i{i2,d} = dat_raw_i{i2,d}.eye_rv.immediate.memresults(:,2);
            old_i{i2,d} = cell2mat(dat_raw_i{i2,d}.eye_rv.immediate.memstimvec(:,4))==1;
            try
            dat_raw_d{i2,d} = load([path_clean 'rd' num2str(d) num2str(VPs(i2)) '_eye_rv200_stim.mat'],'eye_rv');
            accuracies_d{i2,d} = dat_raw_d{i2,d}.eye_rv.delayed.memresults(:,2);
            old_d{i2,d} = cell2mat(dat_raw_d{i2,d}.eye_rv.delayed.memstimvec(:,4))==1;
            catch
            dat_raw_d{i2,d} = {NaN};
            accuracies_d{i2,d} = {NaN};
            old_d{i2,d} = {NaN};
            end
        else
            dat{i2,d} = {NaN};
            dat_raw_i{i2,d} = {NaN};
            dat_raw_d{i2,d} = {NaN};
            accuracies_i{i2,d} = {NaN};
            accuracies_d{i2,d} = {NaN};
            old_i{i2,d} = {NaN};
            old_d{i2,d} = {NaN};
        end
    end
end

iR=[];inR=[];dR=[];dnR=[];iO=[];iN=[];iO_corr=[];iN_corr=[];
dO=[];dN=[];dO_corr=[];dN_corr=[];
iR_1=[];inR_1=[];dR_1=[];dnR_1=[];iO_1=[];iN_1=[];iO_corr_1=[];iN_corr_1=[];
dO_1=[];dN_1=[];dO_corr_1=[];dN_corr_1=[];
iR_2=[];inR_2=[];dR_2=[];dnR_2=[];iO_2=[];iN_2=[];iO_corr_2=[];iN_corr_2=[];
dO_2=[];dN_2=[];dO_corr_2=[];dN_corr_2=[];
for i3 = 1:size(dat,1) % number of IDs
    for d = 1:2
        if tag(i3,d) ~= 0
            if isempty(iR)
                if eval(['d' num2str(d) 'm(i3,1) == 1'])
                iR      = dat{i3,d}.plotdat_preproc.imm.rew;
                inR     = dat{i3,d}.plotdat_preproc.imm.neu;
                iO      = dat{i3,d}.plotdat_preproc.imm.old;
                iN      = dat{i3,d}.plotdat_preproc.imm.new;
                iO_corr = dat{i3,d}.plotdat_preproc.imm.old(accuracies_i{i3,d}(old_i{i3,d}==1)==1,:);
                iN_corr = dat{i3,d}.plotdat_preproc.imm.new(accuracies_i{i3,d}(old_i{i3,d}~=1)==1,:);
                
                eval(['iR_' num2str(d) '= dat{i3,d}.plotdat_preproc.imm.rew;'])
                eval(['inR_' num2str(d) '= dat{i3,d}.plotdat_preproc.imm.neu;'])
                eval(['iO_' num2str(d) '= dat{i3,d}.plotdat_preproc.imm.old;'])
                eval(['iN_' num2str(d) '= dat{i3,d}.plotdat_preproc.imm.new;'])
                eval(['iO_corr_' num2str(d) '= dat{i3,d}.plotdat_preproc.imm.old(accuracies_i{i3,d}(old_i{i3,d}==1)==1,:);'])
                eval(['iN_corr_' num2str(d) '= dat{i3,d}.plotdat_preproc.imm.new(accuracies_i{i3,d}(old_i{i3,d}~=1)==1,:);'])
                end
                if eval(['d' num2str(d) 'm(i3,2) == 2'])
                dR      = dat{i3,d}.plotdat_preproc.del.rew;
                dnR     = dat{i3,d}.plotdat_preproc.del.neu;
                dO      = dat{i3,d}.plotdat_preproc.del.old;
                dN      = dat{i3,d}.plotdat_preproc.del.new;
                dO_corr = dat{i3,d}.plotdat_preproc.del.old(accuracies_d{i3,d}(old_d{i3,d}==1)==1,:);
                dN_corr = dat{i3,d}.plotdat_preproc.del.new(accuracies_d{i3,d}(old_d{i3,d}~=1)==1,:);
                
                eval(['dRr_' num2str(d) '= dat{i3,d}.plotdat_preproc.del.rew;'])
                eval(['dnR_' num2str(d) '= dat{i3,d}.plotdat_preproc.del.neu;'])
                eval(['dO_' num2str(d) '= dat{i3,d}.plotdat_preproc.del.old;'])
                eval(['dN_' num2str(d) '= dat{i3,d}.plotdat_preproc.del.new;'])
                eval(['dO_' num2str(d) 'corr = dat{i3,d}.plotdat_preproc.del.old(accuracies_d{i3,d}(old_d{i3,d}==1)==1,:);'])
                eval(['dN_corr_' num2str(d) '= dat{i3,d}.plotdat_preproc.del.new(accuracies_d{i3,d}(old_d{i3,d}~=1)==1,:);'])
                end
                

                %         iRO     = dat{i3,d}.plotdat_preproc.imm.RO;
                %         iRN     = dat{i3,d}.plotdat_preproc.imm.RN;
                %         dRO     = dat{i3,d}.plotdat_preproc.del.RO;
                %         dRN     = dat{i3,d}.plotdat_preproc.del.RN;
                %
                %         inRO    = dat{i3,d}.plotdat_preproc.imm.neuO;
                %         inRN    = dat{i3,d}.plotdat_preproc.imm.neuN;
                %         dnRO    = dat{i3,d}.plotdat_preproc.del.neuO;
                %         dnRN    = dat{i3,d}.plotdat_preproc.del.neuN;
                
                %         % correct trls
                %         iR_c    = dat{i3,1}.plotdat_preproc.imm_c.R;
                %         inR_c   = dat{i3,1}.plotdat_preproc.imm_c.neu;
                %         dR_c    = dat{i3,1}.plotdat_preproc.del_c.R;
                %         dnR_c   = dat{i3,1}.plotdat_preproc.del_c.neu;
                %
                %         iO_c    = dat{i3,1}.plotdat_preproc.imm_c.old;
                %         iN_c    = dat{i3,1}.plotdat_preproc.imm_c.new;
                %         dO_c    = dat{i3,1}.plotdat_preproc.del_c.old;
                %         dN_c    = dat{i3,1}.plotdat_preproc.del_c.new;
                %
                %         iRO_c   = dat{i3,1}.plotdat_preproc.imm_c.RO;
                %         iRN_c   = dat{i3,1}.plotdat_preproc.imm_c.RN;
                %         dRO_c   = dat{i3,1}.plotdat_preproc.del_c.RO;
                %         dRN_c   = dat{i3,1}.plotdat_preproc.del_c.RN;
                %
                %         inRO_c  = dat{i3,1}.plotdat_preproc.imm_c.neuO;
                %         inRN_c  = dat{i3,1}.plotdat_preproc.imm_c.neuN;
                %         dnRO_c  = dat{i3,1}.plotdat_preproc.del_c.neuO;
                %         dnRN_c  = dat{i3,1}.plotdat_preproc.del_c.neuN;
            else
            try
                if eval(['d' num2str(d) 'm(i3,1) == 1'])
                iR  = vertcat(iR,dat{i3,d}.plotdat_preproc.imm.rew);
                inR = vertcat(inR,dat{i3,d}.plotdat_preproc.imm.neu);
                iO    = vertcat(iO,dat{i3,d}.plotdat_preproc.imm.old);
                iN    = vertcat(iN,dat{i3,d}.plotdat_preproc.imm.new);
                iO_corr = vertcat(iO_corr,dat{i3,d}.plotdat_preproc.imm.old(accuracies_i{i3,d}(old_i{i3,d}==1)==1,:));
                iN_corr = vertcat(iN_corr,dat{i3,d}.plotdat_preproc.imm.new(accuracies_i{i3,d}(old_i{i3,d}~=1)==1,:));
                
                eval(['iR_' num2str(d) '= vertcat(iR_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.rew);'])
                eval(['inR_' num2str(d) '= vertcat(inR_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.neu);'])
                eval(['iO_' num2str(d) '= vertcat(iO_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.old);'])
                eval(['iN_' num2str(d) '= vertcat(iN_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.new);'])
                eval(['iO_corr_' num2str(d) '= vertcat(iO_corr_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.old(accuracies_i{i3,d}(old_i{i3,d}==1)==1,:));'])
                eval(['iN_corr_' num2str(d) '= vertcat(iN_corr_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.new(accuracies_i{i3,d}(old_i{i3,d}~=1)==1,:));'])
                end
                if eval(['d' num2str(d) 'm(i3,2) == 2'])
                dR  = vertcat(dR,dat{i3,d}.plotdat_preproc.del.rew);
                dnR = vertcat(dnR,dat{i3,d}.plotdat_preproc.del.neu);
                dO    = vertcat(dO,dat{i3,d}.plotdat_preproc.del.old);
                dN    = vertcat(dN,dat{i3,d}.plotdat_preproc.del.new);
                dO_corr = vertcat(dO_corr,dat{i3,d}.plotdat_preproc.del.old(accuracies_d{i3,d}(old_d{i3,d}==1)==1,:));
                dN_corr = vertcat(dN_corr,dat{i3,d}.plotdat_preproc.del.new(accuracies_d{i3,d}(old_d{i3,d}~=1)==1,:));
                
                eval(['dR_' num2str(d) '= vertcat(dR_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.rew);'])
                eval(['dnR_' num2str(d) '= vertcat(dnR_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.neu);'])
                eval(['dO_' num2str(d) '= vertcat(dO_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.old);'])
                eval(['dN_' num2str(d) '= vertcat(dN_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.new);'])
                eval(['dO_corr_' num2str(d) '= vertcat(dO_corr_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.old(accuracies_d{i3,d}(old_d{i3,d}==1)==1,:));'])
                eval(['dN_corr_' num2str(d) '= vertcat(dN_corr_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.new(accuracies_d{i3,d}(old_d{i3,d}~=1)==1,:));'])
                end
                
                %         iRO   = vertcat(iRO,dat{i3,1}.plotdat_preproc.imm.RO);
                %         iRN   = vertcat(iRN,dat{i3,1}.plotdat_preproc.imm.RN);
                %         dRO   = vertcat(dRO,dat{i3,1}.plotdat_preproc.del.RO);
                %         dRN   = vertcat(dRN,dat{i3,1}.plotdat_preproc.del.RN);
                %
                %         inRO  = vertcat(inRO,dat{i3,1}.plotdat_preproc.imm.neuO);
                %         inRN  = vertcat(inRN,dat{i3,1}.plotdat_preproc.imm.neuN);
                %         dnRO  = vertcat(dnRO,dat{i3,1}.plotdat_preproc.del.neuO);
                %         dnRN  = vertcat(dnRN,dat{i3,1}.plotdat_preproc.del.neuN);
                
                
                %         % correct trls
                %         iR_c    = vertcat(iR_c,dat{i3,1}.plotdat_preproc.imm_c.R);
                %         inR_c   = vertcat(inR_c,dat{i3,1}.plotdat_preproc.imm_c.neu);
                %         dR_c    = vertcat(dR_c,dat{i3,1}.plotdat_preproc.del_c.R);
                %         dnR_c   = vertcat(dnR_c,dat{i3,1}.plotdat_preproc.del_c.neu);
                %
                %         iO_c    = vertcat(iO_c,dat{i3,1}.plotdat_preproc.imm_c.old);
                %         iN_c    = vertcat(iN_c,dat{i3,1}.plotdat_preproc.imm_c.new);
                %         dO_c    = vertcat(dO_c,dat{i3,1}.plotdat_preproc.del_c.old);
                %         dN_c    = vertcat(dN_c,dat{i3,1}.plotdat_preproc.del_c.new);
                %
                %         iRO_c   = vertcat(iRO_c,dat{i3,1}.plotdat_preproc.imm_c.RO);
                %         iRN_c   = vertcat(iRN_c,dat{i3,1}.plotdat_preproc.imm_c.RN);
                %         dRO_c   = vertcat(dRO_c,dat{i3,1}.plotdat_preproc.del_c.RO);
                %         dRN_c   = vertcat(dRN_c,dat{i3,1}.plotdat_preproc.del_c.RN);
                %
                %         inRO_c  = vertcat(inRO_c,dat{i3,1}.plotdat_preproc.imm_c.neuO);
                %         inRN_c  = vertcat(inRN_c,dat{i3,1}.plotdat_preproc.imm_c.neuN);
                %         dnRO_c  = vertcat(dnRO_c,dat{i3,1}.plotdat_preproc.del_c.neuO);
                %         dnRN_c  = vertcat(dnRN_c,dat{i3,1}.plotdat_preproc.del_c.neuN);
            catch
                if eval(['d' num2str(d) 'm(i3,1) == 1'])
                iR  = vertcat(iR,dat{i3,d}.plotdat_preproc.imm.rew(:,1:2701));
                inR = vertcat(inR,dat{i3,d}.plotdat_preproc.imm.neu(:,1:2701));
                iO    = vertcat(iO,dat{i3,d}.plotdat_preproc.imm.old(:,1:2701));
                iN    = vertcat(iN,dat{i3,d}.plotdat_preproc.imm.new(:,1:2701));
                iO_corr = vertcat(iO_corr,dat{i3,d}.plotdat_preproc.imm.old(accuracies_i{i3,d}(old_i{i3,d}==1)==1,1:2701));
                iN_corr = vertcat(iN_corr,dat{i3,d}.plotdat_preproc.imm.new(accuracies_i{i3,d}(old_i{i3,d}~=1)==1,1:2701));
                
                eval(['iR_' num2str(d) '= vertcat(iR_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.rew(:,1:2701));'])
                eval(['inR_' num2str(d) '= vertcat(inR_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.neu(:,1:2701));'])
                eval(['iO_' num2str(d) '= vertcat(iO_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.old(:,1:2701));'])
                eval(['iN_' num2str(d) '= vertcat(iN_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.new(:,1:2701));'])
                eval(['iO_corr_' num2str(d) '= vertcat(iO_corr_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.old(accuracies_i{i3,d}(old_i{i3,d}==1)==1,1:2701));'])
                eval(['iN_corr_' num2str(d) '= vertcat(iN_corr_' num2str(d) ',dat{i3,d}.plotdat_preproc.imm.new(accuracies_i{i3,d}(old_i{i3,d}~=1)==1,1:2701));'])                
                end
                if eval(['d' num2str(d) 'm(i3,2) == 2'])
                dR  = vertcat(dR,dat{i3,d}.plotdat_preproc.del.rew(:,1:2701));
                dnR = vertcat(dnR,dat{i3,d}.plotdat_preproc.del.neu(:,1:2701));
                dO    = vertcat(dO,dat{i3,d}.plotdat_preproc.del.old(:,1:2701));
                dN    = vertcat(dN,dat{i3,d}.plotdat_preproc.del.new(:,1:2701));
                dO_corr = vertcat(dO_corr,dat{i3,d}.plotdat_preproc.del.old(accuracies_d{i3,d}(old_d{i3,d}==1)==1,1:2701));
                dN_corr = vertcat(dN_corr,dat{i3,d}.plotdat_preproc.del.new(accuracies_d{i3,d}(old_d{i3,d}~=1)==1,1:2701));
                
                eval(['dR_' num2str(d) '= vertcat(dR_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.rew(:,1:2701));'])
                eval(['dnR_' num2str(d) '= vertcat(dnR_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.neu(:,1:2701));'])
                eval(['dO_' num2str(d) '= vertcat(dO_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.old(:,1:2701));'])
                eval(['dN_' num2str(d) '= vertcat(dN_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.new(:,1:2701));'])
                eval(['dO_corr_' num2str(d) '= vertcat(dO_corr_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.old(accuracies_d{i3,d}(old_d{i3,d}==1)==1,1:2701));'])
                eval(['dN_corr_' num2str(d) '= vertcat(dN_corr_' num2str(d) ',dat{i3,d}.plotdat_preproc.del.new(accuracies_d{i3,d}(old_d{i3,d}~=1)==1,1:2701));'])
                end
            end
            end
        else
        end
    end
    
end



cd(path_fig) % saving path for the figures

% fancier: shadedErrorBar

% plot with expectancy: Rewards, stimulus
% stim
maxy = 1; miny = -0.2;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
brightred = [1 , 0 , 0]; % red is first condition
brightblue = [0 , 0 , 1];
darkred = [0.5 , 0 , 0]; 
darkblue = [0 , 0 , 0.5];

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iR_2,{@nanmean,errorbar},'lineprops',{'-','LineWidth',7,'color',brightred,'markerfacecolor',brightred}); hold on; % rare reward vs frequent reward
shadedErrorBar((1:2701),iR_1,{@nanmean,errorbar},'lineprops',{'--','LineWidth',7,'color',brightred,'markerfacecolor',brightred}); hold on;
shadedErrorBar((1:2701),inR_1,{@nanmean,errorbar},'lineprops',{'-','LineWidth',7,'color',brightblue,'markerfacecolor',brightblue}); hold on; % rare reward vs frequent reward
shadedErrorBar((1:2701),inR_2,{@nanmean,errorbar},'lineprops',{'--','LineWidth',7,'color',brightblue,'markerfacecolor',brightblue}); hold on;
% H = ttest2(iR_2,iR_1,'Alpha',0.05); sigbar = H*-0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'m-','LineWidth',7);clear H sigbar
% H = ttest2(inR_1,inR_2,'Alpha',0.05); sigbar = H*0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'c-','LineWidth',7);
% H = ttest2(iR_1,inR_2,'Alpha',0.05); sigbar = H*0.2; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-.','LineWidth',7);
% H = ttest2(iR_2,inR_1,'Alpha',0.05); sigbar = H*0.3; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-','LineWidth',7);
H = ttest(iR_2,iR_1(1:size(iR_2,1),1:size(iR_2,2)),'Alpha',0.05); sigbar = H*-0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'m-','LineWidth',7);clear H sigbar
H = ttest(inR_1,inR_2(1:size(inR_1,1),1:size(inR_1,2)),'Alpha',0.05); sigbar = H*0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'c-','LineWidth',7);
H = ttest(iR_1(1:size(inR_2,1),1:size(inR_2,2)),inR_2,'Alpha',0.05); sigbar = H*0.2; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-.','LineWidth',7);
H = ttest(iR_2,inR_1(1:size(iR_2,1),1:size(iR_2,2)),'Alpha',0.05); sigbar = H*0.3; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-','LineWidth',7);

xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)
shadedErrorBar((1:2701),dR_2,{@nanmean,errorbar},'lineprops',{'-','LineWidth',7,'color',brightred,'markerfacecolor',brightred}); hold on; % rare reward vs frequent reward
shadedErrorBar((1:2701),dR_1,{@nanmean,errorbar},'lineprops',{'--','LineWidth',7,'color',brightred,'markerfacecolor',brightred}); hold on;
shadedErrorBar((1:2701),dnR_1,{@nanmean,errorbar},'lineprops',{'-','LineWidth',7,'color',brightblue,'markerfacecolor',brightblue}); hold on; % rare reward vs frequent reward
shadedErrorBar((1:2701),dnR_2,{@nanmean,errorbar},'lineprops',{'--','LineWidth',7,'color',brightblue,'markerfacecolor',brightblue}); hold on;
% H = ttest2(dR_2,dR_1,'Alpha',0.05); sigbar = H*-0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'m-','LineWidth',7);clear H sigbar
% H = ttest2(dnR_1,dnR_2,'Alpha',0.05); sigbar = H*0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'c-','LineWidth',7);
% H = ttest2(dR_1,dnR_2,'Alpha',0.05); sigbar = H*0.2; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-.','LineWidth',7);
% H = ttest2(dR_2,dnR_1,'Alpha',0.05); sigbar = H*0.3; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-','LineWidth',7);
H = ttest(dR_2,dR_1(1:size(dR_2,1),1:size(dR_2,2)),'Alpha',0.05); sigbar = H*-0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'m-','LineWidth',7);clear H sigbar
H = ttest(dnR_1,dnR_2(1:size(dnR_1,1),1:size(dnR_1,2)),'Alpha',0.05); sigbar = H*0.1; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'c-','LineWidth',7);
H = ttest(dR_1,dnR_2(1:size(dR_1,1),1:size(dR_1,2)),'Alpha',0.05); sigbar = H*0.2; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-.','LineWidth',7);
H = ttest(dR_2(1:size(dnR_1,1),1:size(dnR_1,2)),dnR_1,'Alpha',0.05); sigbar = H*0.3; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-','LineWidth',7);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('delayed','Fontsize',20,'Fontweight','bold')
clear H sigbar


% rew vs. non-rew
maxy = 1; miny = -0.3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iR,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),inR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iR,inR,'Alpha',0.05); sigbar = H*-0.01; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)
shadedErrorBar((1:2701),dR,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),dnR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(dR,dnR,'Alpha',0.05); sigbar = H*-0.01; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('delayed','Fontsize',20,'Fontweight','bold')
suptitle('Reward(red) vs. Non-Reward(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
saveas(gcf,'GroupLvl_memtest_pilot3T_Reward(red)_vs_Non-Reward(blue).fig')
print(['GroupLvl_memtest_pilot3T_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
clear fig H sigbar



% old vs. new
maxy = 1; miny = -0.3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),iN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iO,iN,'Alpha',0.05); sigbar = H*-0.01; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),dO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),dN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(dO,dN,'Alpha',0.05); sigbar = H*-0.01; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('delayed','Fontsize',20,'Fontweight','bold')
suptitle('OLD(red) vs. NEW(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
saveas(gcf,'GroupLvl_memtest_pilot3T_New_vs_Old_pilot.fig')
print(['GroupLvl_memtest_pilot3T_New_vs_Old_pilot'],'-dpdf','-r0')
clear fig H sigbar



% correct trials: old vs. new
maxy = 1; miny = -0.3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (nanstd(data)./sqrt(size(data,1)));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iO_corr,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),iN_corr,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iO_corr,iN_corr,'Alpha',0.05); sigbar = H*-0.01; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),dO_corr,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),dN_corr,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(dO_corr,dN_corr,'Alpha',0.05); sigbar = H*-0.01; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('delayed','Fontsize',20,'Fontweight','bold')
suptitle('Correct Trials: OLD(red) vs. NEW(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
saveas(gcf,'GroupLvl_memtest_pilot3T_New_vs_Old_corrects.fig')
print(['GroupLvl_memtest_pilot3T_New_vs_Old_corrects'],'-dpdf','-r0')
clear fig H sigbar



% 
% % only no-reward trls: old v. new
% 
% maxy = 3; miny = -3;
% onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
% 
% figure;
% % subplot(1,2,1);
% plot(onsetbar,y,'g-','linewidth',2)
% 
% shadedErrorBar((1:2701),inRO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),inRN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(inRO,inRN,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('immediate','Fontsize',20,'Fontweight','bold')
% clear H sigbar
% % % feedback seg: reward vs. non-reward
% % subplot(1,2,2);
% % plot(onsetbar,y,'g-','linewidth',2)
% %
% % shadedErrorBar((1:2701),dnRO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% % shadedErrorBar((1:2701),dnRN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% % H = ttest2(dnRO,dnRN,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% % xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% % xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% % ylabel('zscores','Fontsize',20,'Fontweight','bold');
% % title('delayed','Fontsize',20,'Fontweight','bold')
% % suptitle('no-reward: OLD(red) vs. NEW(blue)')
% 
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_memtest_pilot'  '_noReward_New_vs_Old_pilot'],'-dpdf','-r0')
% clear fig H sigbar




% %% immediate & delayed
% 
% %% preparation
% 
% clear;clc;close all;
% % warning('off')
% addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20180613'))
% 
% path_gen    = '/Users/yeojin/Desktop/';
% 
% path_par    = [ path_gen 'E_data/' ];
% path_dat    = [ path_par 'EA_raw/EAA_pupil/pilot_3T/'];
% path_behav  = [ path_par 'EA_raw/EAC_behav/pilot_3T/' ];
% path_clean  = [ path_par 'EB_cleaned/EBB_pupil/pilot_3T/' ];
% path_fig    = [ path_gen 'C_writings/CB_figures/pilot_3T/zscores/' ];
% path_stim   = [ path_par 'EA_raw/EAA_pupil/stiminfo_all_20181004.mat' ];
% 
% % create file list and VP ID list
% VPs = [1109 1210 1111 1212 1115];
% session = 1;
% 
% %% now run
% 
% for i2 = 1:length(VPs)
%     dat{i2,1} = load([ path_clean 'rr' num2str(session) num2str(VPs(i2)) '_preproc_dat.mat']);
% end
% 
% iR=[];inR=[];dR=[];dnR=[];
% for i3 = 1:size(dat,1) % number of IDs
%     if isempty(iR)
%         iR      = dat{i3,1}.plotdat_preproc.imm.R;
%         inR     = dat{i3,1}.plotdat_preproc.imm.nR;
%         dR      = dat{i3,1}.plotdat_preproc.del.R(:,1:2701);
%         dnR     = dat{i3,1}.plotdat_preproc.del.nR(:,1:2701);
%         
%         iO      = dat{i3,1}.plotdat_preproc.imm.old;
%         iN      = dat{i3,1}.plotdat_preproc.imm.new;
%         dO      = dat{i3,1}.plotdat_preproc.del.old(:,1:2701);
%         dN      = dat{i3,1}.plotdat_preproc.del.new(:,1:2701);
%         
%         iRO     = dat{i3,1}.plotdat_preproc.imm.RO;
%         iRN     = dat{i3,1}.plotdat_preproc.imm.RN;
%         dRO     = dat{i3,1}.plotdat_preproc.del.RO(:,1:2701);
%         dRN     = dat{i3,1}.plotdat_preproc.del.RN(:,1:2701);
%         
%         inRO    = dat{i3,1}.plotdat_preproc.imm.nRO;
%         inRN    = dat{i3,1}.plotdat_preproc.imm.nRN;
%         dnRO    = dat{i3,1}.plotdat_preproc.del.nRO(:,1:2701);
%         dnRN    = dat{i3,1}.plotdat_preproc.del.nRN(:,1:2701);
%         
%         %         correct trls
%         iR_c    = dat{i3,1}.plotdat_preproc.imm_c.R;
%         inR_c   = dat{i3,1}.plotdat_preproc.imm_c.nR;
%         dR_c    = dat{i3,1}.plotdat_preproc.del_c.R(:,1:2701);
%         dnR_c   = dat{i3,1}.plotdat_preproc.del_c.nR(:,1:2701);
%         
%         iO_c    = dat{i3,1}.plotdat_preproc.imm_c.old;
%         iN_c    = dat{i3,1}.plotdat_preproc.imm_c.new;
%         dO_c    = dat{i3,1}.plotdat_preproc.del_c.old(:,1:2701);
%         dN_c    = dat{i3,1}.plotdat_preproc.del_c.new(:,1:2701);
%         
%         iRO_c   = dat{i3,1}.plotdat_preproc.imm_c.RO;
%         iRN_c   = dat{i3,1}.plotdat_preproc.imm_c.RN;
%         dRO_c   = dat{i3,1}.plotdat_preproc.del_c.RO(:,1:2701);
%         dRN_c   = dat{i3,1}.plotdat_preproc.del_c.RN(:,1:2701);
%         
%         inRO_c  = dat{i3,1}.plotdat_preproc.imm_c.nRO;
%         inRN_c  = dat{i3,1}.plotdat_preproc.imm_c.nRN;
%         dnRO_c  = dat{i3,1}.plotdat_preproc.del_c.nRO(:,1:2701);
%         dnRN_c  = dat{i3,1}.plotdat_preproc.del_c.nRN(:,1:2701);
%         
%         
%     else
%         iR  = vertcat(iR,dat{i3,1}.plotdat_preproc.imm.R);
%         inR = vertcat(inR,dat{i3,1}.plotdat_preproc.imm.nR);
%         dR  = vertcat(dR,dat{i3,1}.plotdat_preproc.del.R(:,1:2701));
%         dnR = vertcat(dnR,dat{i3,1}.plotdat_preproc.del.nR(:,1:2701));
%         
%         iO    = vertcat(iO,dat{i3,1}.plotdat_preproc.imm.old);
%         iN    = vertcat(iN,dat{i3,1}.plotdat_preproc.imm.new);
%         dO    = vertcat(dO,dat{i3,1}.plotdat_preproc.del.old(:,1:2701));
%         dN    = vertcat(dN,dat{i3,1}.plotdat_preproc.del.new(:,1:2701));
%         
%         iRO   = vertcat(iRO,dat{i3,1}.plotdat_preproc.imm.RO);
%         iRN   = vertcat(iRN,dat{i3,1}.plotdat_preproc.imm.RN);
%         dRO   = vertcat(dRO,dat{i3,1}.plotdat_preproc.del.RO(:,1:2701));
%         dRN   = vertcat(dRN,dat{i3,1}.plotdat_preproc.del.RN(:,1:2701));
%         
%         inRO  = vertcat(inRO,dat{i3,1}.plotdat_preproc.imm.nRO);
%         inRN  = vertcat(inRN,dat{i3,1}.plotdat_preproc.imm.nRN);
%         dnRO  = vertcat(dnRO,dat{i3,1}.plotdat_preproc.del.nRO(:,1:2701));
%         dnRN  = vertcat(dnRN,dat{i3,1}.plotdat_preproc.del.nRN(:,1:2701));
%         
%         
%         %         correct trls
%         iR_c    = vertcat(iR_c,dat{i3,1}.plotdat_preproc.imm_c.R);
%         inR_c   = vertcat(inR_c,dat{i3,1}.plotdat_preproc.imm_c.nR);
%         dR_c    = vertcat(dR_c,dat{i3,1}.plotdat_preproc.del_c.R(:,1:2701));
%         dnR_c   = vertcat(dnR_c,dat{i3,1}.plotdat_preproc.del_c.nR(:,1:2701));
%         
%         iO_c    = vertcat(iO_c,dat{i3,1}.plotdat_preproc.imm_c.old);
%         iN_c    = vertcat(iN_c,dat{i3,1}.plotdat_preproc.imm_c.new);
%         dO_c    = vertcat(dO_c,dat{i3,1}.plotdat_preproc.del_c.old(:,1:2701));
%         dN_c    = vertcat(dN_c,dat{i3,1}.plotdat_preproc.del_c.new(:,1:2701));
%         
%         iRO_c   = vertcat(iRO_c,dat{i3,1}.plotdat_preproc.imm_c.RO);
%         iRN_c   = vertcat(iRN_c,dat{i3,1}.plotdat_preproc.imm_c.RN);
%         dRO_c   = vertcat(dRO_c,dat{i3,1}.plotdat_preproc.del_c.RO(:,1:2701));
%         dRN_c   = vertcat(dRN_c,dat{i3,1}.plotdat_preproc.del_c.RN(:,1:2701));
%         
%         inRO_c  = vertcat(inRO_c,dat{i3,1}.plotdat_preproc.imm_c.nRO);
%         inRN_c  = vertcat(inRN_c,dat{i3,1}.plotdat_preproc.imm_c.nRN);
%         dnRO_c  = vertcat(dnRO_c,dat{i3,1}.plotdat_preproc.del_c.nRO(:,1:2701));
%         dnRN_c  = vertcat(dnRN_c,dat{i3,1}.plotdat_preproc.del_c.nRN(:,1:2701));
%         
%     end
%     
% end
% 
% 
% cd(path_fig) % saving path for the figures
% 
% % fancier: shadedErrorBar
% % rew vs. non-rew
% maxy = 3; miny = -3;
% onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
% 
% figure;
% subplot(1,2,1);
% plot(onsetbar,y,'g-','linewidth',2)
% 
% shadedErrorBar((1:2701),iR,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),inR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(iR,inR,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('immediate','Fontsize',20,'Fontweight','bold')
% clear H sigbar
% 
% subplot(1,2,2);
% plot(onsetbar,y,'g-','linewidth',2)
% shadedErrorBar((1:2701),dR,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),dnR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(dR,dnR,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('delayed','Fontsize',20,'Fontweight','bold')
% suptitle('Reward(red) vs. Non-Reward(blue)')
% 
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_memtest_pilot'  '_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
% clear fig H sigbar
% 
% 
% % fancier: shadedErrorBar
% 
% % only correct trials: rew vs. non-rew
% maxy = 3; miny = -3;
% onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
% 
% figure;
% subplot(1,2,1);
% plot(onsetbar,y,'g-','linewidth',2)
% 
% shadedErrorBar((1:2701),iR_c,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),inR_c,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(iR_c,inR_c,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('immediate','Fontsize',20,'Fontweight','bold')
% clear H sigbar
% 
% subplot(1,2,2);
% plot(onsetbar,y,'g-','linewidth',2)
% 
% shadedErrorBar((1:2701),dR_c,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),dnR_c,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(dR_c,dnR_c,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('delayed','Fontsize',20,'Fontweight','bold')
% suptitle('Corrects: Reward(red) vs. Non-Reward(blue)')
% 
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_memtest_pilot'  '_CorrectTrials_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
% clear fig H sigbar
% 
% 
% 
% % old vs. new
% maxy = 3; miny = -3;
% onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
% 
% figure;
% subplot(1,2,1);
% plot(onsetbar,y,'g-','linewidth',2)
% 
% shadedErrorBar((1:2701),iO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),iN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(iO,iN,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('immediate','Fontsize',20,'Fontweight','bold')
% clear H sigbar
% 
% subplot(1,2,2);
% plot(onsetbar,y,'g-','linewidth',2)
% 
% shadedErrorBar((1:2701),dO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),dN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(dO,dN,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('delayed','Fontsize',20,'Fontweight','bold')
% suptitle('OLD(red) vs. NEW(blue)')
% 
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_memtest_pilot'  '_New_vs_Old_pilot'],'-dpdf','-r0')
% clear fig H sigbar
% 
% 
% 
% 
% % only no-reward trls: old v. new
% 
% maxy = 3; miny = -3;
% onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
% errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
% 
% figure;
% subplot(1,2,1);
% plot(onsetbar,y,'g-','linewidth',2)
% 
% shadedErrorBar((1:2701),inRO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),inRN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(inRO,inRN,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('immediate','Fontsize',20,'Fontweight','bold')
% clear H sigbar
% % feedback seg: reward vs. non-reward
% subplot(1,2,2);
% plot(onsetbar,y,'g-','linewidth',2)
% 
% shadedErrorBar((1:2701),dnRO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),dnRN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(dnRO,dnRN,'Alpha',0.05); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:1000)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('delayed','Fontsize',20,'Fontweight','bold')
% suptitle('no-reward: OLD(red) vs. NEW(blue)')
% 
% % save
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(['GroupLvl_memtest_pilot'  '_noReward_New_vs_Old_pilot'],'-dpdf','-r0')
% clear fig H sigbar
