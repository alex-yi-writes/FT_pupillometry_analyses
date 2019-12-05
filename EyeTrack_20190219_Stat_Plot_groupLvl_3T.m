%% Plot for all participants altogether

%% work log

%   27_11_2018  created the script

%% Main task

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
VPs = [2201 2202 2203 2204 2205 2206 2207 2208];
tag = [0 0; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 0];

%% now run

for i2 = 1:length(VPs)
    for d = 1:2
        if tag(i2,d) == 0
        else
            dat{i2,d} = load([ path_clean 're' num2str(d) num2str(VPs(i2)) '_preproc_dat.mat']);
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
                
                %         rew_h  = days{i3,1}.plotdat_preproc.stim.hiMem_rew;
                %         pun_h = days{i3,1}.plotdat_preproc.stim.hiMem_pun;
                % %         fR_h  = days{i3,1}.plotdat_preproc.fb.R_mH;
                % %         fnR_h = days{i3,1}.plotdat_preproc.fb.nR_mH;
                %
                %         rew_l  = days{i3,1}.plotdat_preproc.stim.loMem_rew;
                %         pun_l = days{i3,1}.plotdat_preproc.stim.loMem_pun;
                % %         fR_l  = days{i3,1}.plotdat_preproc.fb.R_mL;
                % %         fnR_l = days{i3,1}.plotdat_preproc.fb.nR_mL;
                %
                %         highmem = days{i3,1}.plotdat_preproc.stim.hiMem;
                %         lowmem  = days{i3,1}.plotdat_preproc.stim.loMem;
                %
                %         smH{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mH,days{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mL,days{i3,1}.plotdat_preproc.stim.R_mL);
                
            else
                sR  = vertcat(sR,dat{i3,d}.plotdat_preproc.stim.rew);
                sN = vertcat(sN,dat{i3,d}.plotdat_preproc.stim.neu);
                fR  = vertcat(fR,dat{i3,d}.plotdat_preproc.fb.rew);
                fN = vertcat(fN,dat{i3,d}.plotdat_preproc.fb.neu);
                
                %         rew_h  = vertcat(rew_h,days{i3,1}.plotdat_preproc.stim.hiMem_rew);
                %         pun_h = vertcat(pun_h,days{i3,1}.plotdat_preproc.stim.hiMem_pun);
                % %         fR_h  = vertcat(fR_h,days{i3,1}.plotdat_preproc.fb.R_mH);
                % %         fnR_h = vertcat(fnR_h,days{i3,1}.plotdat_preproc.fb.nR_mH);
                %
                %         rew_l  = vertcat(rew_l,days{i3,1}.plotdat_preproc.stim.loMem_rew);
                %         pun_l = vertcat(pun_l,days{i3,1}.plotdat_preproc.stim.loMem_pun);
                % %         fR_l  = vertcat(fR_l,days{i3,1}.plotdat_preproc.fb.R_mL);
                % %         fnR_l = vertcat(fnR_l,days{i3,1}.plotdat_preproc.fb.nR_mL);
                %
                %         highmem = vertcat(highmem,days{i3,1}.plotdat_preproc.stim.hiMem);
                %         lowmem  = vertcat(lowmem,days{i3,1}.plotdat_preproc.stim.loMem);
                %
                %         smH{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mH,days{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mL,days{i3,1}.plotdat_preproc.stim.R_mL);
                
                
            end
        end
    end
    
end

sR_all = sR; sN_all = sN; fR_all = fR; fN_all = fN;


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
                
                %         rew_h  = days{i3,1}.plotdat_preproc.stim.hiMem_rew;
                %         pun_h = days{i3,1}.plotdat_preproc.stim.hiMem_pun;
                % %         fR_h  = days{i3,1}.plotdat_preproc.fb.R_mH;
                % %         fnR_h = days{i3,1}.plotdat_preproc.fb.nR_mH;
                %
                %         rew_l  = days{i3,1}.plotdat_preproc.stim.loMem_rew;
                %         pun_l = days{i3,1}.plotdat_preproc.stim.loMem_pun;
                % %         fR_l  = days{i3,1}.plotdat_preproc.fb.R_mL;
                % %         fnR_l = days{i3,1}.plotdat_preproc.fb.nR_mL;
                %
                %         highmem = days{i3,1}.plotdat_preproc.stim.hiMem;
                %         lowmem  = days{i3,1}.plotdat_preproc.stim.loMem;
                %
                %         smH{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mH,days{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mL,days{i3,1}.plotdat_preproc.stim.R_mL);
                
            else
                sR  = vertcat(sR,dat{i3,d}.plotdat_preproc.stim.rew);
                sN = vertcat(sN,dat{i3,d}.plotdat_preproc.stim.neu);
                fR  = vertcat(fR,dat{i3,d}.plotdat_preproc.fb.rew);
                fN = vertcat(fN,dat{i3,d}.plotdat_preproc.fb.neu);
                
                %         rew_h  = vertcat(rew_h,days{i3,1}.plotdat_preproc.stim.hiMem_rew);
                %         pun_h = vertcat(pun_h,days{i3,1}.plotdat_preproc.stim.hiMem_pun);
                % %         fR_h  = vertcat(fR_h,days{i3,1}.plotdat_preproc.fb.R_mH);
                % %         fnR_h = vertcat(fnR_h,days{i3,1}.plotdat_preproc.fb.nR_mH);
                %
                %         rew_l  = vertcat(rew_l,days{i3,1}.plotdat_preproc.stim.loMem_rew);
                %         pun_l = vertcat(pun_l,days{i3,1}.plotdat_preproc.stim.loMem_pun);
                % %         fR_l  = vertcat(fR_l,days{i3,1}.plotdat_preproc.fb.R_mL);
                % %         fnR_l = vertcat(fnR_l,days{i3,1}.plotdat_preproc.fb.nR_mL);
                %
                %         highmem = vertcat(highmem,days{i3,1}.plotdat_preproc.stim.hiMem);
                %         lowmem  = vertcat(lowmem,days{i3,1}.plotdat_preproc.stim.loMem);
                %
                %         smH{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mH,days{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mL,days{i3,1}.plotdat_preproc.stim.R_mL);
                
                
            end
        end
    end
    
end

sR_1 = sR; sN_1 = sN; fR_1 = fR; fN_1 = fN;

% session 2 collapsed
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
                
                %         rew_h  = days{i3,1}.plotdat_preproc.stim.hiMem_rew;
                %         pun_h = days{i3,1}.plotdat_preproc.stim.hiMem_pun;
                % %         fR_h  = days{i3,1}.plotdat_preproc.fb.R_mH;
                % %         fnR_h = days{i3,1}.plotdat_preproc.fb.nR_mH;
                %
                %         rew_l  = days{i3,1}.plotdat_preproc.stim.loMem_rew;
                %         pun_l = days{i3,1}.plotdat_preproc.stim.loMem_pun;
                % %         fR_l  = days{i3,1}.plotdat_preproc.fb.R_mL;
                % %         fnR_l = days{i3,1}.plotdat_preproc.fb.nR_mL;
                %
                %         highmem = days{i3,1}.plotdat_preproc.stim.hiMem;
                %         lowmem  = days{i3,1}.plotdat_preproc.stim.loMem;
                %
                %         smH{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mH,days{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mL,days{i3,1}.plotdat_preproc.stim.R_mL);
                
            else
                sR  = vertcat(sR,dat{i3,d}.plotdat_preproc.stim.rew);
                sN = vertcat(sN,dat{i3,d}.plotdat_preproc.stim.neu);
                fR  = vertcat(fR,dat{i3,d}.plotdat_preproc.fb.rew);
                fN = vertcat(fN,dat{i3,d}.plotdat_preproc.fb.neu);
                
                %         rew_h  = vertcat(rew_h,days{i3,1}.plotdat_preproc.stim.hiMem_rew);
                %         pun_h = vertcat(pun_h,days{i3,1}.plotdat_preproc.stim.hiMem_pun);
                % %         fR_h  = vertcat(fR_h,days{i3,1}.plotdat_preproc.fb.R_mH);
                % %         fnR_h = vertcat(fnR_h,days{i3,1}.plotdat_preproc.fb.nR_mH);
                %
                %         rew_l  = vertcat(rew_l,days{i3,1}.plotdat_preproc.stim.loMem_rew);
                %         pun_l = vertcat(pun_l,days{i3,1}.plotdat_preproc.stim.loMem_pun);
                % %         fR_l  = vertcat(fR_l,days{i3,1}.plotdat_preproc.fb.R_mL);
                % %         fnR_l = vertcat(fnR_l,days{i3,1}.plotdat_preproc.fb.nR_mL);
                %
                %         highmem = vertcat(highmem,days{i3,1}.plotdat_preproc.stim.hiMem);
                %         lowmem  = vertcat(lowmem,days{i3,1}.plotdat_preproc.stim.loMem);
                %
                %         smH{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mH,days{i3,1}.plotdat_preproc.stim.R_mH);
                %         smL{i3,1} = vertcat(days{i3,1}.plotdat_preproc.stim.nR_mL,days{i3,1}.plotdat_preproc.stim.R_mL);
                
                
            end
        end
    end
    
end

sR_2 = sR; sN_2 = sN; fR_2 = fR; fN_2 = fN;


cd(path_fig) % saving path for the figures

% fancier: shadedErrorBar

% all sessions
% stim
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_all,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_all,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_all,sN_all); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
H = ttest2(fR_all,fN_all); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('feedback','Fontsize',20,'Fontweight','bold')

suptitle(['Reward(red) vs. neutral(blue): Session 1 and 2'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_all_sessions'],'-dpdf','-r0')
clear fig H sigbar




% each sessions

for d = 1:2

% stim
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

eval(['shadedErrorBar((1:2701),sR_' num2str(d) ',{@nanmean,errorbar},''lineprops'',{''r-o'',''markerfacecolor'',''r''})']); hold on;
eval(['shadedErrorBar((1:2701),sN_' num2str(d) ',{@nanmean,errorbar},''lineprops'',{''b-o'',''markerfacecolor'',''b''})']); hold on;
eval(['H = ttest2(sR_' num2str(d) ',sN_' num2str(d) ');'])
sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
 sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('feedback','Fontsize',20,'Fontweight','bold')

suptitle(['Reward(red) vs. neutral(blue): Session ' num2str(d)])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_session' num2str(d)],'-dpdf','-r0')
clear fig H sigbar


end



% session 1 & 2
% stim_1
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
subplot(2,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_1,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_1,sN_1); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
H = ttest2(fR_1,fN_1); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('session1 feedback','Fontsize',20,'Fontweight','bold')

% stim_2
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

subplot(2,2,3);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_2,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_2,sN_2); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
H = ttest2(fR_2,fN_2); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('session2 feedback','Fontsize',20,'Fontweight','bold')


suptitle(['Reward(red) vs. neutral(blue): Session 1 and 2'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_3T_2sess_maintask_Reward(red)_vs_Neutral(blue)_all_sessions-cross'],'-dpdf','-r0')
clear fig H sigbar



% session 1 vs 2
% stim_reward
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
subplot(2,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sR_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sR_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sR_1,sR_2); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
H = ttest2(fR_1,fR_2); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('Feedback & Reward','Fontsize',17,'Fontweight','bold')

% stim_neutral
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

subplot(2,2,3);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),sN_1,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),sN_2,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(sN_1,sN_2); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
H = ttest2(fN_1,fN_2); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('Feedback & Neutral','Fontsize',17,'Fontweight','bold')

suptitle(['Session 1(red) vs. Session 2(blue)'])

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
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
% H = ttest2(rew_h,rew_l); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
% % H = ttest2(fR_h,fR_l); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
% H = ttest2(pun_h,pun_l); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
% % H = ttest2(fnR_h,fnR_l); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
% H = ttest2(highmem,lowmem); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
%     errorbar = @ (data) (std(data)/sqrt(length(VPs)));
%     
%     subplot(ceil(sqrt(length(VPs))),ceil(sqrt(length(VPs))),i4);
%     plot(onsetbar,y,'g-','linewidth',2)
%     
%     shadedErrorBar((1:2701),dat{i4,1}.plotdat_preproc.stim.rew,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:2701),dat{i4,1}.plotdat_preproc.stim.pun,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(dat{i4,1}.plotdat_preproc.stim.rew,dat{i4,1}.plotdat_preproc.stim.pun); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
%     errorbar = @ (data) (std(data)/sqrt(length(VPs)));
%
%     subplot(ceil(sqrt(length(VPs))),ceil(sqrt(length(VPs))),i4);
%     plot(onsetbar,y,'g-','linewidth',2)
%
%     shadedErrorBar((1:2701),smH{i4,1},{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:2701),smL{i4,1},{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(smH{i4,1},smL{i4,1}); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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
%     errorbar = @ (data) (std(data)/sqrt(length(VPs)));
%
%     subplot(ceil(sqrt(length(VPs))),ceil(sqrt(length(VPs))),i4);
%     plot(onsetbar,y,'g-','linewidth',2)
%
%     shadedErrorBar((1:2701),days{i4,1}.plotdat_preproc.fb.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:2701),days{i4,1}.plotdat_preproc.fb.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(days{i4,1}.plotdat_preproc.fb.R,days{i4,1}.plotdat_preproc.fb.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
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


%% Oddball

%% preparation

clear;clc;close all;
% warning('off')
addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20180613'))

path_gen    = '/Users/yeojin/Desktop/';

path_par    = [ path_gen 'E_data/' ];
path_dat    = [ path_par 'EA_raw/EAA_pupil/pilot_3T/'];
path_behav  = [ path_par 'EA_raw/EAC_behav/pilot_3T/' ];
path_clean  = [ path_par 'EB_cleaned/EBB_pupil/pilot_3T/' ];
path_fig    = [ path_gen 'C_writings/CB_figures/pilot_3T/zscores/' ];
path_stim   = [ path_par 'EA_raw/EAA_pupil/stiminfo_all_20181004.mat' ];

% create file list and VP ID list
VPs = [1204 1107 1109 1210 1111 1212 1113 1214 1216];

%% now run

for i2 = 1:length(VPs)
    dat{i2,1} = load([ path_clean '/ob' num2str(VPs(i2)) '_preproc_dat.mat']);
end

odd=[];stds=[];%fR=[];fnR=[];
for i3 = 1:size(dat,1) % number of IDs
    if isempty(odd)
        odd  = dat{i3,1}.plotdat_preproc.stim.odd;
        stds = dat{i3,1}.plotdat_preproc.stim.std;
        
    else
        odd  = vertcat(odd,dat{i3,1}.plotdat_preproc.stim.odd);
        stds = vertcat(stds,dat{i3,1}.plotdat_preproc.stim.std);
        
    end
    
end


cd(path_fig) % saving path for the figures

% fancier: shadedErrorBar
% stimulus seg: oddball vs. standard

maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),odd,{ @nanmean , errorbar },'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),stds,{ @nanmean , errorbar },'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(odd,stds); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('stimulus(Oddball-red, standard-blue)','Fontsize',20,'Fontweight','bold')
clear H sigbar

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_3T_OB_Oddball(red)_vs_Standard(blue)'],'-dpdf','-r0')
clear fig H sigbar


% single participants in one plot!
% stim, reward vs. stdishment
figure; hold on;
for i4 = 1:length(VPs)
    maxy = 3; miny = -3;
    onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
    errorbar = @ (data) (std(data)/sqrt(length(VPs)));
    
    subplot(ceil(sqrt(length(VPs))),ceil(sqrt(length(VPs))),i4);
    plot(onsetbar,y,'g-','linewidth',2)
    
    shadedErrorBar((1:2701),dat{i4,1}.plotdat_preproc.stim.odd,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:2701),dat{i4,1}.plotdat_preproc.stim.std,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(dat{i4,1}.plotdat_preproc.stim.odd,dat{i4,1}.plotdat_preproc.stim.std); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',8,'Fontweight','bold');
    ylabel('zscores','Fontsize',8,'Fontweight','bold');
    title(num2str(VPs(i4)))
    clear H sigbar
    
end
suptitle('stimulus seg: Oddball(red) vs. Standard(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_3T_maintask_Indiv_stim_Oddball(red)_vs_Standard(blue)'],'-dpdf','-r0')
clear fig H sigbar


%% Memory test

%% immediate

%% preparation

clear;clc;close all;
% warning('off')
addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20180613'))

path_gen    = '/Users/yeojin/Desktop/';

path_par    = [ path_gen 'E_data/' ];
path_dat    = [ path_par 'EA_raw/EAA_pupil/pilot_3T/'];
path_behav  = [ path_par 'EA_raw/EAC_behav/pilot_3T/' ];
path_clean  = [ path_par 'EB_cleaned/EBB_pupil/pilot_3T/' ];
path_fig    = [ path_gen 'C_writings/CB_figures/pilot_3T/zscores/' ];
path_stim   = [ path_par 'EA_raw/EAA_pupil/stiminfo_all_20181004.mat' ];

% create file list and VP ID list
VPs = [1204 1107 1208 1109 1210 1111 1212 1113 1214 1216];
session = 1;

%% now run

for i2 = 1:length(VPs)
    dat{i2,1} = load([ path_clean 'rr' num2str(session) num2str(VPs(i2)) '_preproc_dat.mat']);
end

iR=[];inR=[];%dR=[];dnR=[];
for i3 = 1:size(dat,1) % number of IDs
    if isempty(iR)
        iR      = dat{i3,1}.plotdat_preproc.im.R;
        inR     = dat{i3,1}.plotdat_preproc.im.nR;
        %         dR      = days{i3,1}.plotdat_preproc.de.R;
        %         dnR     = days{i3,1}.plotdat_preproc.de.nR;
        
        iO      = dat{i3,1}.plotdat_preproc.im.old;
        iN      = dat{i3,1}.plotdat_preproc.im.new;
        %         dO      = days{i3,1}.plotdat_preproc.de.old;
        %         dN      = days{i3,1}.plotdat_preproc.de.new;
        
        iRO     = dat{i3,1}.plotdat_preproc.im.RO;
        iRN     = dat{i3,1}.plotdat_preproc.im.RN;
        %         dRO     = days{i3,1}.plotdat_preproc.de.RO;
        %         dRN     = days{i3,1}.plotdat_preproc.de.RN;
        
        inRO    = dat{i3,1}.plotdat_preproc.im.nRO;
        inRN    = dat{i3,1}.plotdat_preproc.im.nRN;
        %         dnRO    = days{i3,1}.plotdat_preproc.de.nRO;
        %         dnRN    = days{i3,1}.plotdat_preproc.de.nRN;
        
        % correct trls
        iR_c    = dat{i3,1}.plotdat_preproc.im_c.R;
        inR_c   = dat{i3,1}.plotdat_preproc.im_c.nR;
        %         dR_c    = days{i3,1}.plotdat_preproc.de_c.R;
        %         dnR_c   = days{i3,1}.plotdat_preproc.de_c.nR;
        
        iO_c    = dat{i3,1}.plotdat_preproc.im_c.old;
        iN_c    = dat{i3,1}.plotdat_preproc.im_c.new;
        %         dO_c    = days{i3,1}.plotdat_preproc.de_c.old;
        %         dN_c    = days{i3,1}.plotdat_preproc.de_c.new;
        
        iRO_c   = dat{i3,1}.plotdat_preproc.im_c.RO;
        iRN_c   = dat{i3,1}.plotdat_preproc.im_c.RN;
        %         dRO_c   = days{i3,1}.plotdat_preproc.de_c.RO;
        %         dRN_c   = days{i3,1}.plotdat_preproc.de_c.RN;
        
        inRO_c  = dat{i3,1}.plotdat_preproc.im_c.nRO;
        inRN_c  = dat{i3,1}.plotdat_preproc.im_c.nRN;
        %         dnRO_c  = days{i3,1}.plotdat_preproc.de_c.nRO;
        %         dnRN_c  = days{i3,1}.plotdat_preproc.de_c.nRN;
        
        
    else
        iR  = vertcat(iR,dat{i3,1}.plotdat_preproc.im.R);
        inR = vertcat(inR,dat{i3,1}.plotdat_preproc.im.nR);
        %         dR  = vertcat(dR,days{i3,1}.plotdat_preproc.de.R);
        %         dnR = vertcat(dnR,days{i3,1}.plotdat_preproc.de.nR);
        
        iO    = vertcat(iO,dat{i3,1}.plotdat_preproc.im.old);
        iN    = vertcat(iN,dat{i3,1}.plotdat_preproc.im.new);
        %         dO    = vertcat(dO,days{i3,1}.plotdat_preproc.de.old);
        %         dN    = vertcat(dN,days{i3,1}.plotdat_preproc.de.new);
        
        iRO   = vertcat(iRO,dat{i3,1}.plotdat_preproc.im.RO);
        iRN   = vertcat(iRN,dat{i3,1}.plotdat_preproc.im.RN);
        %         dRO   = vertcat(dRO,days{i3,1}.plotdat_preproc.de.RO);
        %         dRN   = vertcat(dRN,days{i3,1}.plotdat_preproc.de.RN);
        
        inRO  = vertcat(inRO,dat{i3,1}.plotdat_preproc.im.nRO);
        inRN  = vertcat(inRN,dat{i3,1}.plotdat_preproc.im.nRN);
        %         dnRO  = vertcat(dnRO,days{i3,1}.plotdat_preproc.de.nRO);
        %         dnRN  = vertcat(dnRN,days{i3,1}.plotdat_preproc.de.nRN);
        
        
        % correct trls
        iR_c    = vertcat(iR_c,dat{i3,1}.plotdat_preproc.im_c.R);
        inR_c   = vertcat(inR_c,dat{i3,1}.plotdat_preproc.im_c.nR);
        %         dR_c    = vertcat(dR_c,days{i3,1}.plotdat_preproc.de_c.R);
        %         dnR_c   = vertcat(dnR_c,days{i3,1}.plotdat_preproc.de_c.nR);
        
        iO_c    = vertcat(iO_c,dat{i3,1}.plotdat_preproc.im_c.old);
        iN_c    = vertcat(iN_c,dat{i3,1}.plotdat_preproc.im_c.new);
        %         dO_c    = vertcat(dO_c,days{i3,1}.plotdat_preproc.de_c.old);
        %         dN_c    = vertcat(dN_c,days{i3,1}.plotdat_preproc.de_c.new);
        
        iRO_c   = vertcat(iRO_c,dat{i3,1}.plotdat_preproc.im_c.RO);
        iRN_c   = vertcat(iRN_c,dat{i3,1}.plotdat_preproc.im_c.RN);
        %         dRO_c   = vertcat(dRO_c,days{i3,1}.plotdat_preproc.de_c.RO);
        %         dRN_c   = vertcat(dRN_c,days{i3,1}.plotdat_preproc.de_c.RN);
        
        inRO_c  = vertcat(inRO_c,dat{i3,1}.plotdat_preproc.im_c.nRO);
        inRN_c  = vertcat(inRN_c,dat{i3,1}.plotdat_preproc.im_c.nRN);
        %         dnRO_c  = vertcat(dnRO_c,days{i3,1}.plotdat_preproc.de_c.nRO);
        %         dnRN_c  = vertcat(dnRN_c,days{i3,1}.plotdat_preproc.de_c.nRN);
        
    end
    
end


cd(path_fig) % saving path for the figures

% fancier: shadedErrorBar
% rew vs. non-rew
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
% subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iR,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),inR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iR,inR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar
%
% subplot(1,2,2);
% plot(onsetbar,y,'g-','linewidth',2)
% shadedErrorBar((1:2701),dR,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),dnR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(dR,dnR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('delayed','Fontsize',20,'Fontweight','bold')
% suptitle('Reward(red) vs. Non-Reward(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_memtest_pilot'  '_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
clear fig H sigbar


% fancier: shadedErrorBar

% only correct trials: rew vs. non-rew
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
% subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iR_c,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),inR_c,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iR_c,inR_c); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar
%
% subplot(1,2,2);
% plot(onsetbar,y,'g-','linewidth',2)
%
% shadedErrorBar((1:2701),dR_c,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),dnR_c,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(dR_c,dnR_c); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('delayed','Fontsize',20,'Fontweight','bold')
% suptitle('Corrects: Reward(red) vs. Non-Reward(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_memtest_pilot'  '_CorrectTrials_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
clear fig H sigbar



% old vs. new
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
% subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),iN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iO,iN); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar

% subplot(1,2,2);
% plot(onsetbar,y,'g-','linewidth',2)
%
% shadedErrorBar((1:2701),dO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),dN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(dO,dN); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('delayed','Fontsize',20,'Fontweight','bold')
% suptitle('OLD(red) vs. NEW(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_memtest_pilot'  '_New_vs_Old_pilot'],'-dpdf','-r0')
clear fig H sigbar




% only no-reward trls: old v. new

maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
% subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),inRO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),inRN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(inRO,inRN); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar
% % feedback seg: reward vs. non-reward
% subplot(1,2,2);
% plot(onsetbar,y,'g-','linewidth',2)
%
% shadedErrorBar((1:2701),dnRO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
% shadedErrorBar((1:2701),dnRN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
% H = ttest2(dnRO,dnRN); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
% xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
% xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
% ylabel('zscores','Fontsize',20,'Fontweight','bold');
% title('delayed','Fontsize',20,'Fontweight','bold')
% suptitle('no-reward: OLD(red) vs. NEW(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_memtest_pilot'  '_noReward_New_vs_Old_pilot'],'-dpdf','-r0')
clear fig H sigbar


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
path_fig    = [ path_gen 'C_writings/CB_figures/pilot_3T/zscores/' ];
path_stim   = [ path_par 'EA_raw/EAA_pupil/stiminfo_all_20181004.mat' ];

% create file list and VP ID list
VPs = [1109 1210 1111 1212 1115];
session = 1;

%% now run

for i2 = 1:length(VPs)
    dat{i2,1} = load([ path_clean 'rr' num2str(session) num2str(VPs(i2)) '_preproc_dat.mat']);
end

iR=[];inR=[];dR=[];dnR=[];
for i3 = 1:size(dat,1) % number of IDs
    if isempty(iR)
        iR      = dat{i3,1}.plotdat_preproc.im.R;
        inR     = dat{i3,1}.plotdat_preproc.im.nR;
        dR      = dat{i3,1}.plotdat_preproc.de.R(:,1:2701);
        dnR     = dat{i3,1}.plotdat_preproc.de.nR(:,1:2701);
        
        iO      = dat{i3,1}.plotdat_preproc.im.old;
        iN      = dat{i3,1}.plotdat_preproc.im.new;
        dO      = dat{i3,1}.plotdat_preproc.de.old(:,1:2701);
        dN      = dat{i3,1}.plotdat_preproc.de.new(:,1:2701);
        
        iRO     = dat{i3,1}.plotdat_preproc.im.RO;
        iRN     = dat{i3,1}.plotdat_preproc.im.RN;
        dRO     = dat{i3,1}.plotdat_preproc.de.RO(:,1:2701);
        dRN     = dat{i3,1}.plotdat_preproc.de.RN(:,1:2701);
        
        inRO    = dat{i3,1}.plotdat_preproc.im.nRO;
        inRN    = dat{i3,1}.plotdat_preproc.im.nRN;
        dnRO    = dat{i3,1}.plotdat_preproc.de.nRO(:,1:2701);
        dnRN    = dat{i3,1}.plotdat_preproc.de.nRN(:,1:2701);
        
        %         correct trls
        iR_c    = dat{i3,1}.plotdat_preproc.im_c.R;
        inR_c   = dat{i3,1}.plotdat_preproc.im_c.nR;
        dR_c    = dat{i3,1}.plotdat_preproc.de_c.R(:,1:2701);
        dnR_c   = dat{i3,1}.plotdat_preproc.de_c.nR(:,1:2701);
        
        iO_c    = dat{i3,1}.plotdat_preproc.im_c.old;
        iN_c    = dat{i3,1}.plotdat_preproc.im_c.new;
        dO_c    = dat{i3,1}.plotdat_preproc.de_c.old(:,1:2701);
        dN_c    = dat{i3,1}.plotdat_preproc.de_c.new(:,1:2701);
        
        iRO_c   = dat{i3,1}.plotdat_preproc.im_c.RO;
        iRN_c   = dat{i3,1}.plotdat_preproc.im_c.RN;
        dRO_c   = dat{i3,1}.plotdat_preproc.de_c.RO(:,1:2701);
        dRN_c   = dat{i3,1}.plotdat_preproc.de_c.RN(:,1:2701);
        
        inRO_c  = dat{i3,1}.plotdat_preproc.im_c.nRO;
        inRN_c  = dat{i3,1}.plotdat_preproc.im_c.nRN;
        dnRO_c  = dat{i3,1}.plotdat_preproc.de_c.nRO(:,1:2701);
        dnRN_c  = dat{i3,1}.plotdat_preproc.de_c.nRN(:,1:2701);
        
        
    else
        iR  = vertcat(iR,dat{i3,1}.plotdat_preproc.im.R);
        inR = vertcat(inR,dat{i3,1}.plotdat_preproc.im.nR);
        dR  = vertcat(dR,dat{i3,1}.plotdat_preproc.de.R(:,1:2701));
        dnR = vertcat(dnR,dat{i3,1}.plotdat_preproc.de.nR(:,1:2701));
        
        iO    = vertcat(iO,dat{i3,1}.plotdat_preproc.im.old);
        iN    = vertcat(iN,dat{i3,1}.plotdat_preproc.im.new);
        dO    = vertcat(dO,dat{i3,1}.plotdat_preproc.de.old(:,1:2701));
        dN    = vertcat(dN,dat{i3,1}.plotdat_preproc.de.new(:,1:2701));
        
        iRO   = vertcat(iRO,dat{i3,1}.plotdat_preproc.im.RO);
        iRN   = vertcat(iRN,dat{i3,1}.plotdat_preproc.im.RN);
        dRO   = vertcat(dRO,dat{i3,1}.plotdat_preproc.de.RO(:,1:2701));
        dRN   = vertcat(dRN,dat{i3,1}.plotdat_preproc.de.RN(:,1:2701));
        
        inRO  = vertcat(inRO,dat{i3,1}.plotdat_preproc.im.nRO);
        inRN  = vertcat(inRN,dat{i3,1}.plotdat_preproc.im.nRN);
        dnRO  = vertcat(dnRO,dat{i3,1}.plotdat_preproc.de.nRO(:,1:2701));
        dnRN  = vertcat(dnRN,dat{i3,1}.plotdat_preproc.de.nRN(:,1:2701));
        
        
        %         correct trls
        iR_c    = vertcat(iR_c,dat{i3,1}.plotdat_preproc.im_c.R);
        inR_c   = vertcat(inR_c,dat{i3,1}.plotdat_preproc.im_c.nR);
        dR_c    = vertcat(dR_c,dat{i3,1}.plotdat_preproc.de_c.R(:,1:2701));
        dnR_c   = vertcat(dnR_c,dat{i3,1}.plotdat_preproc.de_c.nR(:,1:2701));
        
        iO_c    = vertcat(iO_c,dat{i3,1}.plotdat_preproc.im_c.old);
        iN_c    = vertcat(iN_c,dat{i3,1}.plotdat_preproc.im_c.new);
        dO_c    = vertcat(dO_c,dat{i3,1}.plotdat_preproc.de_c.old(:,1:2701));
        dN_c    = vertcat(dN_c,dat{i3,1}.plotdat_preproc.de_c.new(:,1:2701));
        
        iRO_c   = vertcat(iRO_c,dat{i3,1}.plotdat_preproc.im_c.RO);
        iRN_c   = vertcat(iRN_c,dat{i3,1}.plotdat_preproc.im_c.RN);
        dRO_c   = vertcat(dRO_c,dat{i3,1}.plotdat_preproc.de_c.RO(:,1:2701));
        dRN_c   = vertcat(dRN_c,dat{i3,1}.plotdat_preproc.de_c.RN(:,1:2701));
        
        inRO_c  = vertcat(inRO_c,dat{i3,1}.plotdat_preproc.im_c.nRO);
        inRN_c  = vertcat(inRN_c,dat{i3,1}.plotdat_preproc.im_c.nRN);
        dnRO_c  = vertcat(dnRO_c,dat{i3,1}.plotdat_preproc.de_c.nRO(:,1:2701));
        dnRN_c  = vertcat(dnRN_c,dat{i3,1}.plotdat_preproc.de_c.nRN(:,1:2701));
        
    end
    
end


cd(path_fig) % saving path for the figures

% fancier: shadedErrorBar
% rew vs. non-rew
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iR,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),inR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iR,inR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)
shadedErrorBar((1:2701),dR,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),dnR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(dR,dnR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('delayed','Fontsize',20,'Fontweight','bold')
suptitle('Reward(red) vs. Non-Reward(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_memtest_pilot'  '_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
clear fig H sigbar


% fancier: shadedErrorBar

% only correct trials: rew vs. non-rew
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iR_c,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),inR_c,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iR_c,inR_c); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),dR_c,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),dnR_c,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(dR_c,dnR_c); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('delayed','Fontsize',20,'Fontweight','bold')
suptitle('Corrects: Reward(red) vs. Non-Reward(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_memtest_pilot'  '_CorrectTrials_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
clear fig H sigbar



% old vs. new
maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),iO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),iN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(iO,iN); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar

subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),dO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),dN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(dO,dN); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('delayed','Fontsize',20,'Fontweight','bold')
suptitle('OLD(red) vs. NEW(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_memtest_pilot'  '_New_vs_Old_pilot'],'-dpdf','-r0')
clear fig H sigbar




% only no-reward trls: old v. new

maxy = 3; miny = -3;
onsetbar=200.*ones(1,2701); y=linspace(miny,maxy,numel(onsetbar));
errorbar = @ (data) (std(data)/sqrt(length(VPs)));

figure;
subplot(1,2,1);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),inRO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),inRN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(inRO,inRN); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('immediate','Fontsize',20,'Fontweight','bold')
clear H sigbar
% feedback seg: reward vs. non-reward
subplot(1,2,2);
plot(onsetbar,y,'g-','linewidth',2)

shadedErrorBar((1:2701),dnRO,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
shadedErrorBar((1:2701),dnRN,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
H = ttest2(dnRO,dnRN); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; sigbar(1:200)=NaN; plot(sigbar,'k-o','LineWidth',3);
xlim([1 2701]); ylim([miny maxy]); %xticklabels([0 500 1000 1500 2000 2500]);
xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
ylabel('zscores','Fontsize',20,'Fontweight','bold');
title('delayed','Fontsize',20,'Fontweight','bold')
suptitle('no-reward: OLD(red) vs. NEW(blue)')

% save
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['GroupLvl_memtest_pilot'  '_noReward_New_vs_Old_pilot'],'-dpdf','-r0')
clear fig H sigbar
