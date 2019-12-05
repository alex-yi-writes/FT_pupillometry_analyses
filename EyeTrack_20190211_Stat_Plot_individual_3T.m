%% Pupillometric data statistics and plotting for 3T pilot : individuals

%% work log

%   11_02_2019  created the script

%% index (search keywords)
%  maintask
%  oddball

%% main task

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
VPs = [1204 1206 1107 1208 1109 1210 1111 1212 1214 1216]; % 1105, 1113 excluded for bad data
% VPs = 1109;

%% plot and create masks

fname_raweye = []; fname_beh =[]; expdat = [];
for i1 = 1:length(VPs)
        
        fname_raweye{i1,1}  = [ 're1' num2str(VPs(i1)) '.asc' ];
        fname_beh{i1,1}     = [ num2str(VPs(i1)) '.mat' ];
                
        expdat{i1,1} = load([path_behav fname_beh{i1,1}]);

        % call in cleaned data
        plotdat_stim  = load([path_clean '/re1' num2str(VPs(i1)) '_eye_rv200_stim.mat']);
        
        fprintf('\n Preparation done \n')
        
        %% housekeeping
        
        % clean the null trials and rejected trials: breaks
        rmrow_stim   = find(plotdat_stim.eye_rv.trl_raus==1); % rejected trials from preprocessing
        breaktrials  = [30 60 90 120 150];
        emptymeans   = cellfun(@nanmean,plotdat_stim.eye_rv.trial);
        emptytrials  = find(emptymeans==0);
                
        % stimulus informations: filename, reward/no-reward, memorability
        stimdat      = expdat{i1,1}.dat.day1.maintask.results.trl;
        stimdat(isemptycell(stimdat)) = {NaN}; % this line, because if you convert data directly from cells to double, it automatically removes empty cells, which in our case are breaks.
        
        % remove breaks and rejected trials
        stimdat(breaktrials,:) = [];
        stimdat(rmrow_stim,:) = []; % and now, we remove the rows of trials we rejected during data cleaning
        stimdat(emptytrials,:) = [];
        
        % indexing
        rewcond  = str2num(expdat{i1,1}.dat.RewardCategory(9)); % rewarding condition, human, was indexed as 2 in our data structure
        if rewcond==1; puncond = 2; else; puncond = 1; end
        
        ind_rew_stim    = cell2mat(stimdat(:,2)) == rewcond; % marking the trials that are rewards, and so on
        ind_pun_stim    = cell2mat(stimdat(:,2)) == puncond;
        ind_hiM_stim    = cell2mat(stimdat(:,3)) == 11 | cell2mat(stimdat(:,3)) == 21;
        ind_loM_stim    = cell2mat(stimdat(:,3)) == 12 | cell2mat(stimdat(:,3)) == 22;
        ind_hiM_rew_stim= cell2mat(stimdat(:,3)) == 11;
        ind_loM_rew_stim= cell2mat(stimdat(:,3)) == 12;
        ind_hiM_pun_stim= cell2mat(stimdat(:,3)) == 21;
        ind_loM_pun_stim= cell2mat(stimdat(:,3)) == 22;
        
        % assign data: now we assign our actual pupill time-series data
        stim_eye        = plotdat_stim.eye_rv.trial'; % raw, all survived trial data
        stim_rew        = stim_eye(ind_rew_stim,1); % only reward trials
        stim_pun        = stim_eye(ind_pun_stim,1); % only non-reward trials
        stim_hiM        = stim_eye(ind_hiM_stim,1);
        stim_loM        = stim_eye(ind_loM_stim,1);
        stim_hiM_rew    = stim_eye(ind_hiM_rew_stim,1);
        stim_loM_rew    = stim_eye(ind_loM_rew_stim,1);
        stim_hiM_pun    = stim_eye(ind_hiM_pun_stim,1);
        stim_loM_pun    = stim_eye(ind_loM_pun_stim,1);
        
        % gather and standardise
        time = plotdat_stim.eye_rv.time;
        timewindow_bsl = [1 200]; % -200ms ~ 0ms
        plotdat        = [];
        
        % first, we squeeze each trial and time-series into one matrix for the sake of convenience
        
        plotdat.stim.all = []; % stim
        for i = 1:size(stim_eye,1) % for the number of trials
            plotdat.stim.all(i,:) = stim_eye{i,:}; % create one big matrix
        end        
        plotdat.stim.all(emptytrials,:) = [];
        
        % stat preprocessing fix 19-07-2018:
        % first, standardise within trials, then baseline correct.
        plotdat_z = []; clear i
        for i = 1:size(plotdat.stim.all,1)
            meandat = nanmean(plotdat.stim.all,2); stddat = nanstd(plotdat.stim.all');
            plotdat_z.stim.all = plotdat.stim.all - repmat(meandat,1,size(plotdat.stim.all,2));
            plotdat_z.stim.all = plotdat.stim.all ./ repmat(stddat',1,size(plotdat.stim.all,2));
            
            bsldat = nanmean(plotdat_z.stim.all(:,timewindow_bsl(1):timewindow_bsl(2))');
            plotdat_preproc.stim.all = plotdat_z.stim.all-repmat(bsldat',1,size(plotdat_z.stim.all,2));
        end
        
        plotdat_preproc.stim.rew     = plotdat_preproc.stim.all(ind_rew_stim,:);
        plotdat_preproc.stim.pun     = plotdat_preproc.stim.all(ind_pun_stim,:);
        plotdat_preproc.stim.hiMem   = plotdat_preproc.stim.all(ind_hiM_stim,:);
        plotdat_preproc.stim.loMem   = plotdat_preproc.stim.all(ind_loM_stim,:);
        plotdat_preproc.stim.hiMem_rew = plotdat_preproc.stim.all(ind_hiM_rew_stim,:);
        plotdat_preproc.stim.loMem_rew = plotdat_preproc.stim.all(ind_loM_rew_stim,:);
        plotdat_preproc.stim.hiMem_pun = plotdat_preproc.stim.all(ind_hiM_pun_stim,:);
        plotdat_preproc.stim.loMem_pun = plotdat_preproc.stim.all(ind_loM_pun_stim,:);
                       
        % average over the corresponding trials
        plotavg.stim.all     = mean(plotdat_preproc.stim.all,1);
        plotavg.stim.rew     = mean(plotdat_preproc.stim.rew,1);
        plotavg.stim.pun     = mean(plotdat_preproc.stim.pun,1);
        plotavg.stim.hiMem   = mean(plotdat_preproc.stim.all(ind_hiM_stim,:));
        plotavg.stim.loMem   = mean(plotdat_preproc.stim.all(ind_loM_stim,:));
        plotavg.stim.hiMem_rew = mean(plotdat_preproc.stim.all(ind_hiM_rew_stim,:));
        plotavg.stim.loMem_rew = mean(plotdat_preproc.stim.all(ind_loM_rew_stim,:));
        plotavg.stim.hiMem_pun = mean(plotdat_preproc.stim.all(ind_hiM_pun_stim,:));
        plotavg.stim.loMem_pun = mean(plotdat_preproc.stim.all(ind_loM_pun_stim,:));
                
        save([ path_clean 're1' num2str(VPs(i1)) '_preproc_dat.mat'], 'plotdat_preproc','plotavg');
        
        
        fprintf('\n Done \n')
        
        %% draw
        % legends do not work -- some errors with cvx, i think
        
        cd(path_fig) % saving path for the figures
        close all
        
        maxy = 3; miny = -3;
        onsetbar=200.*ones(1,length(time{1,1})); y=linspace(miny,maxy,numel(onsetbar));
        errorbar = @ (data) (std(data)/sqrt(length(VPs)));
        
        % fancier: shadedErrorBar
        % if you would like to know about what to put in as arguments,
        % please type: doc shadedErrorBar
        
        % stimulus seg: reward vs. non-reward
        figure;
        plot(onsetbar,y,'g-','linewidth',2)
        shadedErrorBar((1:length(time{1,1})),plotdat_preproc.stim.rew,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
        shadedErrorBar((1:length(time{1,1})),plotdat_preproc.stim.pun,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
        H1 = ttest2(plotdat_preproc.stim.rew,plotdat_preproc.stim.pun); sigbar1 = H1*-0.9; sigbar1(find(sigbar1==0))=NaN; plot(sigbar1,'k-o','LineWidth',3); % this bit is for calculating significance
        xlim([1 length(time{1,1})]); ylim([miny maxy]);% xticklabels([0 500 1000 1500 2000 2500]);
        xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
        ylabel('zscores','Fontsize',20,'Fontweight','bold');
        title('reward vs. non-reward','Fontsize',20,'Fontweight','bold')
        clear H sigbar
        
        % save
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print([num2str(VPs(i1)) '_maintask_Reward(red)_vs_Punishment(blue)'],'-dpdf','-r0')
        clear fig H sigbar
        close all
        
        
        
        fprintf('\n ID %d: PLOTTING DONE\n', VPs(i1))
        
        keep session memtest VPs ...
            path_gen path_par path_dat path_behav path_clean path_fig path_stim ...
            fname_raweye fname_beh expdat i1 onsetbar maxy miny time
        
end

fprintf('\n main task plotting done\n')

%% oddball task

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
VPs = [1204 1206 1107 1109 1210 1111 1212 1113 1214 1216];

% set enviornmental variables
load('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBB_eyetracking/nulltrials_OB.mat')


%% plot and create masks

fname_raweye = []; fname_beh =[]; expdat = [];
for i1 = 1:length(VPs)
        
        fname_raweye{i1,1}  = [ 'ob' num2str(VPs(i1)) '.asc' ];
        fname_beh{i1,1}     = [ num2str(VPs(i1)) '_OB.mat' ];
        
        expdat{i1,1} = load([path_behav fname_beh{i1,1}]);
        
        % call in cleaned data
        plotdat_stim  = load([path_clean '/ob' num2str(VPs(i1)) '_eye_rv200_stim.mat']);
        
        fprintf('\n Preparation done \n')
        
        %% housekeeping
        
        % clean the null trials and rejected trials: breaks
        trl_raus = plotdat_stim.eye_rv.trl_raus;
        null_trls = nulltrials(:,expdat{i1,1}.dat.oddball.StimSchedule);            
        rmrow_stim   = find(trl_raus==1); % rejected trials from preprocessing
        
        % stimulus informations: filename, reward/no-reward, memorability
        stimdat      = expdat{i1,1}.dat.oddball.results.trl;
        stimdat(isemptycell(stimdat)) = {NaN}; % this line, because if you convert data directly from cells to double, it automatically removes empty cells, which in our case are breaks.
        indices = cell2mat(stimdat(:,2));
        
        % remove breaks and rejected trials
        stimdat(null_trls,:) = [];
        stimdat(rmrow_stim,:) = []; % and now, we remove the rows of trials we rejected during data cleaning
        
        % indexing
        oddcond = 1; stdcond = 2; % nullcond = 0;
        
        ind_odd_stim    = cell2mat(stimdat(:,2)) == oddcond; % marking the trials that are rewards, and so on
        ind_std_stim    = cell2mat(stimdat(:,2)) == stdcond;
%         ind_null_stim   = cell2mat(stimdat(:,2)) == nullcond;
                
        % assign data: now we assign our actual pupill time-series data
        stim_eye        = plotdat_stim.eye_rv.trial'; % raw, all survived trial data
        stim_odd        = stim_eye(ind_odd_stim,1); % only reward trials
        stim_std        = stim_eye(ind_std_stim,1); % only non-reward trials
                
        % gather and standardise
        time = plotdat_stim.eye_rv.time;
        timewindow_bsl = [1 200]; % -200ms ~ 0ms
        plotdat        = [];
        
        % first, we squeeze each trial and time-series into one matrix for the sake of convenience
        
        plotdat.stim.all = []; % stim
        for i = 1:size(stim_eye,1) % for the number of trials
            plotdat.stim.all(i,:) = stim_eye{i,:}; % create one big matrix
        end
        
        % stat preprocessing fix 19-07-2018:
        % first, standardise within trials, then baseline correct.
        plotdat_z = []; clear i
        for i = 1:size(plotdat.stim.all,1)
            meandat = nanmean(plotdat.stim.all,2); stddat = nanstd(plotdat.stim.all');
            plotdat_z.stim.all = plotdat.stim.all - repmat(meandat,1,size(plotdat.stim.all,2));
            plotdat_z.stim.all = plotdat.stim.all ./ repmat(stddat',1,size(plotdat.stim.all,2));
            
            bsldat = nanmean(plotdat_z.stim.all(:,timewindow_bsl(1):timewindow_bsl(2))');
            plotdat_preproc.stim.all = plotdat_z.stim.all-repmat(bsldat',1,size(plotdat_z.stim.all,2));
        end
        
        plotdat_preproc.stim.odd     = plotdat_preproc.stim.all(ind_odd_stim,:);
        plotdat_preproc.stim.std     = plotdat_preproc.stim.all(ind_std_stim,:);
%         plotdat_preproc.stim.null    = plotdat_preproc.stim.all(ind_null_stim,:);
                               
        % average over the corresponding trials
        plotavg.stim.all     = mean(plotdat_preproc.stim.all,1);
        plotavg.stim.odd     = mean(plotdat_preproc.stim.odd,1);
        plotavg.stim.std     = mean(plotdat_preproc.stim.std,1);
%         plotavg.stim.null    = mean(plotdat_preproc.stim.null,1);
                        
        save([ path_clean 'ob' num2str(VPs(i1)) '_preproc_dat.mat'], 'plotdat_preproc','plotavg');
        
        
        fprintf('\n Done \n')
        
        %% draw
        % legends do not work -- some errors with cvx, i think
        
        cd(path_fig) % saving path for the figures
        close all
        
        maxy = 3; miny = -3;
        onsetbar=200.*ones(1,length(time{1,1})); y=linspace(miny,maxy,numel(onsetbar));
        errorbar = @ (data) (std(data)/sqrt(length(VPs)));
        
        % fancier: shadedErrorBar
        % if you would like to know about what to put in as arguments,
        % please type: doc shadedErrorBar
        
        % stimulus seg: reward vs. non-reward
        figure;
        plot(onsetbar,y,'g-','linewidth',2)
        shadedErrorBar((1:length(time{1,1})),plotdat_preproc.stim.odd,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
        shadedErrorBar((1:length(time{1,1})),plotdat_preproc.stim.std,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%         shadedErrorBar((1:length(time{1,1})),plotdat_preproc.stim.null,{@nanmean,errorbar},'lineprops',{'y-o','markerfacecolor','y'}); hold on;
        H1 = ttest2(plotdat_preproc.stim.odd,plotdat_preproc.stim.std); sigbar1 = H1*-0.9; sigbar1(find(sigbar1==0))=NaN; plot(sigbar1,'k-o','LineWidth',3); % this bit is for calculating significance
%         H2 = ttest2(plotdat_preproc.stim.null,plotdat_preproc.stim.std); sigbar2 = H2*-0.8; sigbar2(find(sigbar2==0))=NaN; plot(sigbar2,'g-o','LineWidth',3); % this bit is for calculating significance
        xlim([1 length(time{1,1})]); ylim([miny maxy]);% xticklabels([0 500 1000 1500 2000 2500]);
        xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
        ylabel('zscores','Fontsize',20,'Fontweight','bold');
        title('Oddball(red) vs. Standard(blue)','Fontsize',20,'Fontweight','bold')
        clear H sigbar
        
        % save
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print([num2str(VPs(i1)) '_OB_Oddball(red)_vs_Standard(blue)'],'-dpdf','-r0')
        clear fig H1 H2 sigbar1 sigbar2
        close all
        
        
        
        fprintf('\n ID %d: PLOTTING DONE\n', VPs(i1))
        
        keep session memtest VPs ...
            path_gen path_par path_dat path_behav path_clean path_fig path_stim ...
            fname_raweye fname_beh expdat i1 onsetbar maxy miny time nulltrials
        
end

fprintf('\n main task plotting done\n')


%% immediate memory test
%  from this line on, we'll analyse memory data

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
% VPs = [1204 1107 1208 1109 1210];
% VPs = [1204 1107 1208 1109 1111 1212 1214 1216];
VPs = 1115;
session = 1;

%% plot and create masks

fname_raweye = []; fname_beh =[]; expdat = [];
for i1 = 1:length(VPs)
    
    
    fname_raweye{i1,1}  = [ 'ri' num2str(1) num2str(VPs(i1)) '.asc' ];
%     fname_raweye{i1,1}  = [ 'rd' num2str(1) num2str(VPs(i1)) '.asc' ];
    fname_beh{i1,1}     = [ num2str(VPs(i1)) '.mat' ];
    
    expdat{i1,1} = load([path_behav fname_beh{i1,1}]);
    
    
    % call in cleaned data
    plotdat_im = load([ path_clean '/ri' num2str(session) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
%     plotdat_de = load([ path_clean '/rd' num2str(session) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
    
    fprintf('\n Preparation done \n')
    
    %% housekeeping
    
    
    % clean the null trials and rejected trials: breaks at trial 61, 122
    if VPs(i1) == 1204
        plotdat_im.eye_rv.immediate.trl_raus(1) = [];
    end
    
    rmrow_im   = find(plotdat_im.eye_rv.immediate.trl_raus==1);
%     rmrow_de   = find(plotdat_de.eye_rv.delayed.trl_raus==1);
    
    imdat = expdat{i1,1}.dat.day1.memorytest.immediate.config.stim.stimlist_all(:,[2 4]);
%     dedat = expdat{i1,1}.dat.day1.memorytest.delayed.config.stim.stimlist_all(:,[2 4]);    

    behdat_i= expdat{i1,1}.dat.day1.memorytest.immediate.results.all;
%     behdat_d= expdat{i1,1}.dat.day1.memorytest.delayed.results.all;
    
    imdat(rmrow_im,:) = [];
%     dedat(rmrow_de,:) = [];
    behdat_i(rmrow_im,:) = [];
%     behdat_d(rmrow_de,:) = [];
    
    % indexing: session 1
    
    % all trials
    rewcond  = str2num(expdat{i1,1}.dat.RewardCategory(9)); % rewarding condition, human, was indexed as 2 in our data structure
    if rewcond==1; puncond = 2; else; puncond = 1; end
    
    ind_Rew_im      = cell2mat(imdat(:,1) == rewcond);
    ind_nRew_im     = cell2mat(imdat(:,1) == puncond);
%     ind_Rew_de      = cell2mat(dedat(:,1) == rewcond);
%     ind_nRew_de     = cell2mat(dedat(:,1) == puncond);
    
    ind_old_im      = cell2mat(imdat(:,2) == 1);
    ind_new_im      = cell2mat(imdat(:,2) == 2);
%     ind_old_de      = cell2mat(dedat(:,3) == 1);
%     ind_new_de      = cell2mat(dedat(:,3) == 2);
    
    ind_RewO_im     = cell2mat((imdat(:,1) == rewcond) & imdat(:,2) == 1);
    ind_RewN_im     = cell2mat((imdat(:,1) == rewcond) & imdat(:,2) == 2);
    ind_nRewO_im    = cell2mat((imdat(:,1) == puncond) & imdat(:,2) == 1);
    ind_nRewN_im    = cell2mat((imdat(:,1) == puncond) & imdat(:,2) == 2);
    
    
    % correct/incorrect sorted
    
    % the corrects!
    ind_Rew_im_c    = ind_Rew_im == 1 & behdat_i(:,2) == 1;
    ind_nRew_im_c   = ind_nRew_im == 1 & behdat_i(:,2) == 1;
%     ind_Rew_de_c    = ind_Rew_de == 1 & behdat_d(:,2) == 1;
%     ind_nRew_de_c   = ind_nRew_de == 1 & behdat_d(:,2) == 1;
      
    ind_old_im_c    = ind_old_im == 1 & behdat_i(:,2) == 1;
    ind_new_im_c    = ind_new_im == 1 & behdat_i(:,2) == 1;
%     ind_old_de_c    = ind_old_de == 1 & behdat_d(:,2) == 1;
%     ind_new_de_c    = ind_new_de == 1 & behdat_d(:,2) == 1;
    
    ind_RewO_im_c   = ind_RewO_im ==1 & behdat_i(:,2) == 1;
    ind_RewN_im_c   = ind_RewN_im ==1 & behdat_i(:,2) == 1;
%     ind_RewO_de_c   = ind_RewO_de ==1 & behdat_d(:,2) == 1;
%     ind_RewN_de_c   = ind_RewN_de ==1 & behdat_d(:,2) == 1;
    
    ind_nRewO_im_c  = ind_nRewO_im ==1 & behdat_i(:,2) == 1;
    ind_nRewN_im_c  = ind_nRewN_im ==1 & behdat_i(:,2) == 1;
%     ind_nRewO_de_c  = ind_nRewO_de ==1 & behdat_d(:,2) == 1;
%     ind_nRewN_de_c  = ind_nRewN_de ==1 & behdat_d(:,2) == 1;
    
    
    % the incorrects!
    ind_Rew_im_ic    = ind_Rew_im == 1 & behdat_i(:,2) == 0;
    ind_nRew_im_ic   = ind_nRew_im == 1 & behdat_i(:,2) == 0;
%     ind_Rew_de_ic    = ind_Rew_de == 1 & behdat_d(:,2) == 0;
%     ind_nRew_de_ic   = ind_nRew_de == 1 & behdat_d(:,2) == 0;
        
    ind_old_im_ic    = ind_old_im == 1 & behdat_i(:,2) == 0;
    ind_new_im_ic    = ind_new_im == 1 & behdat_i(:,2) == 0;
%     ind_old_de_ic    = ind_old_de == 1 & behdat_d(:,2) == 0;
%     ind_new_de_ic    = ind_new_de == 1 & behdat_d(:,2) == 0;
    
    ind_RewO_im_ic   = ind_RewO_im ==1 & behdat_i(:,2) == 0;
    ind_RewN_im_ic   = ind_RewN_im ==1 & behdat_i(:,2) == 0;
%     ind_RewO_de_ic   = ind_RewO_de ==1 & behdat_d(:,2) == 0;
%     ind_RewN_de_ic   = ind_RewN_de ==1 & behdat_d(:,2) == 0;
    
    ind_nRewO_im_ic  = ind_nRewO_im ==1 & behdat_i(:,2) == 0;
    ind_nRewN_im_ic  = ind_nRewN_im ==1 & behdat_i(:,2) == 0;
%     ind_nRewO_de_ic  = ind_nRewO_de ==1 & behdat_d(:,2) == 0;
%     ind_nRewN_de_ic  = ind_nRewN_de ==1 & behdat_d(:,2) == 0;
    
    % assign data
    im_eye          = plotdat_im.eye_rv.trial';
    im_Rew          = im_eye(ind_Rew_im,1);
    im_nRew         = im_eye(ind_nRew_im,1);
    im_old          = im_eye(ind_old_im,1);
    im_new          = im_eye(ind_new_im,1);
    
%     de_eye          = plotdat_de.eye_rv.trial';
%     de_Rew          = de_eye(ind_Rew_de,1);
%     de_nRew         = de_eye(ind_nRew_de,1);
%     de_old          = de_eye(ind_old_de,1);
%     de_new          = de_eye(ind_new_de,1);
    
    
    % gather and standardise
    time = plotdat_im.eye_rv.time;
    timewindow_bsl = [1 200]; % -200ms ~ 0ms
    plotdat        = [];
    
    plotdat.im.all = []; % immediate
    for i = 1:size(im_eye,1) % for the number of trials
        plotdat.im.all(i,:) = im_eye{i,:}; % create one big matrix
    end
%     plotdat.de.all = []; % delayed
%     for i = 1:size(de_eye,1)
%         plotdat.de.all(i,:) = de_eye{i,:};
%     end
    
    % immediate
    plotdat_z = []; clear i
    for i = 1:size(plotdat.im.all,1)
        meandat = nanmean(plotdat.im.all,2); stddat = nanstd(plotdat.im.all');
        plotdat_z.im.all = plotdat.im.all - repmat(meandat,1,size(plotdat.im.all,2));
        plotdat_z.im.all = plotdat.im.all ./ repmat(stddat',1,size(plotdat.im.all,2));
        
        bsldat = nanmean(plotdat_z.im.all(:,timewindow_bsl(1):timewindow_bsl(2))');
        plotdat_preproc.im.all = plotdat_z.im.all-repmat(bsldat',1,size(plotdat_z.im.all,2));
    end
    plotdat_preproc.im.R      = plotdat_preproc.im.all(ind_Rew_im,:);
    plotdat_preproc.im.nR     = plotdat_preproc.im.all(ind_nRew_im,:);
    plotdat_preproc.im.old    = plotdat_preproc.im.all(ind_old_im,:);
    plotdat_preproc.im.new    = plotdat_preproc.im.all(ind_new_im,:);
    plotdat_preproc.im.RO     = plotdat_preproc.im.all(ind_RewO_im,:);
    plotdat_preproc.im.RN     = plotdat_preproc.im.all(ind_RewN_im,:);
    plotdat_preproc.im.nRO    = plotdat_preproc.im.all(ind_nRewO_im,:);
    plotdat_preproc.im.nRN    = plotdat_preproc.im.all(ind_nRewN_im,:);
    
    plotdat_preproc.im_c.R      = plotdat_preproc.im.all(ind_Rew_im_c,:);
    plotdat_preproc.im_c.nR     = plotdat_preproc.im.all(ind_nRew_im_c,:);
    plotdat_preproc.im_c.old    = plotdat_preproc.im.all(ind_old_im_c,:);
    plotdat_preproc.im_c.new    = plotdat_preproc.im.all(ind_new_im_c,:);
    plotdat_preproc.im_c.RO     = plotdat_preproc.im.all(ind_RewO_im_c,:);
    plotdat_preproc.im_c.RN     = plotdat_preproc.im.all(ind_RewN_im_c,:);
    plotdat_preproc.im_c.nRO    = plotdat_preproc.im.all(ind_nRewO_im_c,:);
    plotdat_preproc.im_c.nRN    = plotdat_preproc.im.all(ind_nRewN_im_c,:);
    
    plotdat_preproc.im_ic.R      = plotdat_preproc.im.all(ind_Rew_im_ic,:);
    plotdat_preproc.im_ic.nR     = plotdat_preproc.im.all(ind_nRew_im_ic,:);
    plotdat_preproc.im_ic.old    = plotdat_preproc.im.all(ind_old_im_ic,:);
    plotdat_preproc.im_ic.new    = plotdat_preproc.im.all(ind_new_im_ic,:);
    plotdat_preproc.im_ic.RO     = plotdat_preproc.im.all(ind_RewO_im_ic,:);
    plotdat_preproc.im_ic.RN     = plotdat_preproc.im.all(ind_RewN_im_ic,:);
    plotdat_preproc.im_ic.nRO    = plotdat_preproc.im.all(ind_nRewO_im_ic,:);
    plotdat_preproc.im_ic.nRN    = plotdat_preproc.im.all(ind_nRewN_im_ic,:);
    
%     % delayed
%     plotdat_z = []; clear i
%     for i = 1:size(plotdat.de.all,1)
%         meandat = nanmean(plotdat.de.all,2); stddat = nanstd(plotdat.de.all');
%         plotdat_z.de.all = plotdat.de.all - repmat(meandat,1,size(plotdat.de.all,2));
%         plotdat_z.de.all = plotdat.de.all ./ repmat(stddat',1,size(plotdat.de.all,2));
%         
%         bsldat = nanmean(plotdat_z.de.all(:,timewindow_bsl(1):timewindow_bsl(2))');
%         plotdat_preproc.de.all = plotdat_z.de.all-repmat(bsldat',1,size(plotdat_z.de.all,2));
%     end
%     plotdat_preproc.de.R      = plotdat_preproc.de.all(ind_Rew_de,:);
%     plotdat_preproc.de.nR     = plotdat_preproc.de.all(ind_nRew_de,:);
%     plotdat_preproc.de.old    = plotdat_preproc.de.all(ind_old_de,:);
%     plotdat_preproc.de.new    = plotdat_preproc.de.all(ind_new_de,:);
%     plotdat_preproc.de.RO     = plotdat_preproc.de.all(ind_RewO_de,:);
%     plotdat_preproc.de.RN     = plotdat_preproc.de.all(ind_RewN_de,:);
%     plotdat_preproc.de.nRO    = plotdat_preproc.de.all(ind_nRewO_de,:);
%     plotdat_preproc.de.nRN    = plotdat_preproc.de.all(ind_nRewN_de,:);
%     
%     plotdat_preproc.de_c.R      = plotdat_preproc.de.all(ind_Rew_de_c,:);
%     plotdat_preproc.de_c.nR     = plotdat_preproc.de.all(ind_nRew_de_c,:);
%     plotdat_preproc.de_c.old    = plotdat_preproc.de.all(ind_old_de_c,:);
%     plotdat_preproc.de_c.new    = plotdat_preproc.de.all(ind_new_de_c,:);
%     plotdat_preproc.de_c.RO     = plotdat_preproc.de.all(ind_RewO_de_c,:);
%     plotdat_preproc.de_c.RN     = plotdat_preproc.de.all(ind_RewN_de_c,:);
%     plotdat_preproc.de_c.nRO    = plotdat_preproc.de.all(ind_nRewO_de_c,:);
%     plotdat_preproc.de_c.nRN    = plotdat_preproc.de.all(ind_nRewN_de_c,:);
%     
%     plotdat_preproc.de_ic.R      = plotdat_preproc.de.all(ind_Rew_de_ic,:);
%     plotdat_preproc.de_ic.nR     = plotdat_preproc.de.all(ind_nRew_de_ic,:);
%     plotdat_preproc.de_ic.old    = plotdat_preproc.de.all(ind_old_de_ic,:);
%     plotdat_preproc.de_ic.new    = plotdat_preproc.de.all(ind_new_de_ic,:);
%     plotdat_preproc.de_ic.RO     = plotdat_preproc.de.all(ind_RewO_de_ic,:);
%     plotdat_preproc.de_ic.RN     = plotdat_preproc.de.all(ind_RewN_de_ic,:);
%     plotdat_preproc.de_ic.nRO    = plotdat_preproc.de.all(ind_nRewO_de_ic,:);
%     plotdat_preproc.de_ic.nRN    = plotdat_preproc.de.all(ind_nRewN_de_ic,:);
    
    
    plotavg.im.R      = mean(plotdat_preproc.im.R,1);
    plotavg.im.nR     = mean(plotdat_preproc.im.nR,1);
    plotavg.im.old    = mean(plotdat_preproc.im.old,1);
    plotavg.im.new    = mean(plotdat_preproc.im.new,1);
    
%     plotavg.de.R      = mean(plotdat_preproc.de.R,1);
%     plotavg.de.nR     = mean(plotdat_preproc.de.nR,1);
%     plotavg.de.old    = mean(plotdat_preproc.de.old,1);
%     plotavg.de.new    = mean(plotdat_preproc.de.new,1);
    
    
    fprintf('\n Done \n')
    
    save([ path_clean '/rr' num2str(session) num2str(VPs(i1)) '_preproc_dat.mat'], 'plotdat_preproc','plotavg');
    
    
    %% draw
    % legends do not work -- some errors with cvx, i think
    
    cd(path_fig) % saving path for the figures
    close all
    
    maxy = 3; miny = -3;
    onsetbar=200.*ones(1,length(time{1,1})); y=linspace(miny,maxy,numel(onsetbar));
    errorbar = @ (data) (std(data)/sqrt(length(VPs)));
    
    % fancier: shadedErrorBar
    % immediate: reward vs. non-reward
    % maxy = max([max(plotdat_preproc.im.R(:)) max(plotdat_preproc.im.nR(:)) max(plotdat_preproc.de.R(:)) max(plotdat_preproc.de.nR(:))]);
    % miny = min([min(plotdat_preproc.im.R(:)) min(plotdat_preproc.im.nR(:)) min(plotdat_preproc.de.R(:)) min(plotdat_preproc.de.nR(:))])
    figure;
%     subplot(1,2,1);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.im.R,plotdat_preproc.im.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.im.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('immediate','Fontsize',20,'Fontweight','bold')
    clear H sigbar
    % delayed seg: reward vs. non-reward
%     subplot(1,2,2);
%     plot(onsetbar,y,'g-','linewidth',2)
%     shadedErrorBar((1:size(plotdat_preproc.de.R,2)),plotdat_preproc.de.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:size(plotdat_preproc.de.nR,2)),plotdat_preproc.de.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(plotdat_preproc.de.R,plotdat_preproc.de.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 size(plotdat_preproc.de.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
%     ylabel('zscores','Fontsize',20,'Fontweight','bold');
%     title('delayed','Fontsize',20,'Fontweight','bold')
%     suptitle('Reward(red) vs. Non-Reward(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    % immediate and correct: reward vs. non-reward
    % maxy = max([max(plotdat_preproc.im.R(:)) max(plotdat_preproc.im.nR(:)) max(plotdat_preproc.de.R(:)) max(plotdat_preproc.de.nR(:))]);
    % miny = min([min(plotdat_preproc.im.R(:)) min(plotdat_preproc.im.nR(:)) min(plotdat_preproc.de.R(:)) min(plotdat_preproc.de.nR(:))])
    figure;
%     subplot(1,2,1);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im_c.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im_c.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.im_c.R,plotdat_preproc.im_c.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.im_c.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('immediate','Fontsize',20,'Fontweight','bold')
    clear H sigbar
    % delayed seg: reward vs. non-reward
%     subplot(1,2,2);
%     plot(onsetbar,y,'g-','linewidth',2)
%     shadedErrorBar((1:size(plotdat_preproc.de_c.R,2)),plotdat_preproc.de.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:size(plotdat_preproc.de_c.nR,2)),plotdat_preproc.de.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(plotdat_preproc.de_c.R,plotdat_preproc.de_c.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 size(plotdat_preproc.de_c.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
%     ylabel('zscores','Fontsize',20,'Fontweight','bold');
%     title('delayed','Fontsize',20,'Fontweight','bold')
%     suptitle('correct trials: Reward(red) vs. Non-Reward(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_CorrectTrls_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    
    % old vs new
    figure;
%     subplot(1,2,1);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:size(plotdat_preproc.im.old,2)),plotdat_preproc.im.old,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:size(plotdat_preproc.im.new,2)),plotdat_preproc.im.new,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.im.old,plotdat_preproc.im.new); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.im.old,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('immediate','Fontsize',20,'Fontweight','bold')
    clear H sigbar
    % delayed seg: reward vs. non-reward
%     subplot(1,2,2);
%     plot(onsetbar,y,'g-','linewidth',2)
%     shadedErrorBar((1:size(plotdat_preproc.de.old,2)),plotdat_preproc.de.old,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:size(plotdat_preproc.de.new,2)),plotdat_preproc.de.new,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(plotdat_preproc.de.old,plotdat_preproc.de.new); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 size(plotdat_preproc.de.old,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
%     ylabel('zscores','Fontsize',20,'Fontweight','bold');
%     title('delayed','Fontsize',20,'Fontweight','bold')
%     suptitle('old(red) vs. new(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_old(red)_vs_new(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    fprintf('\nPLOTTING DONE\n')
    
    keep session memtest VPs ...
            path_gen path_par path_dat path_behav path_clean path_fig path_stim ...
            fname_raweye fname_beh expdat i1 nulltrials
    
end

%% delayed memory test
%  from this line on, we'll analyse memory data

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
% VPs = [1204 1107 1208 1109 1210];
VPs = [1115];
session = 1;

%% plot and create masks

fname_raweye = []; fname_beh =[]; expdat = [];
for i1 = 1:length(VPs)
    
    
    fname_raweye{i1,1}  = [ 'ri' num2str(1) num2str(VPs(i1)) '.asc' ];
    fname_raweye{i1,1}  = [ 'rd' num2str(1) num2str(VPs(i1)) '.asc' ];
    fname_beh{i1,1}     = [ num2str(VPs(i1)) '_delayed.mat' ];
    
    expdat{i1,1} = load([path_behav fname_beh{i1,1}]);
    
    
    % call in cleaned data
    plotdat_im = load([ path_clean '/ri' num2str(session) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
    plotdat_de = load([ path_clean '/rd' num2str(session) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
    
    fprintf('\n Preparation done \n')
    
    %% housekeeping
    
    
    % clean the null trials and rejected trials: breaks at trial 61, 122
%     if VPs(i1) == 1204
%         plotdat_im.eye_rv.immediate.trl_raus(1) = [];
%     end
    
%     rmrow_im   = find(plotdat_im.eye_rv.immediate.trl_raus==1);
    rmrow_de   = find(plotdat_de.eye_rv.delayed.trl_raus==1);
    
%     imdat = expdat{i1,1}.dat.day1.memorytest.immediate.config.stim.stimlist_all(:,[2 4]);
    dedat = expdat{i1,1}.dat.day1.memorytest.delayed.config.stim.stimlist_all(:,[2 4]);
    
%     behdat_i= expdat{i1,1}.dat.day1.memorytest.immediate.results.all;
    behdat_d= expdat{i1,1}.dat.day1.memorytest.delayed.results.all;
    
%     imdat(rmrow_im,:) = [];
    dedat(rmrow_de,:) = [];
%     behdat_i(rmrow_im,:) = [];
    behdat_d(rmrow_de,:) = [];
    
    % indexing: session 1
    
    % all trials
    rewcond  = str2num(expdat{i1,1}.dat.RewardCategory(9)); % rewarding condition, human, was indexed as 2 in our data structure
    if rewcond==1; puncond = 2; else; puncond = 1; end
    
%     ind_Rew_im      = cell2mat(imdat(:,1) == rewcond);
%     ind_nRew_im     = cell2mat(imdat(:,1) == puncond);
    ind_Rew_de      = cell2mat(dedat(:,1) == rewcond);
    ind_nRew_de     = cell2mat(dedat(:,1) == puncond);
    
%     ind_old_im      = cell2mat(imdat(:,2) == 1);
%     ind_new_im      = cell2mat(imdat(:,2) == 2);
    ind_old_de      = cell2mat(dedat(:,2) == 1);
    ind_new_de      = cell2mat(dedat(:,2) == 2);
    
%     ind_RewO_im     = cell2mat((imdat(:,1) == rewcond) & imdat(:,2) == 1);
%     ind_RewN_im     = cell2mat((imdat(:,1) == rewcond) & imdat(:,2) == 2);
%     ind_nRewO_im    = cell2mat((imdat(:,1) == puncond) & imdat(:,2) == 1);
%     ind_nRewN_im    = cell2mat((imdat(:,1) == puncond) & imdat(:,2) == 2);
    
    ind_RewO_de     = cell2mat((dedat(:,1) == rewcond) & dedat(:,2) == 1);
    ind_RewN_de     = cell2mat((dedat(:,1) == rewcond) & dedat(:,2) == 2);
    ind_nRewO_de    = cell2mat((dedat(:,1) == puncond) & dedat(:,2) == 1);
    ind_nRewN_de    = cell2mat((dedat(:,1) == puncond) & dedat(:,2) == 2);
    
    
    % correct/incorrect sorted
    
    % the corrects!
%     ind_Rew_im_c    = ind_Rew_im == 1 & behdat_i(:,2) == 1;
%     ind_nRew_im_c   = ind_nRew_im == 1 & behdat_i(:,2) == 1;
    ind_Rew_de_c    = ind_Rew_de == 1 & behdat_d(:,2) == 1;
    ind_nRew_de_c   = ind_nRew_de == 1 & behdat_d(:,2) == 1;
      
%     ind_old_im_c    = ind_old_im == 1 & behdat_i(:,2) == 1;
%     ind_new_im_c    = ind_new_im == 1 & behdat_i(:,2) == 1;
    ind_old_de_c    = ind_old_de == 1 & behdat_d(:,2) == 1;
    ind_new_de_c    = ind_new_de == 1 & behdat_d(:,2) == 1;
    
%     ind_RewO_im_c   = ind_RewO_im ==1 & behdat_i(:,2) == 1;
%     ind_RewN_im_c   = ind_RewN_im ==1 & behdat_i(:,2) == 1;
    ind_RewO_de_c   = ind_RewO_de ==1 & behdat_d(:,2) == 1;
    ind_RewN_de_c   = ind_RewN_de ==1 & behdat_d(:,2) == 1;
    
%     ind_nRewO_im_c  = ind_nRewO_im ==1 & behdat_i(:,2) == 1;
%     ind_nRewN_im_c  = ind_nRewN_im ==1 & behdat_i(:,2) == 1;
    ind_nRewO_de_c  = ind_nRewO_de ==1 & behdat_d(:,2) == 1;
    ind_nRewN_de_c  = ind_nRewN_de ==1 & behdat_d(:,2) == 1;
    
    
    % the incorrects!
%     ind_Rew_im_ic    = ind_Rew_im == 1 & behdat_i(:,2) == 0;
%     ind_nRew_im_ic   = ind_nRew_im == 1 & behdat_i(:,2) == 0;
    ind_Rew_de_ic    = ind_Rew_de == 1 & behdat_d(:,2) == 0;
    ind_nRew_de_ic   = ind_nRew_de == 1 & behdat_d(:,2) == 0;
        
%     ind_old_im_ic    = ind_old_im == 1 & behdat_i(:,2) == 0;
%     ind_new_im_ic    = ind_new_im == 1 & behdat_i(:,2) == 0;
    ind_old_de_ic    = ind_old_de == 1 & behdat_d(:,2) == 0;
    ind_new_de_ic    = ind_new_de == 1 & behdat_d(:,2) == 0;
    
%     ind_RewO_im_ic   = ind_RewO_im ==1 & behdat_i(:,2) == 0;
%     ind_RewN_im_ic   = ind_RewN_im ==1 & behdat_i(:,2) == 0;
    ind_RewO_de_ic   = ind_RewO_de ==1 & behdat_d(:,2) == 0;
    ind_RewN_de_ic   = ind_RewN_de ==1 & behdat_d(:,2) == 0;
    
%     ind_nRewO_im_ic  = ind_nRewO_im ==1 & behdat_i(:,2) == 0;
%     ind_nRewN_im_ic  = ind_nRewN_im ==1 & behdat_i(:,2) == 0;
    ind_nRewO_de_ic  = ind_nRewO_de ==1 & behdat_d(:,2) == 0;
    ind_nRewN_de_ic  = ind_nRewN_de ==1 & behdat_d(:,2) == 0;
    
    % assign data
%     im_eye          = plotdat_im.eye_rv.trial';
%     im_Rew          = im_eye(ind_Rew_im,1);
%     im_nRew         = im_eye(ind_nRew_im,1);
%     im_old          = im_eye(ind_old_im,1);
%     im_new          = im_eye(ind_new_im,1);
    
    de_eye          = plotdat_de.eye_rv.trial';
    de_Rew          = de_eye(ind_Rew_de,1);
    de_nRew         = de_eye(ind_nRew_de,1);
    de_old          = de_eye(ind_old_de,1);
    de_new          = de_eye(ind_new_de,1);
    
    
    % gather and standardise
    time = plotdat_im.eye_rv.time;
    timewindow_bsl = [1 200]; % -200ms ~ 0ms
    plotdat        = [];
    
%     plotdat.im.all = []; % immediate
%     for i = 1:size(im_eye,1) % for the number of trials
%         plotdat.im.all(i,:) = im_eye{i,:}; % create one big matrix
%     end
    plotdat.de.all = []; % delayed
    for i = 1:size(de_eye,1)
        plotdat.de.all(i,:) = de_eye{i,:};
    end
    
%     % immediate
%     plotdat_z = []; clear i
%     for i = 1:size(plotdat.im.all,1)
%         meandat = nanmean(plotdat.im.all,2); stddat = nanstd(plotdat.im.all');
%         plotdat_z.im.all = plotdat.im.all - repmat(meandat,1,size(plotdat.im.all,2));
%         plotdat_z.im.all = plotdat.im.all ./ repmat(stddat',1,size(plotdat.im.all,2));
%         
%         bsldat = nanmean(plotdat_z.im.all(:,timewindow_bsl(1):timewindow_bsl(2))');
%         plotdat_preproc.im.all = plotdat_z.im.all-repmat(bsldat',1,size(plotdat_z.im.all,2));
%     end
%     plotdat_preproc.im.R      = plotdat_preproc.im.all(ind_Rew_im,:);
%     plotdat_preproc.im.nR     = plotdat_preproc.im.all(ind_nRew_im,:);
%     plotdat_preproc.im.old    = plotdat_preproc.im.all(ind_old_im,:);
%     plotdat_preproc.im.new    = plotdat_preproc.im.all(ind_new_im,:);
%     plotdat_preproc.im.RO     = plotdat_preproc.im.all(ind_RewO_im,:);
%     plotdat_preproc.im.RN     = plotdat_preproc.im.all(ind_RewN_im,:);
%     plotdat_preproc.im.nRO    = plotdat_preproc.im.all(ind_nRewO_im,:);
%     plotdat_preproc.im.nRN    = plotdat_preproc.im.all(ind_nRewN_im,:);
%     
%     plotdat_preproc.im_c.R      = plotdat_preproc.im.all(ind_Rew_im_c,:);
%     plotdat_preproc.im_c.nR     = plotdat_preproc.im.all(ind_nRew_im_c,:);
%     plotdat_preproc.im_c.old    = plotdat_preproc.im.all(ind_old_im_c,:);
%     plotdat_preproc.im_c.new    = plotdat_preproc.im.all(ind_new_im_c,:);
%     plotdat_preproc.im_c.RO     = plotdat_preproc.im.all(ind_RewO_im_c,:);
%     plotdat_preproc.im_c.RN     = plotdat_preproc.im.all(ind_RewN_im_c,:);
%     plotdat_preproc.im_c.nRO    = plotdat_preproc.im.all(ind_nRewO_im_c,:);
%     plotdat_preproc.im_c.nRN    = plotdat_preproc.im.all(ind_nRewN_im_c,:);
%     
%     plotdat_preproc.im_ic.R      = plotdat_preproc.im.all(ind_Rew_im_ic,:);
%     plotdat_preproc.im_ic.nR     = plotdat_preproc.im.all(ind_nRew_im_ic,:);
%     plotdat_preproc.im_ic.old    = plotdat_preproc.im.all(ind_old_im_ic,:);
%     plotdat_preproc.im_ic.new    = plotdat_preproc.im.all(ind_new_im_ic,:);
%     plotdat_preproc.im_ic.RO     = plotdat_preproc.im.all(ind_RewO_im_ic,:);
%     plotdat_preproc.im_ic.RN     = plotdat_preproc.im.all(ind_RewN_im_ic,:);
%     plotdat_preproc.im_ic.nRO    = plotdat_preproc.im.all(ind_nRewO_im_ic,:);
%     plotdat_preproc.im_ic.nRN    = plotdat_preproc.im.all(ind_nRewN_im_ic,:);
    
    % delayed
    plotdat_z = []; clear i
    for i = 1:size(plotdat.de.all,1)
        meandat = nanmean(plotdat.de.all,2); stddat = nanstd(plotdat.de.all');
        plotdat_z.de.all = plotdat.de.all - repmat(meandat,1,size(plotdat.de.all,2));
        plotdat_z.de.all = plotdat.de.all ./ repmat(stddat',1,size(plotdat.de.all,2));
        
        bsldat = nanmean(plotdat_z.de.all(:,timewindow_bsl(1):timewindow_bsl(2))');
        plotdat_preproc.de.all = plotdat_z.de.all-repmat(bsldat',1,size(plotdat_z.de.all,2));
    end
    plotdat_preproc.de.R      = plotdat_preproc.de.all(ind_Rew_de,:);
    plotdat_preproc.de.nR     = plotdat_preproc.de.all(ind_nRew_de,:);
    plotdat_preproc.de.old    = plotdat_preproc.de.all(ind_old_de,:);
    plotdat_preproc.de.new    = plotdat_preproc.de.all(ind_new_de,:);
    plotdat_preproc.de.RO     = plotdat_preproc.de.all(ind_RewO_de,:);
    plotdat_preproc.de.RN     = plotdat_preproc.de.all(ind_RewN_de,:);
    plotdat_preproc.de.nRO    = plotdat_preproc.de.all(ind_nRewO_de,:);
    plotdat_preproc.de.nRN    = plotdat_preproc.de.all(ind_nRewN_de,:);
    
    plotdat_preproc.de_c.R      = plotdat_preproc.de.all(ind_Rew_de_c,:);
    plotdat_preproc.de_c.nR     = plotdat_preproc.de.all(ind_nRew_de_c,:);
    plotdat_preproc.de_c.old    = plotdat_preproc.de.all(ind_old_de_c,:);
    plotdat_preproc.de_c.new    = plotdat_preproc.de.all(ind_new_de_c,:);
    plotdat_preproc.de_c.RO     = plotdat_preproc.de.all(ind_RewO_de_c,:);
    plotdat_preproc.de_c.RN     = plotdat_preproc.de.all(ind_RewN_de_c,:);
    plotdat_preproc.de_c.nRO    = plotdat_preproc.de.all(ind_nRewO_de_c,:);
    plotdat_preproc.de_c.nRN    = plotdat_preproc.de.all(ind_nRewN_de_c,:);
    
    plotdat_preproc.de_ic.R      = plotdat_preproc.de.all(ind_Rew_de_ic,:);
    plotdat_preproc.de_ic.nR     = plotdat_preproc.de.all(ind_nRew_de_ic,:);
    plotdat_preproc.de_ic.old    = plotdat_preproc.de.all(ind_old_de_ic,:);
    plotdat_preproc.de_ic.new    = plotdat_preproc.de.all(ind_new_de_ic,:);
    plotdat_preproc.de_ic.RO     = plotdat_preproc.de.all(ind_RewO_de_ic,:);
    plotdat_preproc.de_ic.RN     = plotdat_preproc.de.all(ind_RewN_de_ic,:);
    plotdat_preproc.de_ic.nRO    = plotdat_preproc.de.all(ind_nRewO_de_ic,:);
    plotdat_preproc.de_ic.nRN    = plotdat_preproc.de.all(ind_nRewN_de_ic,:);
    
    
%     plotavg.im.R      = mean(plotdat_preproc.im.R,1);
%     plotavg.im.nR     = mean(plotdat_preproc.im.nR,1);
%     plotavg.im.old    = mean(plotdat_preproc.im.old,1);
%     plotavg.im.new    = mean(plotdat_preproc.im.new,1);
    
    plotavg.de.R      = mean(plotdat_preproc.de.R,1);
    plotavg.de.nR     = mean(plotdat_preproc.de.nR,1);
    plotavg.de.old    = mean(plotdat_preproc.de.old,1);
    plotavg.de.new    = mean(plotdat_preproc.de.new,1);
    
    
    fprintf('\n Done \n')
    
    save([ path_clean '/rr' num2str(session) num2str(VPs(i1)) '_preproc_dat.mat'], 'plotdat_preproc','plotavg');
    
    
    %% draw
    % legends do not work -- some errors with cvx, i think
    
    cd(path_fig) % saving path for the figures
    close all
    
    maxy = 3; miny = -3;
    onsetbar=200.*ones(1,length(time{1,1})); y=linspace(miny,maxy,numel(onsetbar));
    errorbar = @ (data) (std(data)/sqrt(length(VPs)));
    
    % fancier: shadedErrorBar
    % immediate: reward vs. non-reward
    % maxy = max([max(plotdat_preproc.im.R(:)) max(plotdat_preproc.im.nR(:)) max(plotdat_preproc.de.R(:)) max(plotdat_preproc.de.nR(:))]);
    % miny = min([min(plotdat_preproc.im.R(:)) min(plotdat_preproc.im.nR(:)) min(plotdat_preproc.de.R(:)) min(plotdat_preproc.de.nR(:))])
    figure;
%     subplot(1,2,1);
%     plot(onsetbar,y,'g-','linewidth',2)
%     shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(plotdat_preproc.im.R,plotdat_preproc.im.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 size(plotdat_preproc.im.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
%     ylabel('zscores','Fontsize',20,'Fontweight','bold');
%     title('immediate','Fontsize',20,'Fontweight','bold')
%     clear H sigbar
    % delayed seg: reward vs. non-reward
    subplot(1,2,2);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:size(plotdat_preproc.de.R,2)),plotdat_preproc.de.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:size(plotdat_preproc.de.nR,2)),plotdat_preproc.de.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.de.R,plotdat_preproc.de.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.de.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('delayed','Fontsize',20,'Fontweight','bold')
    suptitle('Reward(red) vs. Non-Reward(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    % immediate and correct: reward vs. non-reward
    % maxy = max([max(plotdat_preproc.im.R(:)) max(plotdat_preproc.im.nR(:)) max(plotdat_preproc.de.R(:)) max(plotdat_preproc.de.nR(:))]);
    % miny = min([min(plotdat_preproc.im.R(:)) min(plotdat_preproc.im.nR(:)) min(plotdat_preproc.de.R(:)) min(plotdat_preproc.de.nR(:))])
%     figure;
%     subplot(1,2,1);
%     plot(onsetbar,y,'g-','linewidth',2)
%     shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im_c.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im_c.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(plotdat_preproc.im_c.R,plotdat_preproc.im_c.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 size(plotdat_preproc.im_c.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
%     ylabel('zscores','Fontsize',20,'Fontweight','bold');
%     title('immediate','Fontsize',20,'Fontweight','bold')
%     clear H sigbar
    % delayed seg: reward vs. non-reward
    subplot(1,2,2);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:size(plotdat_preproc.de_c.R,2)),plotdat_preproc.de.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:size(plotdat_preproc.de_c.nR,2)),plotdat_preproc.de.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.de_c.R,plotdat_preproc.de_c.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.de_c.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('delayed','Fontsize',20,'Fontweight','bold')
    suptitle('correct trials: Reward(red) vs. Non-Reward(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_CorrectTrls_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    
    % old vs new
    figure;
%     subplot(1,2,1);
%     plot(onsetbar,y,'g-','linewidth',2)
%     shadedErrorBar((1:size(plotdat_preproc.im.old,2)),plotdat_preproc.im.old,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%     shadedErrorBar((1:size(plotdat_preproc.im.new,2)),plotdat_preproc.im.new,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%     H = ttest2(plotdat_preproc.im.old,plotdat_preproc.im.new); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
%     xlim([1 size(plotdat_preproc.im.old,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
%     xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
%     ylabel('zscores','Fontsize',20,'Fontweight','bold');
%     title('immediate','Fontsize',20,'Fontweight','bold')
%     clear H sigbar
    % delayed seg: reward vs. non-reward
    subplot(1,2,2);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:size(plotdat_preproc.de.old,2)),plotdat_preproc.de.old,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:size(plotdat_preproc.de.new,2)),plotdat_preproc.de.new,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.de.old,plotdat_preproc.de.new); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.de.old,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('delayed','Fontsize',20,'Fontweight','bold')
    suptitle('old(red) vs. new(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_old(red)_vs_new(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    fprintf('\nPLOTTING DONE\n')
    
    keep session memtest VPs ...
            path_gen path_par path_dat path_behav path_clean path_fig path_stim ...
            fname_raweye fname_beh expdat i1
    
end


%% immediate & delayed memory test
%  from this line on, we'll analyse memory data

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
% VPs = [1204 1107 1208 1109 1210];
VPs = [1109 1210 1111 1212 1115];
session = 1;

%% plot and create masks

fname_raweye = []; fname_beh =[]; expdat = [];
for i1 = 1:length(VPs)
    
    
    fname_raweye{i1,1}  = [ 'ri' num2str(1) num2str(VPs(i1)) '.asc' ];
    fname_raweye{i1,1}  = [ 'rd' num2str(1) num2str(VPs(i1)) '.asc' ];
    fname_beh{i1,1}     = [ num2str(VPs(i1)) '_delayed.mat' ];
    
    expdat{i1,1} = load([path_behav fname_beh{i1,1}]);
    
    
    % call in cleaned data
    plotdat_im = load([ path_clean '/ri' num2str(session) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
    plotdat_de = load([ path_clean '/rd' num2str(session) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
    
    fprintf('\n Preparation done \n')
    
    %% housekeeping
    
    
    % clean the null trials and rejected trials: breaks at trial 61, 122
    if VPs(i1) == 1204
        plotdat_im.eye_rv.immediate.trl_raus(1) = [];
    end
    
    rmrow_im   = find(plotdat_im.eye_rv.immediate.trl_raus==1);
    rmrow_de   = find(plotdat_de.eye_rv.delayed.trl_raus==1);
    
    imdat = expdat{i1,1}.dat.day1.memorytest.immediate.config.stim.stimlist_all(:,[2 4]);
    dedat = expdat{i1,1}.dat.day1.memorytest.delayed.config.stim.stimlist_all(:,[2 4]);
    
    behdat_i= expdat{i1,1}.dat.day1.memorytest.immediate.results.all;
    behdat_d= expdat{i1,1}.dat.day1.memorytest.delayed.results.all;
    
    imdat(rmrow_im,:) = [];
    dedat(rmrow_de,:) = [];
    behdat_i(rmrow_im,:) = [];
    behdat_d(rmrow_de,:) = [];
    
    % indexing: session 1
    
    % all trials
    rewcond  = str2num(expdat{i1,1}.dat.RewardCategory(9)); % rewarding condition, human, was indexed as 2 in our data structure
    if rewcond==1; puncond = 2; else; puncond = 1; end
    
    ind_Rew_im      = cell2mat(imdat(:,1) == rewcond);
    ind_nRew_im     = cell2mat(imdat(:,1) == puncond);
    ind_Rew_de      = cell2mat(dedat(:,1) == rewcond);
    ind_nRew_de     = cell2mat(dedat(:,1) == puncond);
    
    ind_old_im      = cell2mat(imdat(:,2) == 1);
    ind_new_im      = cell2mat(imdat(:,2) == 2);
    ind_old_de      = cell2mat(dedat(:,2) == 1);
    ind_new_de      = cell2mat(dedat(:,2) == 2);
    
    ind_RewO_im     = cell2mat((imdat(:,1) == rewcond) & imdat(:,2) == 1);
    ind_RewN_im     = cell2mat((imdat(:,1) == rewcond) & imdat(:,2) == 2);
    ind_nRewO_im    = cell2mat((imdat(:,1) == puncond) & imdat(:,2) == 1);
    ind_nRewN_im    = cell2mat((imdat(:,1) == puncond) & imdat(:,2) == 2);
    
    ind_RewO_de     = cell2mat((dedat(:,1) == rewcond) & dedat(:,2) == 1);
    ind_RewN_de     = cell2mat((dedat(:,1) == rewcond) & dedat(:,2) == 2);
    ind_nRewO_de    = cell2mat((dedat(:,1) == puncond) & dedat(:,2) == 1);
    ind_nRewN_de    = cell2mat((dedat(:,1) == puncond) & dedat(:,2) == 2);
    
    
    % correct/incorrect sorted
    
    % the corrects!
    ind_Rew_im_c    = ind_Rew_im == 1 & behdat_i(:,2) == 1;
    ind_nRew_im_c   = ind_nRew_im == 1 & behdat_i(:,2) == 1;
    ind_Rew_de_c    = ind_Rew_de == 1 & behdat_d(:,2) == 1;
    ind_nRew_de_c   = ind_nRew_de == 1 & behdat_d(:,2) == 1;
      
    ind_old_im_c    = ind_old_im == 1 & behdat_i(:,2) == 1;
    ind_new_im_c    = ind_new_im == 1 & behdat_i(:,2) == 1;
    ind_old_de_c    = ind_old_de == 1 & behdat_d(:,2) == 1;
    ind_new_de_c    = ind_new_de == 1 & behdat_d(:,2) == 1;
    
    ind_RewO_im_c   = ind_RewO_im ==1 & behdat_i(:,2) == 1;
    ind_RewN_im_c   = ind_RewN_im ==1 & behdat_i(:,2) == 1;
    ind_RewO_de_c   = ind_RewO_de ==1 & behdat_d(:,2) == 1;
    ind_RewN_de_c   = ind_RewN_de ==1 & behdat_d(:,2) == 1;
    
    ind_nRewO_im_c  = ind_nRewO_im ==1 & behdat_i(:,2) == 1;
    ind_nRewN_im_c  = ind_nRewN_im ==1 & behdat_i(:,2) == 1;
    ind_nRewO_de_c  = ind_nRewO_de ==1 & behdat_d(:,2) == 1;
    ind_nRewN_de_c  = ind_nRewN_de ==1 & behdat_d(:,2) == 1;
    
    
    % the incorrects!
    ind_Rew_im_ic    = ind_Rew_im == 1 & behdat_i(:,2) == 0;
    ind_nRew_im_ic   = ind_nRew_im == 1 & behdat_i(:,2) == 0;
    ind_Rew_de_ic    = ind_Rew_de == 1 & behdat_d(:,2) == 0;
    ind_nRew_de_ic   = ind_nRew_de == 1 & behdat_d(:,2) == 0;
        
    ind_old_im_ic    = ind_old_im == 1 & behdat_i(:,2) == 0;
    ind_new_im_ic    = ind_new_im == 1 & behdat_i(:,2) == 0;
    ind_old_de_ic    = ind_old_de == 1 & behdat_d(:,2) == 0;
    ind_new_de_ic    = ind_new_de == 1 & behdat_d(:,2) == 0;
    
    ind_RewO_im_ic   = ind_RewO_im ==1 & behdat_i(:,2) == 0;
    ind_RewN_im_ic   = ind_RewN_im ==1 & behdat_i(:,2) == 0;
    ind_RewO_de_ic   = ind_RewO_de ==1 & behdat_d(:,2) == 0;
    ind_RewN_de_ic   = ind_RewN_de ==1 & behdat_d(:,2) == 0;
    
    ind_nRewO_im_ic  = ind_nRewO_im ==1 & behdat_i(:,2) == 0;
    ind_nRewN_im_ic  = ind_nRewN_im ==1 & behdat_i(:,2) == 0;
    ind_nRewO_de_ic  = ind_nRewO_de ==1 & behdat_d(:,2) == 0;
    ind_nRewN_de_ic  = ind_nRewN_de ==1 & behdat_d(:,2) == 0;
    
    % assign data
    im_eye          = plotdat_im.eye_rv.trial';
    im_Rew          = im_eye(ind_Rew_im,1);
    im_nRew         = im_eye(ind_nRew_im,1);
    im_old          = im_eye(ind_old_im,1);
    im_new          = im_eye(ind_new_im,1);
    
    de_eye          = plotdat_de.eye_rv.trial';
    de_Rew          = de_eye(ind_Rew_de,1);
    de_nRew         = de_eye(ind_nRew_de,1);
    de_old          = de_eye(ind_old_de,1);
    de_new          = de_eye(ind_new_de,1);
    
    
    % gather and standardise
    time = plotdat_im.eye_rv.time;
    timewindow_bsl = [1 200]; % -200ms ~ 0ms
    plotdat        = [];
    
    plotdat.im.all = []; % immediate
    for i = 1:size(im_eye,1) % for the number of trials
        plotdat.im.all(i,:) = im_eye{i,:}; % create one big matrix
    end
    plotdat.de.all = []; % delayed
    for i = 1:size(de_eye,1)
        plotdat.de.all(i,:) = de_eye{i,:};
    end
    
    % immediate
    plotdat_z = []; clear i
    for i = 1:size(plotdat.im.all,1)
        meandat = nanmean(plotdat.im.all,2); stddat = nanstd(plotdat.im.all');
        plotdat_z.im.all = plotdat.im.all - repmat(meandat,1,size(plotdat.im.all,2));
        plotdat_z.im.all = plotdat.im.all ./ repmat(stddat',1,size(plotdat.im.all,2));
        
        bsldat = nanmean(plotdat_z.im.all(:,timewindow_bsl(1):timewindow_bsl(2))');
        plotdat_preproc.im.all = plotdat_z.im.all-repmat(bsldat',1,size(plotdat_z.im.all,2));
    end
    plotdat_preproc.im.R      = plotdat_preproc.im.all(ind_Rew_im,:);
    plotdat_preproc.im.nR     = plotdat_preproc.im.all(ind_nRew_im,:);
    plotdat_preproc.im.old    = plotdat_preproc.im.all(ind_old_im,:);
    plotdat_preproc.im.new    = plotdat_preproc.im.all(ind_new_im,:);
    plotdat_preproc.im.RO     = plotdat_preproc.im.all(ind_RewO_im,:);
    plotdat_preproc.im.RN     = plotdat_preproc.im.all(ind_RewN_im,:);
    plotdat_preproc.im.nRO    = plotdat_preproc.im.all(ind_nRewO_im,:);
    plotdat_preproc.im.nRN    = plotdat_preproc.im.all(ind_nRewN_im,:);
    
    plotdat_preproc.im_c.R      = plotdat_preproc.im.all(ind_Rew_im_c,:);
    plotdat_preproc.im_c.nR     = plotdat_preproc.im.all(ind_nRew_im_c,:);
    plotdat_preproc.im_c.old    = plotdat_preproc.im.all(ind_old_im_c,:);
    plotdat_preproc.im_c.new    = plotdat_preproc.im.all(ind_new_im_c,:);
    plotdat_preproc.im_c.RO     = plotdat_preproc.im.all(ind_RewO_im_c,:);
    plotdat_preproc.im_c.RN     = plotdat_preproc.im.all(ind_RewN_im_c,:);
    plotdat_preproc.im_c.nRO    = plotdat_preproc.im.all(ind_nRewO_im_c,:);
    plotdat_preproc.im_c.nRN    = plotdat_preproc.im.all(ind_nRewN_im_c,:);
    
    plotdat_preproc.im_ic.R      = plotdat_preproc.im.all(ind_Rew_im_ic,:);
    plotdat_preproc.im_ic.nR     = plotdat_preproc.im.all(ind_nRew_im_ic,:);
    plotdat_preproc.im_ic.old    = plotdat_preproc.im.all(ind_old_im_ic,:);
    plotdat_preproc.im_ic.new    = plotdat_preproc.im.all(ind_new_im_ic,:);
    plotdat_preproc.im_ic.RO     = plotdat_preproc.im.all(ind_RewO_im_ic,:);
    plotdat_preproc.im_ic.RN     = plotdat_preproc.im.all(ind_RewN_im_ic,:);
    plotdat_preproc.im_ic.nRO    = plotdat_preproc.im.all(ind_nRewO_im_ic,:);
    plotdat_preproc.im_ic.nRN    = plotdat_preproc.im.all(ind_nRewN_im_ic,:);
    
    % delayed
    plotdat_z = []; clear i
    for i = 1:size(plotdat.de.all,1)
        meandat = nanmean(plotdat.de.all,2); stddat = nanstd(plotdat.de.all');
        plotdat_z.de.all = plotdat.de.all - repmat(meandat,1,size(plotdat.de.all,2));
        plotdat_z.de.all = plotdat.de.all ./ repmat(stddat',1,size(plotdat.de.all,2));
        
        bsldat = nanmean(plotdat_z.de.all(:,timewindow_bsl(1):timewindow_bsl(2))');
        plotdat_preproc.de.all = plotdat_z.de.all-repmat(bsldat',1,size(plotdat_z.de.all,2));
    end
    plotdat_preproc.de.R      = plotdat_preproc.de.all(ind_Rew_de,:);
    plotdat_preproc.de.nR     = plotdat_preproc.de.all(ind_nRew_de,:);
    plotdat_preproc.de.old    = plotdat_preproc.de.all(ind_old_de,:);
    plotdat_preproc.de.new    = plotdat_preproc.de.all(ind_new_de,:);
    plotdat_preproc.de.RO     = plotdat_preproc.de.all(ind_RewO_de,:);
    plotdat_preproc.de.RN     = plotdat_preproc.de.all(ind_RewN_de,:);
    plotdat_preproc.de.nRO    = plotdat_preproc.de.all(ind_nRewO_de,:);
    plotdat_preproc.de.nRN    = plotdat_preproc.de.all(ind_nRewN_de,:);
    
    plotdat_preproc.de_c.R      = plotdat_preproc.de.all(ind_Rew_de_c,:);
    plotdat_preproc.de_c.nR     = plotdat_preproc.de.all(ind_nRew_de_c,:);
    plotdat_preproc.de_c.old    = plotdat_preproc.de.all(ind_old_de_c,:);
    plotdat_preproc.de_c.new    = plotdat_preproc.de.all(ind_new_de_c,:);
    plotdat_preproc.de_c.RO     = plotdat_preproc.de.all(ind_RewO_de_c,:);
    plotdat_preproc.de_c.RN     = plotdat_preproc.de.all(ind_RewN_de_c,:);
    plotdat_preproc.de_c.nRO    = plotdat_preproc.de.all(ind_nRewO_de_c,:);
    plotdat_preproc.de_c.nRN    = plotdat_preproc.de.all(ind_nRewN_de_c,:);
    
    plotdat_preproc.de_ic.R      = plotdat_preproc.de.all(ind_Rew_de_ic,:);
    plotdat_preproc.de_ic.nR     = plotdat_preproc.de.all(ind_nRew_de_ic,:);
    plotdat_preproc.de_ic.old    = plotdat_preproc.de.all(ind_old_de_ic,:);
    plotdat_preproc.de_ic.new    = plotdat_preproc.de.all(ind_new_de_ic,:);
    plotdat_preproc.de_ic.RO     = plotdat_preproc.de.all(ind_RewO_de_ic,:);
    plotdat_preproc.de_ic.RN     = plotdat_preproc.de.all(ind_RewN_de_ic,:);
    plotdat_preproc.de_ic.nRO    = plotdat_preproc.de.all(ind_nRewO_de_ic,:);
    plotdat_preproc.de_ic.nRN    = plotdat_preproc.de.all(ind_nRewN_de_ic,:);
    
    
    plotavg.im.R      = mean(plotdat_preproc.im.R,1);
    plotavg.im.nR     = mean(plotdat_preproc.im.nR,1);
    plotavg.im.old    = mean(plotdat_preproc.im.old,1);
    plotavg.im.new    = mean(plotdat_preproc.im.new,1);
    
    plotavg.de.R      = mean(plotdat_preproc.de.R,1);
    plotavg.de.nR     = mean(plotdat_preproc.de.nR,1);
    plotavg.de.old    = mean(plotdat_preproc.de.old,1);
    plotavg.de.new    = mean(plotdat_preproc.de.new,1);
    
    
    fprintf('\n Done \n')
    
    save([ path_clean '/rr' num2str(session) num2str(VPs(i1)) '_preproc_dat.mat'], 'plotdat_preproc','plotavg');
    
    
    %% draw
    % legends do not work -- some errors with cvx, i think
    
    cd(path_fig) % saving path for the figures
    close all
    
    maxy = 3; miny = -3;
    onsetbar=200.*ones(1,length(time{1,1})); y=linspace(miny,maxy,numel(onsetbar));
    errorbar = @ (data) (std(data)/sqrt(length(VPs)));
    
    % fancier: shadedErrorBar
    % immediate: reward vs. non-reward
    % maxy = max([max(plotdat_preproc.im.R(:)) max(plotdat_preproc.im.nR(:)) max(plotdat_preproc.de.R(:)) max(plotdat_preproc.de.nR(:))]);
    % miny = min([min(plotdat_preproc.im.R(:)) min(plotdat_preproc.im.nR(:)) min(plotdat_preproc.de.R(:)) min(plotdat_preproc.de.nR(:))])
    figure;
    subplot(1,2,1);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.im.R,plotdat_preproc.im.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.im.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('immediate','Fontsize',20,'Fontweight','bold')
    clear H sigbar
    % delayed seg: reward vs. non-reward
    subplot(1,2,2);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:size(plotdat_preproc.de.R,2)),plotdat_preproc.de.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:size(plotdat_preproc.de.nR,2)),plotdat_preproc.de.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.de.R,plotdat_preproc.de.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.de.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('delayed','Fontsize',20,'Fontweight','bold')
    suptitle('Reward(red) vs. Non-Reward(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    % immediate and correct: reward vs. non-reward
    % maxy = max([max(plotdat_preproc.im.R(:)) max(plotdat_preproc.im.nR(:)) max(plotdat_preproc.de.R(:)) max(plotdat_preproc.de.nR(:))]);
    % miny = min([min(plotdat_preproc.im.R(:)) min(plotdat_preproc.im.nR(:)) min(plotdat_preproc.de.R(:)) min(plotdat_preproc.de.nR(:))])
    figure;
    subplot(1,2,1);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im_c.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:length(time{1,1})),plotdat_preproc.im_c.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.im_c.R,plotdat_preproc.im_c.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.im_c.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('immediate','Fontsize',20,'Fontweight','bold')
    clear H sigbar
    % delayed seg: reward vs. non-reward
    subplot(1,2,2);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:size(plotdat_preproc.de_c.R,2)),plotdat_preproc.de.R,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:size(plotdat_preproc.de_c.nR,2)),plotdat_preproc.de.nR,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.de_c.R,plotdat_preproc.de_c.nR); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.de_c.R,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('delayed','Fontsize',20,'Fontweight','bold')
    suptitle('correct trials: Reward(red) vs. Non-Reward(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_CorrectTrls_Reward(red)_vs_Non-Reward(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    
    % old vs new
    figure;
    subplot(1,2,1);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:size(plotdat_preproc.im.old,2)),plotdat_preproc.im.old,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:size(plotdat_preproc.im.new,2)),plotdat_preproc.im.new,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.im.old,plotdat_preproc.im.new); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.im.old,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('immediate','Fontsize',20,'Fontweight','bold')
    clear H sigbar
    % delayed seg: reward vs. non-reward
    subplot(1,2,2);
    plot(onsetbar,y,'g-','linewidth',2)
    shadedErrorBar((1:size(plotdat_preproc.de.old,2)),plotdat_preproc.de.old,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
    shadedErrorBar((1:size(plotdat_preproc.de.new,2)),plotdat_preproc.de.new,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
    H = ttest2(plotdat_preproc.de.old,plotdat_preproc.de.new); sigbar = H*-0.9; sigbar(find(sigbar==0))=NaN; plot(sigbar,'k-o','LineWidth',3);
    xlim([1 size(plotdat_preproc.de.old,2)]); ylim([miny maxy]); xticklabels([0 500 1000 1500 2000 2500]);
    xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
    ylabel('zscores','Fontsize',20,'Fontweight','bold');
    title('delayed','Fontsize',20,'Fontweight','bold')
    suptitle('old(red) vs. new(blue)')
    % save
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print([num2str(VPs(i1)) '_Session' num2str(session) '_memtest_old(red)_vs_new(blue)'],'-dpdf','-r0')
    clear fig H sigbar
    close all
    
    fprintf('\nPLOTTING DONE\n')
    
    keep session memtest VPs ...
            path_gen path_par path_dat path_behav path_clean path_fig path_stim ...
            fname_raweye fname_beh expdat i1
    
end
%% rubbish bin

