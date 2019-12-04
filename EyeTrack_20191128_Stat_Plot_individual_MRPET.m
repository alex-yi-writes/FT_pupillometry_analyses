%% Pupillometric data statistics and plotting for MRPET : individuals

%% work log

%   19_06_2019  created the script

%% index (search keywords)
%  maintask
%  memory test

%% main task

%% preparation

clear;clc;close all;
% warning('off')
% addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20180613'))
path_gen    = '/Users/yeojin/Desktop/E_data/';
path_par    = [ path_gen 'EA_raw/' ];
path_behav  = [ path_par 'EAC_behav/MRPET/' ];
path_clean  = [ path_gen 'EB_cleaned/EBB_pupil/MRPET/' ];
path_fig    = [  '/Users/yeojin/Desktop/C_writings/CB_figures/MRPET/MainTask/zscores/' ];
% path_stim   = [ path_par 'EA_raw/EAA_pupil/stiminfo_all_20181004.mat' ];

% create file list and VP ID list
VPs = [4001 4002 4003 4004];
tag = [0 2; 0 2; 1 0; 1 0];
d1m = [0 2; 0 2; 0 2; 0 2]; % 1=immediate 2=delayed
d2m = [0 2; 0 2; 0 2; 0 2];


%% plot and create masks

fname_raweye = []; fname_beh =[]; expdat = [];
for i1 = 1:length(VPs)
    for d = 1:2
        
        if tag(i1,d) == 0
        else
            
            fname_raweye{i1,1}  = [ 're' num2str(d) num2str(VPs(i1)) '.asc' ];
            fname_beh{i1,1}     = [ num2str(VPs(i1)) '_' num2str(d) '.mat' ];
            
            expdat{i1,d} = load([path_behav fname_beh{i1,1}]);
            
            % call in cleaned data
            plotdat_stim  = load([path_clean 're' num2str(d) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
            plotdat_fb    = load([path_clean 're' num2str(d) num2str(VPs(i1)) '_eye_rv200_fb.mat']);
            
            fprintf('\n Preparation done \n')
            
            %% housekeeping
            
            
            % stim segment!
            % clean the null trials and rejected trials: breaks
            rmrow_nulls  = eval(['expdat{i1,d}.dat.day' num2str(d) '.maintask.config.stim.stimvec(1,find(expdat{i1,d}.dat.day'...
                num2str(d) '.maintask.config.stim.stimvec(2,:)==0))']);
            rmrow_stim   = find(plotdat_stim.eye_rv.trl_raus==1); % rejected trials from preprocessing
            emptymeans   = cellfun(@nanmean,plotdat_stim.eye_rv.trial);
            emptytrials  = find(emptymeans==0);
            
            % stimulus informations: filename, reward/no-reward, memorability
            eval(['stimdat = expdat{i1,d}.dat.day' num2str(d) '.maintask.results.trl;']);
            stimdat(isempty(stimdat)) = {NaN}; % this line, because if you convert data directly from cells to double, it automatically removes empty cells, which in our case are breaks.
            
            % remove breaks and rejected trials
            stimdat(rmrow_nulls,:) = [];
            stimdat(rmrow_stim,:)  = []; % and now, we remove the rows of trials we rejected during data cleaning
            stimdat(emptytrials,:) = [];
            
            % fb segment!
            % clean the null trials and rejected trials: breaks
            rmrow_fb   = find(plotdat_fb.eye_rv.trl_raus==1); % rejected trials from preprocessing
            emptymeans   = cellfun(@nanmean,plotdat_fb.eye_rv.trial);
            emptytrialf  = find(emptymeans==0);
            
            % stimulus informations: filename, reward/no-reward, memorability
            eval(['fbdat = expdat{i1,d}.dat.day' num2str(d) '.maintask.results.trl;']);
            fbdat(isempty(fbdat)) = {NaN}; % this line, because if you convert data directly from cells to double, it automatically removes empty cells, which in our case are breaks.
            
            % remove breaks and rejected trials
            fbdat(rmrow_nulls,:) = [];
            fbdat(rmrow_fb,:)    = []; % and now, we remove the rows of trials we rejected during data cleaning
            fbdat(emptytrialf,:) = [];
            
            % indexing
            eval(['rewcond = str2num(expdat{i1,d}.dat.day' num2str(d) '.RewardCategory(9));']); % rewarding condition, human, was indexed as 2 in our data structure
            if rewcond==1
                neucond = 2;
            elseif rewcond == 2
                neucond = 1;
            elseif rewcond == 3
                neucond = 4;
            elseif rewcond == 4
                neucond = 3;
            end
            
            ind_rew_stim    = cell2mat(stimdat(:,2)) == rewcond; % marking the trials that are rewards, and so on
            ind_neu_stim    = cell2mat(stimdat(:,2)) == neucond;
            
            ind_rew_fb    = cell2mat(fbdat(:,2)) == rewcond; % marking the trials that are rewards, and so on
            ind_neu_fb    = cell2mat(fbdat(:,2)) == neucond;
            
            % assign data: now we assign our actual pupill time-series data
            stim_eye        = plotdat_stim.eye_rv.trial'; % raw, all survived trial data
            stim_rew        = stim_eye(ind_rew_stim,1); % only reward trials
            stim_neu        = stim_eye(ind_neu_stim,1); % only non-reward trials
            
            fb_eye        = plotdat_fb.eye_rv.trial'; % raw, all survived trial data
            fb_rew        = fb_eye(ind_rew_fb,1); % only reward trials
            fb_neu        = fb_eye(ind_neu_fb,1); % only non-reward trials
            
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
            
            plotdat.fb.all = []; % stim
            for i = 1:size(fb_eye,1) % for the number of trials
                plotdat.fb.all(i,:) = fb_eye{i,:}; % create one big matrix
            end
            plotdat.fb.all(emptytrialf,:) = [];
            
            % stat preprocessing fix 19-07-2018:
            % first, standardise within trials, then baseline correct.
            plotdat_z = []; clear i
            for i = 1:size(plotdat.stim.all,1)
                mean_fb_rewdat = nanmean(plotdat.stim.all,2); stddat_neu = nanstd(plotdat.stim.all');
                plotdat_z.stim.all = plotdat.stim.all - repmat(mean_fb_rewdat,1,size(plotdat.stim.all,2));
                plotdat_z.stim.all = plotdat.stim.all ./ repmat(stddat_neu',1,size(plotdat.stim.all,2));
                
                bsldat = nanmean(plotdat_z.stim.all(:,timewindow_bsl(1):timewindow_bsl(2))');
                plotdat_preproc.stim.all = plotdat_z.stim.all-repmat(bsldat',1,size(plotdat_z.stim.all,2));
            end
            
            for i = 1:size(plotdat.fb.all,1)
                mean_fb_rewdat = nanmean(plotdat.fb.all,2); stddat_neu = nanstd(plotdat.fb.all');
                plotdat_z.fb.all = plotdat.fb.all - repmat(mean_fb_rewdat,1,size(plotdat.fb.all,2));
                plotdat_z.fb.all = plotdat.fb.all ./ repmat(stddat_neu',1,size(plotdat.fb.all,2));
                
                bsldat = nanmean(plotdat_z.fb.all(:,timewindow_bsl(1):timewindow_bsl(2))');
                plotdat_preproc.fb.all = plotdat_z.fb.all-repmat(bsldat',1,size(plotdat_z.fb.all,2));
            end
            
            plotdat_preproc.stim.rew     = plotdat_preproc.stim.all(ind_rew_stim,:);
            plotdat_preproc.stim.neu     = plotdat_preproc.stim.all(ind_neu_stim,:);
            
            plotdat_preproc.fb.rew     = plotdat_preproc.fb.all(ind_rew_fb,:);
            plotdat_preproc.fb.neu     = plotdat_preproc.fb.all(ind_neu_fb,:);
            
            % average over the corresponding trials
            plotavg.stim.all     = mean(plotdat_preproc.stim.all,1);
            plotavg.stim.rew     = mean(plotdat_preproc.stim.rew,1);
            plotavg.stim.neu     = mean(plotdat_preproc.stim.neu,1);
            
            plotavg.fb.all       = mean(plotdat_preproc.fb.all,1);
            plotavg.fb.rew       = mean(plotdat_preproc.fb.rew,1);
            plotavg.fb.neu       = mean(plotdat_preproc.fb.neu,1);
            
            save([ path_clean 're' num2str(d) num2str(VPs(i1)) '_preproc_dat.mat'], 'plotdat_preproc','plotavg');
            
            
            fprintf('\n Done \n')
            
            %% draw
            % legends do not work -- some errors with cvx, i think
            
                      
            fb_rewdat = plotdat_preproc.fb.rew;
            fb_neudat = plotdat_preproc.fb.neu;
            
            mean_fb_rewdat = nanmean(fb_rewdat,1);stddat_rew = nanstd(fb_rewdat);
            fb_rewdat = fb_rewdat-repmat(mean_fb_rewdat,1,size(fb_rewdat,2));fb_rewdat = fb_rewdat./repmat(stddat_rew',1,size(fb_rewdat,2));
            
            for trl = 1:size(fb_rewdat,1)
                fb_rewdat(trl,:) = fb_rewdat(trl,:)-nanmean(fb_rewdat(trl,1:200));
            end
           
            mean_fb_neudat = nanmean(fb_neudat,1);stddat_neu = nanstd(fb_neudat);
            fb_neudat = fb_neudat-repmat(mean_fb_neudat,1,size(fb_neudat,2));fb_neudat = fb_neudat./repmat(stddat_neu',1,size(fb_neudat,2));
           
            for trl = 1:size(fb_neudat,1)
                fb_neudat(trl,:) = fb_neudat(trl,:)-nanmean(fb_neudat(trl,1:200));
            end
            
%             cd(path_fig) % saving path for the figures
%             close all
            
            maxy = 3; miny = -3;
            onsetbar=200.*ones(1,length(time{1,1})); y=linspace(miny,maxy,numel(onsetbar));
            errorbar = @ (data) (nanstd(data)./sqrt(length(data)));
            
            % fancier: shadedErrorBar
            % if you would like to know about what to put in as arguments,
            % please type: doc shadedErrorBar
            
            % stimulus seg
            figure;
%             subplot(1,2,1);
            plot(onsetbar,y,'g-','linewidth',2)
            shadedErrorBar((1:length(time{1,1})),fb_rewdat,stddat_rew,'lineprops',{'r-o','markerfacecolor','r'}); hold on;
            shadedErrorBar((1:length(time{1,1})),fb_neudat,stddat_neu,'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%             H1 = ttest2(plotdat_preproc.stim.rew(:,201:end),plotdat_preproc.stim.neu(:,201:end)); sigbar1 = H1*-0.9; sigbar1(find(sigbar1==0))=NaN; plot(sigbar1,'k-o','LineWidth',3); % this bit is for calculating significance
            xlim([1 length(time{1,1})]); ylim([miny maxy]);% xticklabels([0 500 1000 1500 2000 2500]);
            xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
            ylabel('zscores','Fontsize',20,'Fontweight','bold'); clear H sigbar
            title('stimulus','Fontsize',20,'Fontweight','bold')
%             
%             % feedback seg
%             subplot(1,2,2);
%             plot(onsetbar,y,'g-','linewidth',2)
%             shadedErrorBar((1:length(time{1,1})),plotdat_preproc.fb.rew,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
%             shadedErrorBar((1:length(time{1,1})),plotdat_preproc.fb.neu,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
%             H1 = ttest2(plotdat_preproc.fb.rew(:,201:end),plotdat_preproc.fb.neu(:,201:end)); sigbar1 = H1*-0.9; sigbar1(find(sigbar1==0))=NaN; plot(sigbar1,'k-o','LineWidth',3); % this bit is for calculating significance
%             xlim([1 length(time{1,1})]); ylim([miny maxy]);% xticklabels([0 500 1000 1500 2000 2500]);
%             xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
%             ylabel('zscores','Fontsize',20,'Fontweight','bold'); clear H sigbar
%             title('feedback','Fontsize',20,'Fontweight','bold')
%             
%             suptitle(['reward(red) vs. neutral(blue) : session' num2str(d)])
%             
% %             
%             % save
%             fig = gcf;
%             fig.PaperPositionMode = 'auto';
%             print([num2str(VPs(i1)) '_maintask_Reward(red)_vs_Neutral(blue)_session' num2str(d) ],'-dpdf','-r0')
%             clear fig H sigbar
%             close all
%             
%             
            
            fprintf('\n ID %d: PLOTTING DONE\n', VPs(i1))
            
            keep session memtest VPs ...
                path_gen path_par path_dat path_behav path_clean path_fig path_stim ...
                fname_raweye fname_beh expdat i1 onsetbar maxy miny time tag d
            
        end
    end
end

fprintf('\n main task plotting done\n')


%% memory test

%% preparation

clear;clc;close all;
% warning('off')
% addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20180613'))
path_gen    = '/Users/yeojin/Desktop/E_data/';
path_par    = [ path_gen 'EA_raw/' ];
path_behav  = [ path_par 'EAC_behav/pilot_3T_2sess/tmp/' ];
path_clean  = [ path_gen 'EB_cleaned/EBB_pupil/pilot_3T/' ];
path_fig    = [ '/Users/yeojin/Desktop/C_writings/CB_figures/pilot_3T_2sess/MainTask/zscores/' ];
% path_stim   = [ path_par 'EA_raw/EAA_pupil/stiminfo_all_20181004.mat' ];

% create file list and VP ID list
VPs = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221]; % 2113(D1), 2207(D2), 2116(D1) 2217(D2) 2218(D1&D2) -> bad data / 2201, 2111 -> dropout
tag = [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 0 2; 1 0; 1 2; 1 2; 1 2; 1 2];
d1m = [1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

fname_raweye_main = []; fname_raweye_immediate = []; fname_raweye_delayed = []; fname_beh =[]; expdat = [];
for i1 = 1:length(VPs)
    for d = 1:2
        if tag(i1,d) == 0
            fname_raweye_main{i1,d}  = []; % rewardtask-encoding (main task)
            fname_raweye_immediate{i1,d}  = []; % rewardtask-immediate recall
            fname_raweye_delayed{i1,d}  = []; % rewardtask-delayed recall
            fname_beh{i1,d} = [];
            expdat{i1,d} = [];
        elseif tag(i1,d) == 1
            fname_raweye_main{i1,d}  = [ 're' num2str(d) num2str(VPs(i1)) '.asc' ]; % rewardtask-encoding (main task)
            if d1m(i1,1) == 0
                fname_raweye_immediate{i1,d}  = [];
            else
                fname_raweye_immediate{i1,d}  = [ 'ri' num2str(d) num2str(VPs(i1)) '.asc' ]; % rewardtask-immediate recall
            end
            if d1m(i1,2) == 0
                fname_raweye_delayed{i1,d} = [];
            else
                fname_raweye_delayed{i1,d}  = [ 'rd' num2str(d) num2str(VPs(i1)) '.asc' ]; % rewardtask-delayed recall
            end
            fname_beh{i1,d}     = [ num2str(VPs(i1)) '_' num2str(d) '.mat' ];
            expdat{i1,d} = load([path_behav fname_beh{i1,d}]);
        elseif tag(i1,d) == 2
            fname_raweye_main{i1,d}  = [ 're' num2str(d) num2str(VPs(i1)) '.asc' ]; % rewardtask-encoding (main task)
            if d2m(i1,1) == 0
                fname_raweye_immediate{i1,d}  = [];
            else
                fname_raweye_immediate{i1,d}  = [ 'ri' num2str(d) num2str(VPs(i1)) '.asc' ]; % rewardtask-immediate recall
            end
            if d2m(i1,2) == 0
                fname_raweye_delayed{i1,d} = [];
            else
                fname_raweye_delayed{i1,d}  = [ 'rd' num2str(d) num2str(VPs(i1)) '.asc' ]; % rewardtask-delayed recall
            end
            fname_beh{i1,d} = [ num2str(VPs(i1)) '_' num2str(d) '.mat' ];
            expdat{i1,d}    = load([path_behav fname_beh{i1,d}]);
        end
    end
    
end

%% plot and create masks

for i1 = 1:length(VPs)
    
    for d = 1:2
        
        if tag(i1,d) == 0
            fprintf('\nID %d Day %d Skipped\n',VPs(i1),d) 
        else
            
            for m = 1:2
                % call in cleaned data
                if eval(['d' num2str(d) 'm(i1,m) == 1'])
                    fprintf('\n***immediate***\n')
                    plotdat_imm     = load([path_clean 'ri' num2str(d) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
                    
                    % immediate segment!
                    % clean the null trials and rejected trials: breaks
                    rmrow_imm    = find(plotdat_imm.eye_rv.immediate.trl_raus==1); % rejected trials from preprocessing
                    emptymeans   = cellfun(@nanmean,plotdat_imm.eye_rv.trial);
                    emptytrials  = find(emptymeans==0);
                    
                    % stimulus informations: filename, reward/no-reward, memorability
                    eval(['immdat = expdat{i1,d}.dat.day' num2str(d) '.memorytest.immediate.results.trl;']);
                    immdat(isempty(immdat)) = {NaN}; % this line, because if you convert data directly from cells to double, it automatically removes empty cells, which in our case are breaks.
                    
                    % remove breaks and rejected trials
                    immdat(rmrow_imm,:) = []; % and now, we remove the rows of trials we rejected during data cleaning
                    immdat(emptytrials,:) = [];
            
                elseif eval(['d' num2str(d) 'm(i1,m) == 2'])
                    fprintf('\n***delayed***\n')
                    plotdat_del     = load([path_clean 'rd' num2str(d) num2str(VPs(i1)) '_eye_rv200_stim.mat']);
                    
                    % delayed segment!
                    % clean the null trials and rejected trials: breaks
                    rmrow_del   = find(plotdat_del.eye_rv.delayed.trl_raus==1); % rejected trials from preprocessing
                    emptymeans   = cellfun(@nanmean,plotdat_del.eye_rv.trial);
                    emptytrialf  = find(emptymeans==0);
                    
                    % stimulus informations: filename, reward/no-reward, memorability
                    eval(['deldat = expdat{i1,d}.dat.day' num2str(d) '.memorytest.delayed.results.trl;']);
                    deldat(isempty(deldat)) = {NaN}; % this line, because if you convert data directly from cells to double, it automatically removes empty cells, which in our case are breaks.
                    
                    % remove breaks and rejected trials
                    deldat(rmrow_del,:) = []; % and now, we remove the rows of trials we rejected during data cleaning
                    deldat(emptytrialf,:) = [];
                    
                else
                    if m == 1; errorname = 'immediate'; else errorname = 'delayed';end
                    fprintf(['no data for day %d for ID %d %s test\n'], d, VPs(i1),errorname)
                end
            end
            
            fprintf('\n Preparation done \n')
            
            %% housekeeping
            
            % indexing
            eval(['rewcond = str2num(expdat{i1,d}.dat.day' num2str(d) '.RewardCategory(9))']); % rewarding condition, human, was indexed as 2 in our data structure
            if rewcond==1
                neucond = 2;
            elseif rewcond == 2
                neucond = 1;
            elseif rewcond == 3
                neucond = 4;
            elseif rewcond == 4
                neucond = 3;
            end
            
            for m = 1:2
                % call in cleaned data
                if eval(['d' num2str(d) 'm(i1,m) == 1'])
                    ind_rew_imm    = cell2mat(immdat(:,2)) == rewcond; % marking the trials that are rewards, and so on
                    ind_neu_imm    = cell2mat(immdat(:,2)) == neucond;
                    ind_old_imm    = cell2mat(immdat(:,4)) == 1; % marking the trials that are rewards, and so on
                    ind_new_imm    = cell2mat(immdat(:,4)) == 2;
                    
                    % assign data: now we assign our actual pupill time-series data
                    imm_eye        = plotdat_imm.eye_rv.trial'; % raw, all survived trial data
                    imm_rew        = imm_eye(ind_rew_imm,1); % only reward trials
                    imm_neu        = imm_eye(ind_neu_imm,1); % only non-reward trials
                    imm_old        = imm_eye(ind_old_imm,1);
                    imm_new        = imm_eye(ind_new_imm,1);
                    
                    % gather and standardise
                    time = plotdat_imm.eye_rv.time;
                    timewindow_bsl = [1 200]; % -200ms ~ 0ms
                    plotdat        = [];
                    
                    % first, we squeeze each trial and time-series into one matrix for the sake of convenience
                    
                    plotdat.imm.all = []; % stim
                    for i = 1:size(imm_eye,1) % for the number of trials
                        plotdat.imm.all(i,:) = imm_eye{i,:}; % create one big matrix
                    end
                    plotdat.imm.all(emptytrials,:) = [];
                    
                    % stat preprocessing fix 19-07-2018:
                    % first, standardise within trials, then baseline correct.
                    plotdat_z = []; clear i
                    for i = 1:size(plotdat.imm.all,1)
                        meandat = nanmean(plotdat.imm.all,2); stddat = nanstd(plotdat.imm.all');
                        plotdat_z.imm.all = plotdat.imm.all - repmat(meandat,1,size(plotdat.imm.all,2));
                        plotdat_z.imm.all = plotdat.imm.all ./ repmat(stddat',1,size(plotdat.imm.all,2));
                        
                        bsldat = nanmean(plotdat_z.imm.all(:,timewindow_bsl(1):timewindow_bsl(2))');
                        plotdat_preproc.imm.all = plotdat_z.imm.all-repmat(bsldat',1,size(plotdat_z.imm.all,2));
                    end
                    
                    plotdat_preproc.imm.rew     = plotdat_preproc.imm.all(ind_rew_imm,:);
                    plotdat_preproc.imm.neu     = plotdat_preproc.imm.all(ind_neu_imm,:);
                    plotdat_preproc.imm.old     = plotdat_preproc.imm.all(ind_old_imm,:);
                    plotdat_preproc.imm.new     = plotdat_preproc.imm.all(ind_new_imm,:);
                    
                    % average over the corresponding trials
                    plotavg.imm.all     = mean(plotdat_preproc.imm.all,1);
                    plotavg.imm.rew     = mean(plotdat_preproc.imm.rew,1);
                    plotavg.imm.neu     = mean(plotdat_preproc.imm.neu,1);
                    plotavg.imm.old     = mean(plotdat_preproc.imm.old,1);
                    plotavg.imm.new     = mean(plotdat_preproc.imm.new,1);
                    
                elseif eval(['d' num2str(d) 'm(i1,m) == 2'])
                    ind_rew_del    = cell2mat(deldat(:,2)) == rewcond; % marking the trials that are rewards, and so on
                    ind_neu_del    = cell2mat(deldat(:,2)) == neucond;
                    ind_old_del    = cell2mat(deldat(:,4)) == 1; % marking the trials that are rewards, and so on
                    ind_new_del    = cell2mat(deldat(:,4)) == 2;
                    
                    del_eye        = plotdat_del.eye_rv.trial'; % raw, all survived trial data
                    del_rew        = del_eye(ind_rew_del,1); % only reward trials
                    del_neu        = del_eye(ind_neu_del,1); % only non-reward trials
                    del_old        = del_eye(ind_old_del,1);
                    del_new        = del_eye(ind_new_del,1);
                    
                    plotdat.del.all = []; % stim
                    for i = 1:size(del_eye,1) % for the number of trials
                        plotdat.del.all(i,:) = del_eye{i,:}; % create one big matrix
                    end
                    plotdat.del.all(emptytrialf,:) = [];
                    
                    for i = 1:size(plotdat.del.all,1)
                        meandat = nanmean(plotdat.del.all,2); stddat = nanstd(plotdat.del.all');
                        plotdat_z.del.all = plotdat.del.all - repmat(meandat,1,size(plotdat.del.all,2));
                        plotdat_z.del.all = plotdat.del.all ./ repmat(stddat',1,size(plotdat.del.all,2));
                        
                        bsldat = nanmean(plotdat_z.del.all(:,timewindow_bsl(1):timewindow_bsl(2))');
                        plotdat_preproc.del.all = plotdat_z.del.all-repmat(bsldat',1,size(plotdat_z.del.all,2));
                    end
                    plotdat_preproc.del.rew     = plotdat_preproc.del.all(ind_rew_del,:);
                    plotdat_preproc.del.neu     = plotdat_preproc.del.all(ind_neu_del,:);
                    plotdat_preproc.del.old     = plotdat_preproc.del.all(ind_old_del,:);
                    plotdat_preproc.del.new     = plotdat_preproc.del.all(ind_new_del,:);
                    
                    plotavg.del.all       = mean(plotdat_preproc.del.all,1);
                    plotavg.del.rew       = mean(plotdat_preproc.del.rew,1);
                    plotavg.del.neu       = mean(plotdat_preproc.del.neu,1);
                    plotavg.del.old       = mean(plotdat_preproc.del.old,1);
                    plotavg.del.new       = mean(plotdat_preproc.del.new,1);
                    
                else
                    if m == 1; errorname = 'immediate'; else errorname = 'delayed';end
                    fprintf(['no data for day %d for ID %d %s test\n'], d, VPs(i1),errorname)
                end
            end
            
            save([ path_clean 'rr' num2str(d) num2str(VPs(i1)) '_preproc_dat.mat'], 'plotdat_preproc','plotavg');
            
            fprintf('\n Done \n')
            
            %% draw
            % legends do not work -- some errors with cvx, i think
            
            cd(path_fig) % saving path for the figures
            close all
            
            maxy = 3; miny = -3;
            onsetbar=200.*ones(1,length(time{1,1})); y=linspace(miny,maxy,numel(onsetbar));
            errorbar = @ (data) (nanstd(data)./sqrt(length(VPs)));
            
            % fancier: shadedErrorBar
            % if you would like to know about what to put in as arguments,
            % please type: doc shadedErrorBar
            
            % immediate seg
            figure;
            if eval(['d' num2str(d) 'm(i1,1) == 1'])
                fprintf(['drawing immediate memory test graph\n'])
            subplot(1,2,1);
            plot(onsetbar,y,'g-','linewidth',2)
            shadedErrorBar((1:length(time{1,1})),plotdat_preproc.imm.rew,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
            shadedErrorBar((1:length(time{1,1})),plotdat_preproc.imm.neu,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
            H1 = ttest2(plotdat_preproc.imm.rew(:,201:end),plotdat_preproc.imm.neu(:,201:end)); sigbar1 = H1*-0.9; sigbar1(find(sigbar1==0))=NaN; plot(sigbar1,'k-o','LineWidth',3); % this bit is for calculating significance
            xlim([1 length(time{1,1})]); ylim([miny maxy]);% xticklabels([0 500 1000 1500 2000 2500]);
            xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
            ylabel('zscores','Fontsize',20,'Fontweight','bold'); clear H sigbar
            title('immediate','Fontsize',20,'Fontweight','bold')
            end
            
            % delayed seg
            if eval(['d' num2str(d) 'm(i1,2) == 2'])
                fprintf(['drawing delayed memory test graph\n'])
            subplot(1,2,2);
            plot(onsetbar,y,'g-','linewidth',2)
            shadedErrorBar((1:length(time{1,1})),plotdat_preproc.del.rew,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
            shadedErrorBar((1:length(time{1,1})),plotdat_preproc.del.neu,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
            H1 = ttest2(plotdat_preproc.del.rew(:,201:end),plotdat_preproc.del.neu(:,201:end)); sigbar1 = H1*-0.9; sigbar1(find(sigbar1==0))=NaN; plot(sigbar1,'k-o','LineWidth',3); % this bit is for calculating significance
            xlim([1 length(time{1,1})]); ylim([miny maxy]);% xticklabels([0 500 1000 1500 2000 2500]);
            xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
            ylabel('zscores','Fontsize',20,'Fontweight','bold'); clear H sigbar
            title('delayed','Fontsize',20,'Fontweight','bold')
            end
            
            suptitle(['reward(red) vs. neutral(blue) : session' num2str(d)])
            
            
            % save
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print([num2str(VPs(i1)) '_memorytest_Reward(red)_vs_Neutral(blue)_session' num2str(d) ],'-dpdf','-r0')
            clear fig H sigbar
            close all
            
            
            
            % immediate seg
            figure;
            if eval(['d' num2str(d) 'm(i1,1) == 1'])
                fprintf(['drawing immediate memory test graph\n'])
            subplot(1,2,1);
            plot(onsetbar,y,'g-','linewidth',2)
            shadedErrorBar((1:length(time{1,1})),plotdat_preproc.imm.old,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
            shadedErrorBar((1:length(time{1,1})),plotdat_preproc.imm.new,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
            H1 = ttest2(plotdat_preproc.imm.old(:,201:end),plotdat_preproc.imm.new(:,201:end)); sigbar1 = H1*-0.9; sigbar1(find(sigbar1==0))=NaN; plot(sigbar1,'k-o','LineWidth',3); % this bit is for calculating significance
            xlim([1 length(time{1,1})]); ylim([miny maxy]);% xticklabels([0 500 1000 1500 2000 2500]);
            xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
            ylabel('zscores','Fontsize',20,'Fontweight','bold'); clear H sigbar
            title('immediate','Fontsize',20,'Fontweight','bold')
            end
            
            % delayed seg
            if eval(['d' num2str(d) 'm(i1,2) == 2'])
                fprintf(['drawing delayed memory test graph\n'])
            subplot(1,2,2);
            plot(onsetbar,y,'g-','linewidth',2)
            shadedErrorBar((1:length(time{1,1})),plotdat_preproc.del.old,{@nanmean,errorbar},'lineprops',{'r-o','markerfacecolor','r'}); hold on;
            shadedErrorBar((1:length(time{1,1})),plotdat_preproc.del.new,{@nanmean,errorbar},'lineprops',{'b-o','markerfacecolor','b'}); hold on;
            H1 = ttest2(plotdat_preproc.del.old(:,201:end),plotdat_preproc.del.new(:,201:end)); sigbar1 = H1*-0.9; sigbar1(find(sigbar1==0))=NaN; plot(sigbar1,'k-o','LineWidth',3); % this bit is for calculating significance
            xlim([1 length(time{1,1})]); ylim([miny maxy]);% xticklabels([0 500 1000 1500 2000 2500]);
            xlabel('time(msec)','Fontsize',20,'Fontweight','bold');
            ylabel('zscores','Fontsize',20,'Fontweight','bold'); clear H sigbar
            title('delayed','Fontsize',20,'Fontweight','bold')
            end
            
            suptitle(['old(red) vs. new(blue) : session' num2str(d)])
            
            
            % save
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print([num2str(VPs(i1)) '_memorytest_Old(red)_vs_New(blue)_session' num2str(d) ],'-dpdf','-r0')
            clear fig H sigbar
            close all
            
            
            
            fprintf('\n ID %d: PLOTTING DONE\n', VPs(i1))
            
            keep session memtest VPs ...
                path_gen path_par path_dat path_behav path_clean path_fig path_stim ...
                fname_raweye fname_beh expdat i1 onsetbar maxy miny time tag d m d1m d2m 
            
        end
    end
end



fprintf('\n memory test plotting done\n')
