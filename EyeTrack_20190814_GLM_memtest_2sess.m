%% Run GLM analysis on pupil data

% description from Dorothea's script

% ----------------------------------------------------------------------- %
% This function computes a GLM regressig behavioural parameters on pupil
% data.
%
% @pupil: pupil data in the same structure as fieldtrip
% @beh: structure with fields containg regressors (e.g. beh.condition =
% randn(n trials,1)
% @reg_names: names of regressors to be used. Must be the same names as the
% fields in the beh structure (e.g. {'condition'})
% @contrast: contrast matrix with each row defining a new contrast and the
% columns define the regressor number (as specified in reg_names +
% intercept as last column). If empty, a separate contrast for each
% regressor will be computed.
% @con_name: name of each of the contrast - used for the output structure
% options:
%       .rejBad: [0,1] reject bad trials specified in pupil.bad (default: 1)
%       .bl_corr: [0,1] baseline-correct each trial (default: 0)
%       .bl_win: [start end] timewindow (in s) to correct for baseline
%       effect. Default: every datapoint before t=0
%       .verbose: plot design matrix (default: 1)
% ----------------------------------------------------------------------- %



%% work log

%   25-02-2019  created the script

%% stimulus bits

%% preparation

clear;clc;close all;
% rng(0)
addpath(genpath('/Users/yeojin/Desktop/B_scripts/BC_functions'))

pilotNr     = 1;

% paths
path_gen    = '/Users/yeojin/Desktop/';
path_par    = [ path_gen 'E_data/' ];
path        = [ path_par 'EA_raw/EAC_behav/pilot_3T_2sess/tmp/' ]; % maybe the mat file with the behavioural results
pathdata    = [ path_par 'EB_cleaned/EBB_pupil/pilot_3T/' ];
pathsave    = [ path_gen 'C_writings/CB_figures/pilot_3T_2sess/MainTask/GLMs/' ]; % where do you save your figures?

% variables
IDs     = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221]; % 2113(D1), 2207(D2), 2116(D1) 2117(D2) -> bad data / 2201, 2111 -> dropout
sessions= [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 0 2; 1 0; 1 2; 1 2; 1 2; 1 2];
d1m     = [1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m     = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
% maxtr = 120;

% load behavioural results for each participants
expdat = [];
for ids = 1:length(IDs)
    for d = 1:2
        if sessions(ids,d) == 0
            expdat{ids,d} = {};
        else
            expdat{ids,d}    = load([path num2str(IDs(ids)) '_' num2str(d) '.mat']);
        end
    end
end

% stat presets
randpermute = 1;  % do the random permutation? (1 = yes! / 0 = nope)
length_randpermute_indiv = 100; length_randpermute = 1000;
sign_randpermute = .05; % less than 5 percent of randpermute results have to be learger than age effect
signcut = 5; % significance level for group contrast
tw = [201:2701];

% call in indice
% indice = load([pathdata 'pilot' num2str(pilotNr) '_index.mat']);


%% run: immediate recall tests, INDIVIDUAL level

storedat_i = [];

for id = 1:length(IDs)
    
    for d = 1:2
        if sessions(id,d) == 0
            fprintf(['\nno data for %d for day %d\n'], IDs(id), d)
            storedat_i(id,d).glm = NaN;
            eval(['contrasts_i_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
            eval(['residuals_i_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
            eval(['tstats_i_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
            
        else
            
            if eval(['d' num2str(d) 'm(id,1) == 0'])
                
                fprintf(['\nno data for %d for day %d\n'], IDs(id), d)
                storedat_d(id,d).glm = NaN;
                eval(['contrasts_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
                eval(['residuals_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
                eval(['tstats_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
                
            else
                
            load([pathdata 'ri' num2str(d) num2str(IDs(id)) '_eye_rv200_stim.mat'],'eye_rv');
            
            % make missing trials NaNs
            alltrials   = eye_rv.cfg.previous.previous.trl(:,1);
            newtrials   = eye_rv.sampleinfo(:,1);
            trlnum      = find(ismember(alltrials,newtrials));
            
            %bsl correct (200 to 0 ms time window)
            for trl = 1:length(eye_rv.trial)
                eye_rv.trial{trl} = eye_rv.trial{trl}-nanmean(eye_rv.trial{trl}(1:200));
            end
            
            % define markers
            eval(['rewcond = str2num(expdat{id,d}.dat.day' num2str(d) '.RewardCategory(9))']) % rewarding condition, human, was indexed as 2 in our data structure
            if rewcond==1
                neucond = 2
            elseif rewcond == 2
                neucond = 1
            elseif rewcond == 3
                neucond = 4
            elseif rewcond == 4
                neucond = 3
            end
            oldcond = 1; newcond = 2;
            imdat = cell2mat(eye_rv.immediate.memstimvec(:,[2 4])); % first column in/out, second old/new
            
            % indexing trial properties (for me)
            indRew      = find(imdat(:,1) == rewcond);
            indNeu      = find(imdat(:,1) == neucond);
            indOld      = find(imdat(:,2) == oldcond);
            indNew      = find(imdat(:,2) == newcond);
            indOldCorr  = find(imdat(:,2) == oldcond & eye_rv.immediate.memresults(:,2) == 1);
            
            % assign dummy variables for GLM analysis
            
%             glmname     = {'REWARD versus NoREWARD'};
%             numContr    = 1;
%             
%             dummie      = zeros(length(trlnum),1);
%             dummie(indRew) = 1;
%             beh.Reward  = zscore(dummie); clear dummie; %%%%%
%             
%             storebehav  = beh;
%             storepupil  = eye_rv;
%             
%                 glmname     = {'Old versus New'};
%                 numContr    = 1;
%             
%                 dummie      = zeros(length(trlnum),1);
%                 dummie(indOld) = 1;
%                 beh.Old  = zscore(dummie); clear dummie; %%%%%
%             
%                 storebehav  = beh;
%                 storepupil  = eye_rv;
                
                glmname     = {'corrects: Old versus New'};
                numContr    = 1;
            
                dummie      = zeros(length(trlnum),1);
                dummie(indOldCorr) = 1;
                beh.Oldcorr  = zscore(dummie); clear dummie; %%%%%
            
                storebehav  = beh;
                storepupil  = eye_rv;
            
            %     glmname     = {'HiMemVsLoMem'};
            %     numContr    = 1;
            %
            %     dummie      = zeros(length(trlnum),1);
            %     dummie(indmH) = 1;
            %     beh.HiMem  = zscore(dummie); clear dummie; %%%%%
            %
            %     storebehav  = beh;
            %     storepupil  = eye_rv;
            
            
            %     glmname     = {'Reward and Memorability:'};
            %     numContr    = 1; % 1: Reward_hiMem 2: Reward_loMem 3: NoReward_hiMem 4: NoReward_lomem
            %
            %     dummie      = zeros(length(trlnum),1);
            %     dummie(indRew_mH) = 1;
            %     beh.RewardHiMem  = zscore(dummie); clear dummie;
            %
            %     dummie      = zeros(length(trlnum),1);
            %     dummie(indRew_mL) = 1;
            %     beh.RewardLoMem= zscore(dummie); clear dummie;
            %
            %     dummie      = zeros(length(trlnum),1);
            %     dummie(indnRew_mH) = 1;
            %     beh.NoRewardHiMem  = zscore(dummie); clear dummie;
            %
            %     dummie      = zeros(length(trlnum),1);
            %     dummie(indnRew_mL) = 1;
            %     beh.NoRewardLoMem= zscore(dummie); clear dummie;
            %
            %     storebehav  = beh;
            %     storepupil  = eye_rv;
            
            
            % build structures for eye_glm (beh is already done above)
            
            pupil       = eye_rv;
            pupil.trial = []; pupil.trial = storepupil.trial;
            pupil.time  = []; pupil.time = storepupil.time;
            pupil.sampleinfo = []; pupil.sampleinfo = storepupil.sampleinfo;
            
            reg_names   = {'Oldcorr'}; % regressor names
            %     reg_names   = {'RewardHiMem'; 'RewardLoMem'; 'NoRewardHiMem'; 'NoRewardLoMem';};
            %     contrast    = [1 -1 0; 1 0 0; 0 1 0; 0 0 1]; d
            con_name    = {'Oldcorr pos'};
            contrast    = 1;
            %     con_name    = {'RewardHiMem'; 'RewardLoMem'; 'NoRewardHiMem'; 'NoRewardLoMem';};
            
            options.rejBad = 0;
            options.bl_corr = 0;
            options.bl_win = [-.2 0];
            
            if id ~= length(IDs)
                clear eye_rv
            elseif id == length(IDs)
            end
            
            % run GLM
            storedat_i(id,d).glm = eye_glm(pupil,beh,reg_names,contrast,con_name,options)
            
            
            % make randpermute for each person
            if randpermute
                
                for permcount = 1:length_randpermute_indiv
                    
                    eval(['lengthbehav = length(beh.' reg_names{1} ');'])
                    ind = randperm(lengthbehav);
                    
                    for c1 = 1:size(reg_names,1)
                        eval(['beh.' reg_names{c1} '= beh.' reg_names{c1} '(ind);'])
                    end
                    
                    dummie = eye_glm(pupil,beh,reg_names,contrast,con_name,options);
                    for stats = 1:size(reg_names,1)+1
                        try
                            eval(['contrasts_i_' num2str(d) '(id,stats,permcount,:) = dummie.con{1,stats};'])
                            eval(['residuals_i_' num2str(d) '(id,stats,permcount,:) = dummie.res{1,stats};'])
                            eval(['tstats_i_' num2str(d) '(id,stats,permcount,:) = dummie.tstat{1,stats};'])
                        catch
                            eval(['contrasts_i_' num2str(d) '(id,stats,permcount,:) = dummie.con{1,stats}(:,2701);'])
                            eval(['residuals_i_' num2str(d) '(id,stats,permcount,:) = dummie.res{1,stats}(:,2701);'])
                            eval(['tstats_i_' num2str(d) '(id,stats,permcount,:) = dummie.tstat{1,stats}(:,2701);'])
                        end
                        
                    end
                    
                    clear ind dummie
                    
                end % close random permutation loop
                
            end % close random permutation calculation
            
            fprintf(['\n*********\n' num2str(IDs(id)) ' done\n*********\n'])
        end
    end
    
    end
end


groupGLM_i_d1.contrasts = contrasts_i_1;
groupGLM_i_d1.residuals = residuals_i_1;
groupGLM_i_d1.tstats    = tstats_i_1;

groupGLM_i_d2.contrasts = contrasts_i_2;
groupGLM_i_d2.residuals = residuals_i_2;
groupGLM_i_d2.tstats    = tstats_i_2;

save([pathdata 'groupGLM_3T_2sess_immediate_' char(glmname) '.mat'],'groupGLM_i_d1', 'groupGLM_i_d2')

% cleaning
groupGLM_i_d1.contrasts(find(sessions(:,1)==0),:,:,:) = [];
groupGLM_i_d1.residuals(find(sessions(:,1)==0),:,:,:) = [];
groupGLM_i_d1.tstats(find(sessions(:,1)==0),:,:,:) = [];

groupGLM_i_d2.contrasts(find(sessions(:,2)==0),:,:,:) = [];
groupGLM_i_d2.residuals(find(sessions(:,2)==0),:,:,:) = [];
groupGLM_i_d2.tstats(find(sessions(:,2)==0),:,:,:) = [];


%% run: immediate task, stimulus bits, GROUP level

% clean up the storedat
storedat_i_d1 = storedat_i(:,1); storedat_i_d1(sessions(:,1)==0 | d1m(:,1)==0) = [];
storedat_i_d2 = storedat_i(:,2); storedat_i_d2(sessions(:,2)==0 | d2m(:,1)==0) = [];
storedat_all = vertcat(storedat_i_d1,storedat_i_d2);

contrasts = vertcat(groupGLM_i_d1.contrasts,groupGLM_i_d2.contrasts);
residuals = vertcat(groupGLM_i_d1.residuals,groupGLM_i_d2.residuals);
tstats = vertcat(groupGLM_i_d1.tstats,groupGLM_i_d2.tstats);

con_name{length(con_name)+1} = 'intercept'; % last contrast is intercept

dm = ones(length(IDs),1); % design matrix for groups (which will come later)

nLC = size(storedat_all(1).glm.con,2); % number of lower level contrasts (each experimental condition)
nSamp = size(storedat_all(1).glm.con{1},2); % number of samples in timeseries

% cg = nan(nHC,nSamp); vg = nan(nHC,nSamp); tg = nan(nHC,nSamp);
% cgb = nan(length(con_name),length_randpermute,nHC,nSamp); vgb = nan(length(con_name),length_randpermute,nHC,nSamp); tgb = nan(length(con_name),length_randpermute,nHC,nSamp);


for con_counter = 1:length(con_name)
    
    for i = 1:length(IDs)
        % save contrast estimates for regressors
        dat1(i,:) = storedat_all(i).glm.con{1,con_counter};
    end
    con_estimates = dat1;
    
    
    % contrast data from each age group for plotting (not the randpermute data)
    time    = pupil.time{1};
    range   = 1:length(time);
    means   = nanmean(dat1(:,:));
    SEs     = nanstd(dat1(:,:))./sqrt(length(IDs));
    
    % 5%, expecting, higher than 0 / assess maximal and minimal value in particular time window
    datdum = squeeze(nanmean(contrasts(:,con_counter,:,:),1));
    for rep = 1:size(contrasts,3)
        datmin(rep) = min(squeeze(datdum(rep,min(tw):max(tw))));
        datmax(rep) = max(squeeze(datdum(rep,min(tw):max(tw))));
    end
    plotsign_indiv = nan(1,size(con_estimates,2));
    plotsign_indiv(find(means(min(tw):max(tw)) < prctile(datmin,10))+min(tw))  = 0;
    plotsign_indiv(find(means(min(tw):max(tw)) > prctile(datmax,90))+min(tw)) = 0;
    
    clear datmin datmax rep datdum
    
    
    %  plot with significance
    darkred = [.9 , 0.1 , 0.1]; darkblue = [.1 , 0.1 , .9];
    figure;
    shadedErrorBar(time(range),means(range),SEs(range),'lineprops',{'k-o','markerfacecolor','k'}); hold on
    plot(time(range),plotsign_indiv(range),'','LineWidth',6,'Color',darkblue);hold off
    xlabel('time(sec)','FontSize',20);ylabel('Effect size (mean +/- SE)','FontSize',20);
    title([glmname{1} ' ' con_name{con_counter}],'FontSize',12);
    set(gca,'FontSize',20);
    axis([min(time(range)) 2.5 -.1 .1]);
    %     axis([min(time(range)) 2.5 min(means(:))-.3  max(means(:))+.3]);
    
    % save the figures
    cd(pathsave)
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    
    saveas(gcf,['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_immediate.fig'])
    print(['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_immediate'],'-dpdf','-fillpage') % save as pdf
    print(['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_immediate'],'-depsc','-tiff','-r0') % save as eps
    
    
    clear dat1 con_estimates fig
    
    
end % close contrast calculation loop

clear id con_counter


%% run: delayed recall tests, INDIVIDUAL level
%% preparation

clear;clc;close all;
% rng(0)
addpath(genpath('/Users/yeojin/Desktop/B_scripts/BC_functions'))

pilotNr     = 1;

% paths
path_gen    = '/Users/yeojin/Desktop/';
path_par    = [ path_gen 'E_data/' ];
path        = [ path_par 'EA_raw/EAC_behav/pilot_3T_2sess/tmp/' ]; % maybe the mat file with the behavioural results
pathdata    = [ path_par 'EB_cleaned/EBB_pupil/pilot_3T/' ];
pathsave    = [ path_gen 'C_writings/CB_figures/pilot_3T_2sess/MainTask/GLMs/' ]; % where do you save your figures?

% variables
IDs     = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221]; % 2113(D1), 2207(D2), 2116(D1) 2117(D2) -> bad data / 2201, 2111 -> dropout
sessions= [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 0 2; 1 0; 1 2; 1 2; 1 2; 1 2];
d1m     = [1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m     = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
% maxtr = 120;

% load behavioural results for each participants
expdat = [];
for ids = 1:length(IDs)
    for d = 1:2
        if sessions(ids,d) == 0
            expdat{ids,d} = {};
        else
            expdat{ids,d}    = load([path num2str(IDs(ids)) '_' num2str(d) '.mat']);
        end
    end
end

% stat presets
randpermute = 1;  % do the random permutation? (1 = yes! / 0 = nope)
length_randpermute_indiv = 100; length_randpermute = 1000;
sign_randpermute = .05; % less than 5 percent of randpermute results have to be learger than age effect
signcut = 10; % significance level for group contrast
tw = [201:2701];

% call in indice
% indice = load([pathdata 'pilot' num2str(pilotNr) '_index.mat']);


%% run: delayed recall tests, INDIVIDUAL level

storedat_d = [];

for id = 1:length(IDs)
    
    for d = 1:2
        if sessions(id,d) == 0
            fprintf(['\nno data for %d for day %d\n'], IDs(id), d)
            storedat_d(id,d).glm = NaN;
            eval(['contrasts_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
            eval(['residuals_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
            eval(['tstats_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
            
        else
            
            if eval(['d' num2str(d) 'm(id,2) == 0'])
                
                fprintf(['\nno data for %d for day %d\n'], IDs(id), d)
                storedat_d(id,d).glm = NaN;
                eval(['contrasts_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
                eval(['residuals_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
                eval(['tstats_d_' num2str(d) '(id,:,:,:) = zeros(1,2,100,2701);'])
                
            else
                
                
                load([pathdata 'rd' num2str(d) num2str(IDs(id)) '_eye_rv200_stim.mat'],'eye_rv');
                
                % make missing trials NaNs
                alltrials   = eye_rv.cfg.previous.previous.trl(:,1);
                newtrials   = eye_rv.sampleinfo(:,1);
                trlnum      = find(ismember(alltrials,newtrials));
                
                %bsl correct (200 to 0 ms time window)
                for trl = 1:length(eye_rv.trial)
                    eye_rv.trial{trl} = eye_rv.trial{trl}-nanmean(eye_rv.trial{trl}(1:200));
                end
                
                % define markers
                eval(['rewcond = str2num(expdat{id,d}.dat.day' num2str(d) '.RewardCategory(9))']) % rewarding condition, human, was indexed as 2 in our data structure
                if rewcond==1
                    neucond = 2
                elseif rewcond == 2
                    neucond = 1
                elseif rewcond == 3
                    neucond = 4
                elseif rewcond == 4
                    neucond = 3
                end
                oldcond = 1; newcond = 2;
                dedat = cell2mat(eye_rv.delayed.memstimvec(:,[2 4])); % first column in/out, second old/new
                
                % indexing trial properties (for me)
                indRew      = find(dedat(:,1) == rewcond);
                indNeu      = find(dedat(:,1) == neucond);
                indOld      = find(dedat(:,2) == oldcond);
                indNew      = find(dedat(:,2) == newcond);
                indOldCorr  = find(dedat(:,2) == oldcond & eye_rv.delayed.memresults(:,2) == 1);
                
                % assign dummy variables for GLM analysis
                
%                 glmname     = {'REWARD versus NoREWARD'};
%                 numContr    = 1;
%                 
%                 dummie      = zeros(length(trlnum),1);
%                 dummie(indRew) = 1;
%                 beh.Reward  = zscore(dummie); clear dummie; %%%%%
%                 
%                 storebehav  = beh;
%                 storepupil  = eye_rv;
                
%                 glmname     = {'Old versus New'};
%                 numContr    = 1;
%             
%                 dummie      = zeros(length(trlnum),1);
%                 dummie(indOld) = 1;
%                 beh.Old  = zscore(dummie); clear dummie; %%%%%
%             
%                 storebehav  = beh;
%                 storepupil  = eye_rv;
                
                glmname     = {'corrects: Old versus New'};
                numContr    = 1;
            
                dummie      = zeros(length(trlnum),1);
                dummie(indOldCorr) = 1;
                beh.Oldcorr  = zscore(dummie); clear dummie; %%%%%
            
                storebehav  = beh;
                storepupil  = eye_rv;
                
                %     glmname     = {'HiMemVsLoMem'};
                %     numContr    = 1;
                %
                %     dummie      = zeros(length(trlnum),1);
                %     dummie(indmH) = 1;
                %     beh.HiMem  = zscore(dummie); clear dummie; %%%%%
                %
                %     storebehav  = beh;
                %     storepupil  = eye_rv;
                
                
                %     glmname     = {'Reward and Memorability:'};
                %     numContr    = 1; % 1: Reward_hiMem 2: Reward_loMem 3: NoReward_hiMem 4: NoReward_lomem
                %
                %     dummie      = zeros(length(trlnum),1);
                %     dummie(indRew_mH) = 1;
                %     beh.RewardHiMem  = zscore(dummie); clear dummie;
                %
                %     dummie      = zeros(length(trlnum),1);
                %     dummie(indRew_mL) = 1;
                %     beh.RewardLoMem= zscore(dummie); clear dummie;
                %
                %     dummie      = zeros(length(trlnum),1);
                %     dummie(indnRew_mH) = 1;
                %     beh.NoRewardHiMem  = zscore(dummie); clear dummie;
                %
                %     dummie      = zeros(length(trlnum),1);
                %     dummie(indnRew_mL) = 1;
                %     beh.NoRewardLoMem= zscore(dummie); clear dummie;
                %
                %     storebehav  = beh;
                %     storepupil  = eye_rv;
                
                
                % build structures for eye_glm (beh is already done above)
                
                pupil       = eye_rv;
                pupil.trial = []; pupil.trial = storepupil.trial;
                pupil.time  = []; pupil.time = storepupil.time;
                pupil.sampleinfo = []; pupil.sampleinfo = storepupil.sampleinfo;
                
                reg_names   = {'Oldcorr'}; % regressor names
                %     reg_names   = {'RewardHiMem'; 'RewardLoMem'; 'NoRewardHiMem'; 'NoRewardLoMem';};
                %     contrast    = [1 -1 0; 1 0 0; 0 1 0; 0 0 1]; d
                con_name    = {'Oldcorr pos'};
                contrast    = 1;
                %     con_name    = {'RewardHiMem'; 'RewardLoMem'; 'NoRewardHiMem'; 'NoRewardLoMem';};
                
                options.rejBad = 0;
                options.bl_corr = 0;
                options.bl_win = [-.2 0];
                
                if id ~= length(IDs)
                    clear eye_rv
                elseif id == length(IDs)
                end
                
                % run GLM
                storedat_d(id,d).glm = eye_glm(pupil,beh,reg_names,contrast,con_name,options)
                
                
                % make randpermute for each person
                if randpermute
                    
                    for permcount = 1:length_randpermute_indiv
                        
                        eval(['lengthbehav = length(beh.' reg_names{1} ');'])
                        ind = randperm(lengthbehav);
                        
                        for c1 = 1:size(reg_names,1)
                            eval(['beh.' reg_names{c1} '= beh.' reg_names{c1} '(ind);'])
                        end
                        
                        dummie = eye_glm(pupil,beh,reg_names,contrast,con_name,options);
                        for stats = 1:size(reg_names,1)+1
                            try
                                eval(['contrasts_d_' num2str(d) '(id,stats,permcount,:) = dummie.con{1,stats};'])
                                eval(['residuals_d_' num2str(d) '(id,stats,permcount,:) = dummie.res{1,stats};'])
                                eval(['tstats_d_' num2str(d) '(id,stats,permcount,:) = dummie.tstat{1,stats};'])
                            catch
                                eval(['contrasts_d_' num2str(d) '(id,stats,permcount,:) = dummie.con{1,stats}(:,2701);'])
                                eval(['residuals_d_' num2str(d) '(id,stats,permcount,:) = dummie.res{1,stats}(:,2701);'])
                                eval(['tstats_d_' num2str(d) '(id,stats,permcount,:) = dummie.tstat{1,stats}(:,2701);'])
                            end
                            
                        end
                        
                        clear ind dummie
                        
                    end % close random permutation loop
                    
                end % close random permutation calculation
                
                fprintf(['\n*********\n' num2str(IDs(id)) ' done\n*********\n'])
            end
        end
    end
end


groupGLM_d_d1.contrasts = contrasts_d_1;
groupGLM_d_d1.residuals = residuals_d_1;
groupGLM_d_d1.tstats    = tstats_d_1;

groupGLM_d_d2.contrasts = contrasts_d_2;
groupGLM_d_d2.residuals = residuals_d_2;
groupGLM_d_d2.tstats    = tstats_d_2;

save([pathdata 'groupGLM_3T_2sess_delayed_' char(glmname) '.mat'],'groupGLM_d_d1', 'groupGLM_d_d2')

% cleaning
groupGLM_d_d1.contrasts(find(sessions(:,1)==0),:,:,:) = [];
groupGLM_d_d1.residuals(find(sessions(:,1)==0),:,:,:) = [];
groupGLM_d_d1.tstats(find(sessions(:,1)==0),:,:,:) = [];

groupGLM_d_d2.contrasts(find(sessions(:,2)==0),:,:,:) = [];
groupGLM_d_d2.residuals(find(sessions(:,2)==0),:,:,:) = [];
groupGLM_d_d2.tstats(find(sessions(:,2)==0),:,:,:) = [];


%% run: delayed task, stimulus bits, GROUP level

% clean up the storedat
storedat_d_d1 = storedat_d(:,1); storedat_d_d1(sessions(:,1)==0 | d1m(:,2)==0) = [];
storedat_d_d2 = storedat_d(:,2); storedat_d_d2(sessions(:,2)==0 | d2m(:,2)==0) = [];
storedat_all = vertcat(storedat_d_d1,storedat_d_d2);

contrasts = vertcat(groupGLM_d_d1.contrasts,groupGLM_d_d2.contrasts);
residuals = vertcat(groupGLM_d_d1.residuals,groupGLM_d_d2.residuals);
tstats = vertcat(groupGLM_d_d1.tstats,groupGLM_d_d2.tstats);

con_name{length(con_name)+1} = 'intercept'; % last contrast is intercept

dm = ones(length(IDs),1); % design matrix for groups (which will come later)

nLC = size(storedat_all(1).glm.con,2); % number of lower level contrasts (each experimental condition)
nSamp = size(storedat_all(1).glm.con{1},2); % number of samples in timeseries

% cg = nan(nHC,nSamp); vg = nan(nHC,nSamp); tg = nan(nHC,nSamp);
% cgb = nan(length(con_name),length_randpermute,nHC,nSamp); vgb = nan(length(con_name),length_randpermute,nHC,nSamp); tgb = nan(length(con_name),length_randpermute,nHC,nSamp);


for con_counter = 1:length(con_name)
    
    for i = 1:length(storedat_all)
        % save contrast estimates for regressors
        try
        dat1(i,:) = storedat_all(i).glm.con{1,con_counter};
        catch
            dat1(i,:) = storedat_all(i).glm.con{1,con_counter}(1,1:2701);
        end
    end
    con_estimates = dat1;
    
    
    % contrast data from each age group for plotting (not the randpermute data)
    time    = pupil.time{1};
    range   = 1:length(time);
    means   = nanmean(dat1(:,:));
    SEs     = nanstd(dat1(:,:))./sqrt(length(IDs));
    
    % 5%, expecting, higher than 0 / assess maximal and minimal value in particular time window
    datdum = squeeze(nanmean(contrasts(:,con_counter,:,:),1));
    for rep = 1:size(contrasts,3)
        datmin(rep) = min(squeeze(datdum(rep,min(tw):max(tw))));
        datmax(rep) = max(squeeze(datdum(rep,min(tw):max(tw))));
    end
    plotsign_indiv = nan(1,size(con_estimates,2));
    plotsign_indiv(find(means(min(tw):max(tw)) < prctile(datmin,10))+min(tw))  = 0;
    plotsign_indiv(find(means(min(tw):max(tw)) > prctile(datmax,90))+min(tw)) = 0;
    
    clear datmin datmax rep datdum
    
    
    %  plot with significance
    darkred = [.9 , 0.1 , 0.1]; darkblue = [.1 , 0.1 , .9];
    figure;
    shadedErrorBar(time(range),means(range),SEs(range),'lineprops',{'k-o','markerfacecolor','k'}); hold on
    plot(time(range),plotsign_indiv(range),'','LineWidth',6,'Color',darkblue);hold off
    xlabel('time(sec)','FontSize',20);ylabel('Effect size (mean +/- SE)','FontSize',20);
    title([glmname{1} ' ' con_name{con_counter}],'FontSize',12);
    set(gca,'FontSize',20);
    axis([min(time(range)) 2.5 -.1 .1]);
    %     axis([min(time(range)) 2.5 min(means(:))-.3  max(means(:))+.3]);
    
    % save the figures
    cd(pathsave)
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    
    saveas(gcf,['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_delayed.fig'])
    print(['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_delayed'],'-dpdf','-fillpage') % save as pdf
    print(['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_delayed'],'-depsc','-tiff','-r0') % save as eps
    
    
    clear dat1 con_estimates fig
    
    
end % close contrast calculation loop

clear id con_counter


%% run: delayed task, stimulus bits, GROUP level, LESS DAY ONLY

% clean up the storedat
% storedat_d_d1 = storedat_d(:,1); storedat_d_d1(sessions(:,1)==0 | d1m(:,2)==0) = [];
storedat_d_d2 = storedat_d(:,2); storedat_d_d2(sessions(:,2)==0 | d2m(:,2)==0) = [];
storedat_all = storedat_d_d2;

contrasts = groupGLM_d_d2.contrasts;
residuals = groupGLM_d_d2.residuals;
tstats = groupGLM_d_d2.tstats;

con_name{length(con_name)+1} = 'intercept'; % last contrast is intercept

dm = ones(length(IDs),1); % design matrix for groups (which will come later)

nLC = size(storedat_all(1).glm.con,2); % number of lower level contrasts (each experimental condition)
nSamp = size(storedat_all(1).glm.con{1},2); % number of samples in timeseries

% cg = nan(nHC,nSamp); vg = nan(nHC,nSamp); tg = nan(nHC,nSamp);
% cgb = nan(length(con_name),length_randpermute,nHC,nSamp); vgb = nan(length(con_name),length_randpermute,nHC,nSamp); tgb = nan(length(con_name),length_randpermute,nHC,nSamp);


for con_counter = 1:length(con_name)
    
    for i = 1:length(storedat_all)
        % save contrast estimates for regressors
        try
        dat1(i,:) = storedat_all(i).glm.con{1,con_counter};
        catch
            dat1(i,:) = storedat_all(i).glm.con{1,con_counter}(1,1:2701);
        end
    end
    con_estimates = dat1;
    
    
    % contrast data from each age group for plotting (not the randpermute data)
    time    = pupil.time{1};
    range   = 1:length(time);
    means   = nanmean(dat1(:,:));
    SEs     = nanstd(dat1(:,:))./sqrt(length(IDs));
    
    % 5%, expecting, higher than 0 / assess maximal and minimal value in particular time window
    datdum = squeeze(nanmean(contrasts(:,con_counter,:,:),1));
    for rep = 1:size(contrasts,3)
        datmin(rep) = min(squeeze(datdum(rep,min(tw):max(tw))));
        datmax(rep) = max(squeeze(datdum(rep,min(tw):max(tw))));
    end
    plotsign_indiv = nan(1,size(con_estimates,2));
    plotsign_indiv(find(means(min(tw):max(tw)) < prctile(datmin,10))+min(tw))  = 0;
    plotsign_indiv(find(means(min(tw):max(tw)) > prctile(datmax,90))+min(tw)) = 0;
    
    clear datmin datmax rep datdum
    
    
    %  plot with significance
    darkred = [.9 , 0.1 , 0.1]; darkblue = [.1 , 0.1 , .9];
    figure;
    shadedErrorBar(time(range),means(range),SEs(range),'lineprops',{'k-o','markerfacecolor','k'}); hold on
    plot(time(range),plotsign_indiv(range),'','LineWidth',6,'Color',darkblue);hold off
    xlabel('time(sec)','FontSize',20);ylabel('Effect size (mean +/- SE)','FontSize',20);
    title([glmname{1} ' ' con_name{con_counter} ' LessDay Only'],'FontSize',12);
    set(gca,'FontSize',20);
    axis([min(time(range)) 2.5 -.1 .1]);
    %     axis([min(time(range)) 2.5 min(means(:))-.3  max(means(:))+.3]);
    
    % save the figures
    cd(pathsave)
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    
    saveas(gcf,['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_delayed_LessdayOnly.fig'])
    print(['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_delayed_LessdayOnly'],'-dpdf','-fillpage') % save as pdf
    print(['GLM_memtest_randpermute_' con_name{con_counter} '_stim_pilot_3T_memtest_delayed_LessdayOnly'],'-depsc','-tiff','-r0') % save as eps
    
    
    clear dat1 con_estimates fig
    
    
end % close contrast calculation loop

clear id con_counter

%% session 1 vs session 2: delayed only

cd(path_gen)

d1l = length(storedat_d_d1); d2l = length(storedat_d_d2);
if d1l > d2l
    datalength = d2l;
    storedat_group    = vertcat(storedat_d_d1(1:datalength),storedat_d_d2);
elseif d1l < d2l
    datalength = d1l;
    storedat_group    = vertcat(storedat_d_d1,storedat_d_d2(1:datalength));
else
    datalength = d1l;
    storedat_group    = vertcat(storedat_d_d1,storedat_d_d2);
end


% last contrast always intercept

group_dm    = [[ones(datalength,1);zeros(datalength,1)],[zeros(datalength,1); ones(datalength,1)]];
group_con   = [1 -1 0 ; 0 0 1]';% looking for day1>day2

nLC     = size(storedat_group(1).glm.con,2); %number of lower level contrasts
nSamp   = size(storedat_group(1).glm.con{1},2); %number of samples in timeseries
nHC     = size(group_con,1); %number of higher level contrasts

cg = nan(nHC,nSamp); vg = nan(nHC,nSamp); tg = nan(nHC,nSamp);
cgb = nan(length(con_name),length_randpermute,nHC,nSamp); vgb = nan(length(con_name),length_randpermute,nHC,nSamp); tgb = nan(length(con_name),length_randpermute,nHC,nSamp);

for contrIND = 1:length(con_name)
    for i = 1:length(storedat_group)%length(IDs)*2
        try
            dat1(i,:) = storedat_group(i).glm.con{1,contrIND}; % save contrast estimates for regressors
        catch
            dat1(i,:) = storedat_group(i).glm.con{1,contrIND}(1,1:2701);
        end
    end
    Y = dat1;
    
    % contrast data from each age group for plotting (not the randpermute data)
    time = pupil.time{1};
    range = 1:length(time);
    mean1more = nanmean(dat1(1:datalength,:));
    se1more = nanstd(dat1(1:datalength,:))./sqrt(datalength);
    mean1less = nanmean(dat1(datalength+1:datalength*2,:));
    se1less = nanstd(dat1(datalength+1:datalength*2,:))./sqrt(datalength);
    
    % calculate group contrast
    group_dm = [];
    group_dm = [[ones(datalength,1);zeros(datalength,1)],[zeros(datalength,1); ones(datalength,1)]];
    %     group_dm = [zscore(group_dm) zscore(group_dm.*repmat(behav_group',1,2))  zscore(behav_group')]; % with interaction effects
    %     group_dm = [zscore(group_dm) zscore(behav_group')];% add behavioral indicator on group level, zscored
    
    [cg,vg,tg] = ols(Y,group_dm,group_con);
    
    % store contrasteffects from group
    store_group_glm{contrIND}.cg = cg;
    store_group_glm{contrIND}.vg = vg;
    store_group_glm{contrIND}.tg = tg;
    
    if randpermute  % shuffle age labels around
        for r = 1:length_randpermute
            ind = randperm(datalength*2);
            Y = Y(ind,:);
            r
            [cgb(contrIND,r,:,:),vgb(contrIND,r,:,:),tgb(contrIND,r,:,:)] = ols(Y,group_dm,group_con);
            % gives one row per group contrast, store per indiviudal contrast (contr)
        end
    end
    
    % store contrasteffects from group, randpermute
    store_group_glm{contrIND}.cgb = squeeze(cgb(contrIND,:,:,:)); % error previously, takinng always age-group difference
    store_group_glm{contrIND}.vgb = squeeze(vgb(contrIND,:,:,:));
    store_group_glm{contrIND}.tgb = squeeze(tgb(contrIND,:,:,:));
    
    % GROUP comparision (5%, expecting immediate > delayed): assess maximal and minimal value in particular time window
    plotsign = nan(2,size(Y,2));
    for grcon = 1:size(group_con,1)% group contrasts
        for rep = 1:size(cgb,2)
            datmin(rep) = min(squeeze(cgb(contrIND,rep,grcon,min(tw):max(tw))));
            datmax(rep) = max(squeeze(cgb(contrIND,rep,grcon,min(tw):max(tw))));
        end
        plotsign(grcon,find(cg(grcon,min(tw):max(tw)) < prctile(datmin,signcut))+min(tw)) = -.2+grcon*.05; % immediate > delayed
        plotsign(grcon,find(cg(grcon,min(tw):max(tw)) > prctile(datmax,100-signcut))+min(tw)) = -.2+grcon*.05;
    end
    
    % get at t value for regression:
    timesign = find(cg(grcon,min(tw):max(tw)) > prctile(datmax,100-signcut))+min(tw);
    max(tg(2,timesign)); % group level second effect
    clear datmin datmax rep
    
    
    % more day: (5%, expecting, higher than 0) assess maximal and minimal value in particular time window
    datdum = squeeze(nanmean(groupGLM_d_d1.contrasts(1:datalength,contrIND,:,:),1));
    for rep = 1:size(groupGLM_d_d1.contrasts,3)
        datmin(rep) = min(squeeze(datdum(rep,min(tw):max(tw))));
        datmax(rep) = max(squeeze(datdum(rep,min(tw):max(tw))));
    end
    plotsignindivMore = nan(1,size(Y,2)); plotsignindivMore(find(mean1more(min(tw):max(tw)) < prctile(datmin,5))+min(tw)) = -.3;
    plotsignindivMore(find(mean1more(min(tw):max(tw)) > prctile(datmax,95))+min(tw)) = -.3;clear datmin datmax rep datdum
    
    % less day: (5%, expecting, higher than 0) assess maximal and minimal value in particular time window
    datdum = squeeze(nanmean(groupGLM_d_d2.contrasts(1:datalength,contrIND,:,:),1));
    for rep = 1:size(groupGLM_d_d2.contrasts,3)
        datmin(rep) = min(squeeze(datdum(rep,min(tw):max(tw))));
        datmax(rep) = max(squeeze(datdum(rep,min(tw):max(tw))));
    end
    plotsignindivLess = nan(1,size(Y,2)); plotsignindivLess(find(mean1less(min(tw):max(tw)) < prctile(datmin,5))+min(tw)) = -.4;
    plotsignindivLess(find(mean1less(min(tw):max(tw)) > prctile(datmax,95))+min(tw)) = -.4;clear datmin datmax rep datdum
    
    
    %  plot with significance
    darkred = [.9 , 0.1 , 0.1]; darkblue = [.1 , 0.1 , .9];
    
    figure,
    shadedErrorBar(time(range),mean1more(range),se1more(range),'lineprops',{'b-o','markerfacecolor','b'});hold on % yas
    shadedErrorBar(time(range),mean1less(range),se1less(range),'lineprops',{'r-o','markerfacecolor','r'});hold on% oas
    plot(time(range),plotsign(1,range),'k','LineWidth',6);hold on
    plot(time(range),plotsign(2,range),'gr','LineWidth',6);hold on
    plot(time(range),plotsignindivLess(range),'','LineWidth',6,'Color',darkred);hold on
    plot(time(range),plotsignindivMore(range),'','LineWidth',6,'Color',darkblue);hold off
    xlabel('time in sec','FontSize',20);ylabel('Effect size (mean +/- SE)','FontSize',20);
    title([glmname ' ' con_name{contrIND} ' MoreDay(blue)&LessDay(red)'],'FontSize',15);
    set(gca,'FontSize',20);axis([min(time(range)) 2.5 min([mean1more(:);mean1less(:)])-.5  max([mean1more(:);mean1less(:)])+.1]);
    h = gcf;
    %     saveas(h,[pathsave glmname{1} ' ' con_name{contrIND} 'permsign'],'fig');
    %     saveas(h,[pathsave glmname{1} ' ' con_name{contrIND} 'permsign'],'eps');
    
    cd(pathsave)
    
    saveas(gcf,['GLM_memtest_randpermute_' con_name{contrIND} '_delayed_MoreDay&LessDay.fig'])
    print(['GLM_memtest_randpermute_' con_name{contrIND} '_delayed_MoreDay&LessDay'],'-dpdf','-fillpage') % save as pdf
    print(['GLM_memtest_randpermute_' con_name{contrIND} '_delayed_MoreDay&LessDay'],'-depsc','-tiff','-r0') % save as eps
    
    clear dat1 cg vg tg c
    
    store_plotsign{contrIND} = plotsign;clear plotsign
end

%% assess peak in different contrasts
%
% tw = [1000:2700];
%
% maxval = []; maxtime = [];
% meanNode = 20;
%
% % extract peak values for each ID and each contrast
% for con_counter = 1:numContr
%     for id = 1:length(IDs)
%         [maxval(id,con_counter) maxtime(id,con_counter)] = max(storedat(id).glm.con{con_counter}(tw));
%     end
% end


%% feedback bits : only for the pilot 1

% %% preparation
%
%
% clear;clc;close all;
% rng(0)
% addpath(genpath('/Users/yeojin/Desktop/B_scripts/functions'))
%
%
% % paths
% path_par    = '/Users/yeojin/Desktop/';
% path        = [ path_par 'E_data/EA_raw/' ]; % maybe the mat file with the behavioural results
% pathdata    = [ path_par 'E_data/EB_cleaned/' ];
% pathsave    = [ path_par 'C_writings/CB_figures/GLMs/' ]; % where do you save your figures?
%
% % variables
% rewcond = 2; nrewcond = 1; maxtr = 120;
% IDs = [1006 1007 1008 1011 1012 1013 1014 1015 1016 1017]; % this variable used to be rangeids
%
% % load behavioural results for each participants
% fname_beh =[]; expdat = [];
% for ids = 1:length(IDs)
%     fname_beh{ids,1} = [ num2str(IDs(ids)) '.mat' ];
%     expdat{ids,1}    = load([path fname_beh{ids,1}]);
% end
%
% % stat presets
% randpermute = 1;  % do the random permutation? (1 = yes! / 0 = nope)
% length_randpermute_indiv = 100; length_randpermute = 1000;
% sign_randpermute = .05; % less than 5 percent of randpermute results have to be learger than age effect
% signcut = 10; % significance level for group contrast
% tw = [500:2500];
%
%
% %% run: main task, INDIVIDUAL level
%
% storedat = [];
%
% for id = 1:length(IDs)
%
%     load([pathdata 're1' num2str(IDs(id)) '_eye_rv200_fb.mat'],'eye_rv');
%
%     % make missing trials NaNs
%     alltrials   = eye_rv.cfg.previous.previous.trl(:,1);
%     newtrials   = eye_rv.sampleinfo(:,1);
%     trlnum      = find(ismember(alltrials,newtrials));
%
%     %bsl correct (200 to 0 ms time window)
%     for trl = 1:length(eye_rv.trial)
%         eye_rv.trial{trl} = eye_rv.trial{trl}-nanmean(eye_rv.trial{trl}(1:200));
%     end
%
%
%     % indexing trial properties (for me)
%     indRew      = find(cell2mat(eye_rv.trlord(:,2)) == rewcond);
%     indnRew     = find(cell2mat(eye_rv.trlord(:,2)) == nrewcond);
%     indRew_mH   = find(cell2mat(eye_rv.trlord(:,2)) == rewcond & cell2mat(eye_rv.trlord(:,3))>=0.85);
%     indRew_mL   = find(cell2mat(eye_rv.trlord(:,2)) == rewcond & cell2mat(eye_rv.trlord(:,3))<=0.68);
%     indnRew_mH  = find(cell2mat(eye_rv.trlord(:,2)) == nrewcond & cell2mat(eye_rv.trlord(:,3))>=0.85);
%     indnRew_mL  = find(cell2mat(eye_rv.trlord(:,2)) == nrewcond & cell2mat(eye_rv.trlord(:,3))<=0.68);
%     indmH       = find(cell2mat(eye_rv.trlord(:,3))>=0.85);
%     indmL       = find(cell2mat(eye_rv.trlord(:,3))<=0.68);
%
%     % assign dummy variables for GLM analysis
%
% %     glmname     = {'REWARD versus NoREWARD'};
% %     numContr    = 1;
% %
% %     dummie      = zeros(length(trlnum),1);
% %     dummie(indRew) = 1;
% %     beh.Reward  = zscore(dummie); clear dummie;
% %
% %     storebehav  = beh;
% %     storepupil  = eye_rv;
%
% %     memorabilities = cell2mat(expdat{id,1}.dat.SESSION.Session_1.memtest.results.trlord(:,3));
% %     dummie = memorabilities(trlnum);
% %     beh.Memorability = dummie; clear dummie;
% %
% %     beh.RewMemorability = zscore(beh.Reward.*beh.Memorability);
%
%     glmname     = {'HiMemVsLoMem'};
%     numContr    = 1;
%
%     dummie      = zeros(length(trlnum),1);
%     dummie(indmH) = 1;
%     beh.HiMem  = zscore(dummie); clear dummie;
%
%     storebehav  = beh;
%     storepupil  = eye_rv;
%
%
% %
% %     dummie      = zeros(length(trlnum),1);
% %     dummie(indnRew) = 1;
% %     beh.NoReward= zscore(dummie); clear dummie;
% %
% %     storebehav  = beh;
% %     storepupil  = eye_rv;
% %
% %     glmname     = {'Reward and Memorability:'};
% %     numContr    = 4; % 1: Reward_hiMem 2: Reward_loMem 3: NoReward_hiMem 4: NoReward_lomem
% %
% %     dummie      = zeros(length(trlnum),1);
% %     dummie(indRew_mH) = 1;
% %     beh.RewardHiMem  = zscore(dummie); clear dummie;
% %
% %     dummie      = zeros(length(trlnum),1);
% %     dummie(indRew_mL) = 1;
% %     beh.RewardLoMem= zscore(dummie); clear dummie;
% %
% %     dummie      = zeros(length(trlnum),1);
% %     dummie(indnRew_mH) = 1;
% %     beh.NoRewardHiMem  = zscore(dummie); clear dummie;
% %
% %     dummie      = zeros(length(trlnum),1);
% %     dummie(indnRew_mL) = 1;
% %     beh.NoRewardLoMem= zscore(dummie); clear dummie;
%
%     storebehav  = beh;
%     storepupil  = eye_rv;
%
%
%     % build structures for eye_glm (beh is already done above)
%
%     pupil       = eye_rv;
%     pupil.trial = []; pupil.trial = storepupil.trial;
%     pupil.time  = []; pupil.time = storepupil.time;
%     pupil.sampleinfo = []; pupil.sampleinfo = storepupil.sampleinfo;
%
%     reg_names   = {'HiMem'}; % regressor names
% %     reg_names   = {'RewardHiMem'; 'RewardLoMem'; 'NoRewardHiMem'; 'NoRewardLoMem';};
% %     contrast    = [1 -1 0; 1 0 0; 0 1 0; 0 0 1]; d
%     con_name    = {'HiMem pos'};
%     contrast    = 1;
% %     con_name    = {'RewardHiMem'; 'RewardLoMem'; 'NoRewardHiMem'; 'NoRewardLoMem';};
%
%     options.rejBad = 0;
%     options.bl_corr = 0;
%     options.bl_win = [-.2 0];
%
%     clear eye_rv
%
%
%     % run GLM
%     storedat(id).glm = eye_glm(pupil,beh,reg_names,contrast,con_name,options)
%
%
%     % make randpermute for each person
%     if randpermute
%
%         for permcount = 1:length_randpermute_indiv
%
%             eval(['lengthbehav = length(beh.' reg_names{1} ');'])
%             ind = randperm(lengthbehav);
%
%             for c1 = 1:size(reg_names,1)
%                eval(['beh.' reg_names{c1} '= beh.' reg_names{c1} '(ind);'])
%             end
%
%             dummie = eye_glm(pupil,beh,reg_names,contrast,con_name,options);
%             for stats = 1:size(reg_names,1)+1
%                 contrasts(id,stats,permcount,:) = dummie.con{1,stats};
%                 residuals(id,stats,permcount,:) = dummie.res{1,stats};
%                 tstats(id,stats,permcount,:) = dummie.tstat{1,stats};
%             end
%
%             clear ind dummie
%
%         end % close random permutation loop
%
%     end % close random permutation calculation
%
%     fprintf(['\n*********\n' num2str(IDs(id)) ' done\n*********\n'])
%
% end
%
% %% run: main task, stimulus bits, GROUP level
%
% con_name{length(con_name)+1} = 'intercept'; % last contrast is intercept
%
% dm = ones(length(IDs),1); % design matrix for groups (which will come later)
%
% nLC = size(storedat(1).glm.con,2); % number of lower level contrasts (each experimental condition)
% nSamp = size(storedat(1).glm.con{1},2); % number of samples in timeseries
%
% % cg = nan(nHC,nSamp); vg = nan(nHC,nSamp); tg = nan(nHC,nSamp);
% % cgb = nan(length(con_name),length_randpermute,nHC,nSamp); vgb = nan(length(con_name),length_randpermute,nHC,nSamp); tgb = nan(length(con_name),length_randpermute,nHC,nSamp);
%
%
% for con_counter = 1:length(con_name)
%
%     for i = 1:length(IDs)
%         % save contrast estimates for regressors
%         dat1(i,:) = storedat(i).glm.con{1,con_counter};
%     end
%     con_estimates = dat1;
%
%
%     % contrast data from each age group for plotting (not the randpermute data)
%     time    = pupil.time{1};
%     range   = 1:length(time);
%     means   = nanmean(dat1(:,:));
%     SEs     = nanstd(dat1(:,:))./sqrt(length(IDs));
%
%     % 5%, expecting, higher than 0 / assess maximal and minimal value in particular time window
%     datdum = squeeze(nanmean(contrasts(:,con_counter,:,:),1));
%     for rep = 1:size(contrasts,3)
%         datmin(rep) = min(squeeze(datdum(rep,min(tw):max(tw))));
%         datmax(rep) = max(squeeze(datdum(rep,min(tw):max(tw))));
%     end
%     plotsign_indiv = nan(1,size(con_estimates,2));
%     plotsign_indiv(find(means(min(tw):max(tw)) < prctile(datmin,5))+min(tw))  = 0;
%     plotsign_indiv(find(means(min(tw):max(tw)) > prctile(datmax,95))+min(tw)) = 0;
%
%     clear datmin datmax rep datdum
%
%
%     %  plot with significance
%     darkred = [.9 , 0.1 , 0.1]; darkblue = [.1 , 0.1 , .9];
%     figure;
%     shadedErrorBar(time(range),means(range),SEs(range),'lineprops',{'k-o','markerfacecolor','k'}); hold on
%     plot(time(range),plotsign_indiv(range),'','LineWidth',6,'Color',darkblue);hold off
%     xlabel('time(sec)','FontSize',20);ylabel('Effect size (mean +/- SE)','FontSize',20);
%     title([glmname{1} ' ' con_name{con_counter} ' bstrsign'],'FontSize',15);
%     set(gca,'FontSize',20);
%     axis([min(time(range)) 2.5 -.3 .3]);
% %     axis([min(time(range)) 2.5 min(means(:))-.3  max(means(:))+.3]);
%
%     % save the figures
%     cd(pathsave)
%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
%     print(['GLM_memtest_randpermute_' con_name{con_counter} '_fb_pilot1'],'-dpdf','-fillpage') % save as pdf
%     print(['GLM_memtest_randpermute_' con_name{con_counter} '_fb_pilot1'],'-depsc','-tiff','-r0') % save as eps
%
%     clear dat1 con_estimates fig
%
%
% end % close contrast calculation loop


%% plot beta signs



