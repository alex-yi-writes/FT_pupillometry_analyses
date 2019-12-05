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
rng(0)
addpath(genpath('/Users/yeojin/Desktop/B_scripts/BC_functions'))

pilotNr     = 1;

% paths
path_gen    = '/Users/yeojin/Desktop/';
path_par    = [ path_gen 'E_data/' ];
path        = [ path_par 'EA_raw/EAC_behav/pilot_3T/' ]; % maybe the mat file with the behavioural results
pathdata    = [ path_par 'EB_cleaned/EBB_pupil/pilot_3T/' ];
pathsave    = [ path_gen 'C_writings/CB_figures/pilot_3T/GLMs/' ]; % where do you save your figures?

% variables
IDs = [1204 1107 1109 1210 1111 1212 1113 1214 1216]; % this variable used to be rangeids
% IDs = [1107 1212 1216]
%maxtr = 120;
oddcond = 1; stdcond = 2;

% load behavioural results for each participants
fname_beh =[]; expdat = [];
for ids = 1:length(IDs)
    fname_beh{ids,1} = [ num2str(IDs(ids)) '_OB.mat' ];
    expdat{ids,1}    = load([path fname_beh{ids,1}]);
end
load('/Users/yeojin/Desktop/B_scripts/BB_analyses/BBB_eyetracking/nulltrials_OB.mat')

% stat presets
randpermute = 1;  % do the random permutation? (1 = yes! / 0 = nope)
length_randpermute_indiv = 300; length_randpermute = 1000;
sign_randpermute = .05; % less than 5 percent of randpermute results have to be learger than age effect
signcut = 10; % significance level for group contrast
tw = [201:2701];

% call in indice
% indice = load([pathdata 'pilot' num2str(pilotNr) '_index.mat']);


%% run: main task, INDIVIDUAL level

storedat = [];

for id = 1:length(IDs)
    
    load([pathdata 'ob' num2str(IDs(id)) '_eye_rv200_stim.mat'],'eye_rv');
    
    % make missing trials NaNs
    alltrials   = eye_rv.cfg.previous.previous.trl(:,1);
    newtrials   = eye_rv.sampleinfo(:,1);
    trlnum      = find(ismember(alltrials,newtrials));
    
    %bsl correct (200 to 0 ms time window)
    for trl = 1:length(eye_rv.trial)
        eye_rv.trial{trl} = eye_rv.trial{trl}-nanmean(eye_rv.trial{trl}(1:200));
    end
    
    % clean the null trials and rejected trials: breaks
    trl_raus  = eye_rv.trl_raus;
    null_trls = nulltrials(:,expdat{id,1}.dat.oddball.StimSchedule);
    rmrow_stim   = find(trl_raus==1); % rejected trials from preprocessing
    
    % stimulus informations: filename, reward/no-reward, memorability
    stimdat      = expdat{id,1}.dat.oddball.results.trl;
    stimdat(isemptycell(stimdat)) = {NaN}; % this line, because if you convert data directly from cells to double, it automatically removes empty cells, which in our case are breaks.
    indices = cell2mat(stimdat(:,2)); 
    
    % remove breaks and rejected trials
    stimdat(null_trls,:) = [];
    stimdat(rmrow_stim,:) = []; % and now, we remove the rows of trials we rejected during data cleaning
    
    % indexing trial properties (for me)
    indOdd      = find(cell2mat(stimdat(:,2)) == oddcond);
    indStd     = find(cell2mat(stimdat(:,2)) == stdcond);
    
    % assign dummy variables for GLM analysis
    
    glmname     = {'Oddball versus Standard'};
    numContr    = 1;
    
    dummie      = zeros(length(trlnum),1);
    dummie(indOdd) = 1;
    beh.Oddball  = zscore(dummie); clear dummie; %%%%%
    
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
    %     dummie(indStd_mH) = 1;
    %     beh.NoRewardHiMem  = zscore(dummie); clear dummie;
    %
    %     dummie      = zeros(length(trlnum),1);
    %     dummie(indStd_mL) = 1;
    %     beh.NoRewardLoMem= zscore(dummie); clear dummie;
    %
    %     storebehav  = beh;
    %     storepupil  = eye_rv;
    
    
    % build structures for eye_glm (beh is already done above)
    
    pupil       = eye_rv;
    pupil.trial = []; pupil.trial = storepupil.trial;
    pupil.time  = []; pupil.time = storepupil.time;
    pupil.sampleinfo = []; pupil.sampleinfo = storepupil.sampleinfo;
    
    reg_names   = {'Oddball'}; % regressor names
    %     reg_names   = {'RewardHiMem'; 'RewardLoMem'; 'NoRewardHiMem'; 'NoRewardLoMem';};
    %     contrast    = [1 -1 0; 1 0 0; 0 1 0; 0 0 1]; d
    con_name    = {'Oddball pos'};
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
    storedat(id).glm = eye_glm(pupil,beh,reg_names,contrast,con_name,options)
    
    
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
                contrasts(id,stats,permcount,:) = dummie.con{1,stats};
                residuals(id,stats,permcount,:) = dummie.res{1,stats};
                tstats(id,stats,permcount,:) = dummie.tstat{1,stats};
            end
            
            clear ind dummie
            
        end % close random permutation loop
        
    end % close random permutation calculation
    
    fprintf(['\n*********\n' num2str(IDs(id)) ' done\n*********\n'])
    
end

%% run: main task, stimulus bits, GROUP level

con_name{length(con_name)+1} = 'intercept'; % last contrast is intercept

dm = ones(length(IDs),1); % design matrix for groups (which will come later)

nLC = size(storedat(1).glm.con,2); % number of lower level contrasts (each experimental condition)
nSamp = size(storedat(1).glm.con{1},2); % number of samples in timeseries

% cg = nan(nHC,nSamp); vg = nan(nHC,nSamp); tg = nan(nHC,nSamp);
% cgb = nan(length(con_name),length_randpermute,nHC,nSamp); vgb = nan(length(con_name),length_randpermute,nHC,nSamp); tgb = nan(length(con_name),length_randpermute,nHC,nSamp);


for con_counter = 1:length(con_name)
    
    for i = 1:length(IDs)
        % save contrast estimates for regressors
        dat1(i,:) = storedat(i).glm.con{1,con_counter};
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
    axis([min(time(range)) 2.5 -.3 .3]);
    %     axis([min(time(range)) 2.5 min(means(:))-.3  max(means(:))+.3]);
    
    % save the figures
    cd(pathsave)
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    
    print(['GLM_maintask_randpermute_' con_name{con_counter} '_stim_pilot_3T_oddball'],'-dpdf','-fillpage') % save as pdf
    print(['GLM_maintask_randpermute_' con_name{con_counter} '_stim_pilot_3T_oddball'],'-depsc','-tiff','-r0') % save as eps
    
    
    clear dat1 con_estimates fig
    
    
end % close contrast calculation loop

clear id con_counter

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
%     indStd     = find(cell2mat(eye_rv.trlord(:,2)) == nrewcond);
%     indRew_mH   = find(cell2mat(eye_rv.trlord(:,2)) == rewcond & cell2mat(eye_rv.trlord(:,3))>=0.85);
%     indRew_mL   = find(cell2mat(eye_rv.trlord(:,2)) == rewcond & cell2mat(eye_rv.trlord(:,3))<=0.68);
%     indStd_mH  = find(cell2mat(eye_rv.trlord(:,2)) == nrewcond & cell2mat(eye_rv.trlord(:,3))>=0.85);
%     indStd_mL  = find(cell2mat(eye_rv.trlord(:,2)) == nrewcond & cell2mat(eye_rv.trlord(:,3))<=0.68);
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
% %     memorabilities = cell2mat(expdat{id,1}.dat.SESSION.Session_1.maintask.results.trlord(:,3));
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
% %     dummie(indStd) = 1;
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
% %     dummie(indStd_mH) = 1;
% %     beh.NoRewardHiMem  = zscore(dummie); clear dummie;
% %
% %     dummie      = zeros(length(trlnum),1);
% %     dummie(indStd_mL) = 1;
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
%     print(['GLM_maintask_randpermute_' con_name{con_counter} '_fb_pilot1'],'-dpdf','-fillpage') % save as pdf
%     print(['GLM_maintask_randpermute_' con_name{con_counter} '_fb_pilot1'],'-depsc','-tiff','-r0') % save as eps
%
%     clear dat1 con_estimates fig
%
%
% end % close contrast calculation loop


%% plot beta signs



