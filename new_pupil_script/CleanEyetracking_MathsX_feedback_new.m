% clean eyetracker data

clear all
close all
clc

% addpath('/Users/dorothea/Dropbox/matlab/toolboxes/fieldtrip-20170119'); % add your path to fieldtrip folderpathdat = ['\\imed\fme\kliniken\iknd\rakshit\Eigene Dateien\Dorothea\Pupil_clining\EM_raw\'];
addpath('/Users/dorothea/Dropbox/work/Lehre_Innsbruck/Veranstaltungen_SS22/Sem_Forschung_II/MathsX/Eyetracker_Skripte/');
pathdat = '/Users/dorothea/Dropbox/work/Lehre_Innsbruck/Veranstaltungen_SS22/Sem_Forschung_II/MathsX/Daten/';


perc_cutoff = 33; % if more than perc_cutoff% interpolated per trial, remove whole trial
rejectlatency = [-.2 3];

ids = [25:34];

for n = 1:length(ids)
    data = []; trls = []; samples = []; all_trial = [];
    
    for p = 1:2 % concatenate across two runs
        
        load([pathdat num2str(ids(n)) '/math' num2str(ids(n)) num2str(p) '_eyedat_feedback.mat']);
        
        datatemp = cell2mat(eye.trial');
        samples = [samples;eye.sampleinfo+(p-1)*1*10^8];%%% newly added, should be done with all future datasets
        trls = [trls;eye.cfg.trl+(p-1)*1*10^8];%%% newly added, should be done with all future datasets
        data = [data; datatemp];clear datatemp
        all_trial = [all_trial;eye.trial];
    end
    
    % cleaning according to Mathot's paper
    dat = f_Clean_Mathots_way_v2(data,perc_cutoff);
    
    datai = dat.datai; % interpolated automatically cleaned data
    datac = dat.removed; % documentation of removed artefacts per trial
    
    % trl_raus1 == autom. approach
    trl_raus = dat.trl_raus; %% 1=rejected trials, 0=selected trials
    auto_select = dat.auto_select; %% Index of selected trails from automatic correction
    
    
    for i = 1:size(auto_select,1)
        subplot(12,12,i);plot(data(auto_select(i),:));hold on
        plot(datai(i,:),'r');hold off
    end
    h = gcf;set(h,'Position',[100 400 1300 800]);
    %%%%%%%%%%%%% can be commented out
    %% make new data structure from interpolated data for manual rejection
    eyeinterp = eye;
    
    eyeinterp.cfg.trl = trls;
    eyeinterp.sampleinfo = samples;
    
    %%Deleting rejected trails' info
    
    eyeinterp.trial = all_trial;
    eyeinterp.trial(find(trl_raus==1)) = [];
    %eyeinterp.time(find(trl_raus==1)) = [];
    eyeinterp.time = [];
    
    for tr = 1: size(auto_select,1)
        temp = zscore(datai(tr,:)); % z-scoring
        eyeinterp.trial{1,tr} = temp ; clear temp %atai(tr,:);%clear temp
        eyeinterp.time{1,tr} = eye.time{1};
    end
    
    eyeinterp.cfg.trl(find(trl_raus == 1),:)=[];
    eyeinterp.sampleinfo(find(trl_raus == 1),:)=[];
    
    clear eye
    
    %% reject individual trials manually
    
    cfg.method = 'trial'; % then individual trials excluded manually
    cfg.alim = 4; % scaling in z range
    cfg.latency = rejectlatency;
    eye_rv = ft_rejectvisual(cfg, eyeinterp)
    
    % determine which trials out based on manual rejection
    trlnew = eye_rv.sampleinfo(:,1);
    trlold = eye_rv.cfg.previous.trl(:,1);
    
    % trl_raus == add manual
    trl_raus_manu = ~ismember(trlold,trlnew);%  - index based on whole dataset
    
    ['trials removed = ' num2str(sum(trl_raus))]
    
    trl_raus = [trl_raus trl_raus]; %1st col = automatic ; 2nd col = Manual
    trl_raus(auto_select(find(trl_raus_manu==1)),2) = 1;
    
    %% load behavioral dataset, and add behavioral data to cleaned pupil data (keep behavioral data in for deleted trials, but add which trials have been excluded)
    pics = []; rts_task = []; resp_task = []; corrrespside = []; corrresp = []; cue = []; rts_diffi = []; resp_diffi = [];
    for p = 1:2
        load([pathdat num2str(ids(n)) '/mathsx_' num2str(ids(n)) '_' num2str(p) '_behav.mat']);
        exp.dat.response(end+1:length(exp.dat.picnums)) = 0; % if missings at the end, fill up with zeros, same for corrresp
        exp.dat.response_corr(end+1:length(exp.dat.picnums)) = 0; % if missings at the end, fill up with zeros, same for corrresp
        cueDum = zeros(1,length(exp.dat.picnums));
        cueDum(find(exp.dat.picnums <33)) = 1;
        cueDum(find(exp.dat.picnums > 96)) = 4;
        cueDum(find(exp.dat.picnums <65 & exp.dat.picnums > 32)) = 2;
        cueDum(find(exp.dat.picnums <96 & exp.dat.picnums > 64)) = 3;
        pics = [pics exp.dat.picnums];
        rts_task = [rts_task exp.dat.RTtask];
        resp_task = [resp_task exp.dat.response];
        corrrespside = [corrrespside exp.dat.correct_response_side];
        corrresp = [corrresp exp.dat.response_corr];
        cue = [cue cueDum]; % cuename was unreliable, so sorted after picnums
        rts_diffi = [rts_diffi exp.dat.RTdiffi];
        resp_diffi = [resp_diffi exp.dat.responsediffi];
    end
    
    eye_rv.sex = exp.dat.sex;
    eye_rv.age = exp.dat.age;
    
    eye_rv.trlraus = trl_raus(:,2); % 2 columns: 1 = auto, 2 = auto+manual; 1 if trial excluded, 0 if in
    eye_rv.trlnum = 1:length(pics); % just consecutive trial number
    eye_rv.rts_task = rts_task;
    eye_rv.resp_task = resp_task;
    eye_rv.pics = pics;
    eye_rv.rts_diffi = rts_diffi;
    eye_rv.resp_diffi = resp_diffi;
    eye_rv.cue = cue;
    eye_rv.corrrespside = corrrespside;
    eye_rv.corrresp = corrresp;
    eye_rv.solutionvec = exp.dat.solutionvec;
    eye_rv.taskvec = exp.dat.taskvec;
    
    mkdir([pathdat num2str(ids(n)) '/cleaned']); % make directory to save cleaned data
    save([pathdat  num2str(ids(n)) '/cleaned/' num2str(ids(n)) '_eye_cleaned_feedback.mat'],'eye_rv'); % save cleaned data and behavioral data
    
    keep pathdat n ids rejectlatency perc_cutoff
    close all
end