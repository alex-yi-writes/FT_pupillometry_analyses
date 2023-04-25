% clean eyetracker data

clear all
close all
clc

% addpath('/Users/dorothea/Dropbox/matlab/toolboxes/fieldtrip-20170119'); % add your path to fieldtrip folderpathdat = ['\\imed\fme\kliniken\iknd\rakshit\Eigene Dateien\Dorothea\Pupil_clining\EM_raw\'];
% addpath('/Users/dorothea/Dropbox/work/Lehre_Innsbruck/Veranstaltungen_SS22/Sem_Forschung_II/FACES/Eyetracker_Skripte/');
%pathdat = 'C:\Users\FriedaLubre\Documents\DH_Innsbruck\Faces\Test Final\';
pathdat = ['/Users/yeojin/Downloads/immediate_segs3/'];
savepathdat = ['/Users/yeojin/Downloads/immediate_segs3/'];
%addpath ('C:\Users\FriedaLubre\Documents\MATLAB\spm12\external\fieldtrip\')

perc_cutoff = 20; % if more than perc_cutoff% interpolated per trial, remove whole trial
rejectlatency = [-.2 3];

all_files = dir (strcat(pathdat,'*mem_eyedat.mat')); % for stim_eyedat
%all_files = dir (strcat(pathdat,'*fb_eyedat.mat')); % for stim_eyedat
%ids = ["rd12110mem_eyedat.mat","rd12112mem_eyedat.mat","rd12114mem_eyedat.mat","rd12115mem_eyedat.mat"]; % list of ids


for n = 1:length(all_files)
    ids = all_files(n).name;
    data = []; trls = []; samples = [];
    
    for p = 2 % concatenate across two runs
        
       % load([pathdat num2str(ids(n)) '/face' num2str(ids(n)) num2str(p) '_eyedat.mat']);
        %not using the below code as the file naming convention is
        %different
        %load([pathdat 'e' num2str(ids(n)) num2str(p) '_eyedat.mat']);
        load([strcat(pathdat,ids)])
        %eye_mem = eye_fb;%when it is eye_fb.mat
        datatemp = cell2mat(eye_mem.trial');
        samples = [samples;eye_mem.sampleinfo+(p-1)*1*10^8];%%% newly added, should be done with all future datasets
        trls = [trls;eye_mem.cfg.trl+(p-1)*1*10^8];%%% newly added, should be done with all future datasets
        data = [data; datatemp];clear datatemp
    end
    
    % cleaning according to Mathot's paper
    dat = f_Clean_Mathots_way_v2(data,perc_cutoff);
    
    datai = dat.datai; % interpolated automatically cleaned data
    datac = dat.removed; % documentation of removed artefacts per trial
    
    % trl_raus1 == autom. approach
    trl_raus = dat.trl_raus; %% 1=rejected trials, 0=selected trials
    auto_select = dat.auto_select; %% Index of selected trails from automatic correction
    
    %% %%%%%%%%%%% can be commented out, is just showing old (blue) and red (interpolated) data
     %  figure
   
    % datai == selected data !!!
    
%     for i = 1:size(datai,1)
%         subplot(15,15,i);plot(data(i,:));hold on
%         plot(datai(i,:),'r');hold off
%     end

%%%%%%%plot Start
    for i = 1:size(auto_select,1)
        subplot(10,9,i);plot(data(auto_select(i),:));hold on
        plot(datai(i,:),'r');hold off
        if i == 90
            break;
        end
    end
%%%%%%%plot end
    h = gcf;set(h,'Position',[100 400 1300 800]);
    %%%%%%%%%%%%% can be commented out
    %% make new data structure from interpolated data for manual rejection
    eyeinterp = eye_mem;
    
    %%Deleting rejected trails' info
    
    eyeinterp.trial(find(trl_raus==1)) = [];
    eyeinterp.time(find(trl_raus==1)) = [];
    eyeinterp.trl_raus = trl_raus;
    eyeinterp.dat = dat;
    eyeinterp.data_before_corr = data;

    %%Save the corrected data
    %mkdir([pathdat num2str(ids(n)) '/cleaned']); % make directory to save cleaned data
    save([strcat(savepathdat,replace(ids,".mat","_cln.mat"))],'eyeinterp'); % save cleaned data and behavioral data
    clear eye_mem
    clear eyeinterp



    % make data structure new, and baseline correct before rejection trials
    % individually
%     for tr = 1: size(datai,1)
%         temp = zscore(datai(tr,:)); % z-scoring
%         eyeinterp.trial{1,tr} = temp;clear temp
%         eyeinterp.time{1,tr} = eye.time{1};
%     end
%     for tr = 1: size(auto_select,1)
%         temp = zscore(datai(tr,:)); % z-scoring
%         eyeinterp.trial{1,tr} = temp ; clear temp %atai(tr,:);%clear temp
%         eyeinterp.time{1,tr} = eye.time{1};
%     end
%     eyeinterp.cfg.trl = trls;
%     eyeinterp.sampleinfo = samples;
%     
%     eyeinterp.cfg.trl(find(trl_raus == 1),:)=[];
%     eyeinterp.sampleinfo(find(trl_raus == 1),:)=[];
%     
%     clear eye
%     
%     %% reject individual trials manually
%         
%     cfg.method = 'trial'; % then individual trials excluded manually
%     cfg.alim = 4; % scaling in z range
%     cfg.latency = rejectlatency;
%     eye_rv = ft_rejectvisual(cfg, eyeinterp)
%     
%     % determine which trials out based on manual rejection
%     trlnew = eye_rv.sampleinfo(:,1);
%     trlold = eye_rv.cfg.previous.trl(:,1);
%     
%     % trl_raus == add manual 
%     trl_raus_manu = ~ismember(trlold,trlnew);%  - index based on whole dataset
%     
%     ['trials removed = ' num2str(sum(trl_raus))]
%     
%     trl_raus = [trl_raus trl_raus]; %1st col = automatic ; 2nd col = Manual
%     trl_raus(auto_select(find(trl_raus_manu==1)),2) = 1;
%     %% load behavioral dataset, and add behavioral data to cleaned pupil data (keep behavioral data in for deleted trials, but add which trials have been excluded)
%     pics = []; rts = []; resp = [];
%     for p = 1:2
%         load([pathdat num2str(ids(n)) '/faces_' num2str(ids(n)) '_' num2str(p) '_behav.mat']);
%         pics = [pics exp.dat.picnums];
%         rts = [rts exp.dat.RT];
%         resp = [resp exp.dat.responsequest];
%     end
%     
%     eye_rv.sex = exp.dat.sex;
%     eye_rv.age = exp.dat.age;
%     
%     eye_rv.trlraus = trl_raus; % 2 columns: 1 = auto, 2 = auto+manual; 1 if trial excluded, 0 if in
%     eye_rv.trlnum = 1:length(pics); % just consecutive trial number
%     eye_rv.rts = rts; % clear rts
%     eye_rv.resp = resp; % clear resp
%     eye_rv.pics = pics; % clear resp
%     eye_rv.intpupil_final = datai(find(trl_raus_manu==0),:);
%     
%     
%     mkdir([pathdat num2str(ids(n)) '/cleaned']); % make directory to save cleaned data
%     save([pathdat  num2str(ids(n)) '/cleaned/' num2str(ids(n)) '_eye_cleaned.mat'],'eye_rv'); % save cleaned data and behavioral data
%     
%     keep pathdat n ids rejectlatency
%     close all
end