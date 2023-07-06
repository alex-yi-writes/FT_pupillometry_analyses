%% Pupillometric data preprocessing: encoding phase
%   main task eyedata
 
%% work log
 
%   27_11_2018 created the script
%              everything works fine
%   13_06_2019 added two session system now
 
%% preparation
 
clear;clc;close all;
warning off
 
% path setup
addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20220104')); % change those two
addpath(genpath('/Users/yeojin/Desktop/B_scripts/BA_preprocessing/BAA_pupillometry/preproc_functions/'))
 
% which data will you be processing? (for now, it's all 1)
answ        = inputdlg({'which session? 1 or 2'});
days        = str2num(answ{1,1}); clear answ
% answ        = inputdlg({'which pilot?'});
% pilotNr     = str2num(answ{1,1}); clear answ
pilotNr     = 2;
 
path_parent = '/Users/yeojin/Desktop/E_data/EA_raw/EAB_pupil/MRPETpilot/ASC/maintask/';
path_dat    = [ path_parent ]; % tailor it to your original path
path_behav  = '/Users/yeojin/Desktop/E_data/EA_raw/EAC_behav/MRPETpilot/';
path_save   = '/Users/yeojin/Desktop/E_data/EB_cleaned/EBB_pupil/MRPETpilot/segmented/task_segs/';
 
% variables
if days == 1
    VPs = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 ...
     2127 2128 2129 2130 2131 2132 2233 2234 2235 2236 2237 2238 2239 2240 2142 2143 2144 2145 2146 2147 2148 2249 2250];
else
    VPs = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 ...
     2127 2128 2129 2130 2131 2132 2233 2234 2235 2236 2237 2238 2239 2240 2142 2143 2144 2145 2146 2147 2148 2249 2250]; % 2203: empty
end
 
fname_raweye = []; fname_beh =[]; expdat = [];
for i1 = 1:length(VPs)
    fname_raweye{i1,1}  = [ 're' num2str(days) num2str(VPs(i1)) '.asc' ]; % rewardtask-encoding (main task)
    fname_raweye{i1,2}  = [ 'ri' num2str(days) num2str(VPs(i1)) '.asc' ]; % rewardtask-immediate recall
    fname_raweye{i1,3}  = [ 'rd' num2str(days) num2str(VPs(i1)) '.asc' ]; % rewardtask-delayed recall
    fname_beh{i1,1}     = [ num2str(VPs(i1)) '_' num2str(days) '.mat' ];
    expdat{i1,1} = load([path_behav fname_beh{i1,1}]);
end
 
% script specific variables 
pathdat       = fullfile(path_dat);
pathdatbehav  = fullfile(path_behav,fname_beh);
 
fprintf('\n Preparation done \n')
%% import eyetracker data
 
close all
 
for vps = 38:length(VPs)
    
    fprintf('\n importing and segmenting ID: %d\n', VPs(vps))
    
    name = fname_raweye{vps,1};
        
    %%%%%%%%%%%%% alex added this bit to solve the problem where this function
    %%%%%%%%%%%%% could not read in the calibration bit altogether
    
    % file setup : you may need to tailor here to your config accordingly
    pathfile = pathdat;
    filename = name;
    
    cd(fullfile(pathdat))
    copyfile(fullfile(pathdat,filename), ['backup_' filename]); % back up the original file
    
    fid =fopen(filename); % read in the raw data with calibration info
    C=textscan(fid,'%s','delimiter','\n');
    fclose(fid);

    try
        i1 = 0; i3 = 0; clear startline endline
        for k=1:numel(C{1,1})
            tmp1 = regexp(C{1,1}(k),'CAMERA_CONFIG'); % find where header ends
            if ~isemptycell(tmp1)
                i1 = i1+1;
                startline(i1) = k+2; % mark the start of lines deleted
            end
            tmp2 = regexp(C{1,1}(k),'RECCFG CR 1000 2 1 R'); % find where the recording starts
            if ~isemptycell(tmp2)
                i3 = i3+1;
                endline(i3) = k-1; % mark the end of lines deleted
            else
                tmp2 = regexp(C{1,1}(k),'RECCFG CR 1000 2 1 L');
                if ~isemptycell(tmp2)
                    i3 = i3+1;
                    endline(i3) = k-1; % mark the end of lines deleted
                end
            end

        end
        fprintf('\n importing done \n')

        newtext = [];
        if ~isempty(startline) && ~isempty(endline)

            newtext{1,1} = [C{1,1}(1:startline-1);C{1,1}(endline+1:end)];

            fprintf('\n file modified \n')

            % print new file
            fName = name;
            fid = fopen(fName,'w'); % open the file
            for k=1:numel(newtext{1,1})
                fprintf(fid,'%s\r\n',newtext{1,1}{k,1});
            end
            fclose(fid);
            fprintf('\n saved \n')
            clear C newtext
        else
            clear C newtext
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch
    end
    
    eyeAdjustTrigNam([pathdat '/' name],'MSG','INPUT');
    [type] = ft_filetype([pathdat '/' name])
    
    % stimulus segment
    cfg                     = [];
    cfg.trialfun            = 'my_trialfunction_eyetracker_rewardtask_2sess_stim'; % for the second pilot, which didn't have feedbacks  
    cfg.dataset             = [pathdat '/' name]; % ascii converted eyelink filename
    cfg.headerformat        = 'eyelink_asc';
    cfg.dataformat          = 'eyelink_asc';
    cfg.trialdef.eventtype  = 'MSG';
    cfg.trialdef.eventvalue = [1 2 3 4 5]; % triggers
    cfg.trialdef.prestim    = 0.2; % add start trigger fixation cross before the whole thing starts
    cfg.trialdef.poststim   = 2.5;
    cfg.channel             = {'4'}; % channel 2 is the x-coordinate
    % channel 3 is the y-coordinate
    % channel 4 is the pupil dilation
    cfg.dataformat          = 'eyelink_asc';
    cfg.headerformat        = 'eyelink_asc';
    cfg                     = ft_definetrial(cfg);
    eye_stim                = ft_preprocessing(cfg);
    
    % fb(reward) segment
    cfg                     = [];
    cfg.trialfun            = 'my_trialfunction_eyetracker_rewardtask_2sess_fb';
    cfg.dataset             = [pathdat '/' name]; % ascii converted eyelink filename
    cfg.headerformat        = 'eyelink_asc';
    cfg.dataformat          = 'eyelink_asc';
    cfg.trialdef.eventtype  = 'MSG';
    cfg.trialdef.eventvalue = [1 2 3 4 5]; % triggers
    cfg.trialdef.prestim    = 0.2; % add start trigger fixation cross before the whole thing starts
    cfg.trialdef.poststim   = 2.5;
    cfg.channel             = {'4'}; % channel 2 is the x-coordinate
    % channel 3 is the y-coordinate
    % channel 4 is the pupil dilation
    cfg.dataformat          = 'eyelink_asc';
    cfg.headerformat        = 'eyelink_asc';
    cfg                     = ft_definetrial(cfg);
    eye_fb                  = ft_preprocessing(cfg);
    
    %% pupil data
    
    save([pathdat '/' name(1:end-4) 'stim_eyedat.mat'],'eye_stim');
    fname_StimSeg = [pathdat '/' name(1:end-4) 'stim_eyedat.mat'];
    
    data_stim = cell2mat(eye_stim.trial');
    figure,errorbar(nanmedian(data_stim),nanstd(data_stim)./sqrt(size(data_stim,1))),title(name);
    h = gcf;
    saveas(h,[pathdat '/' name(1:end-4) 'stim_eyemedian'],'fig')
    figure,plot(nanmedian(data_stim)),title(name);
    
 

    save([pathdat '/' name(1:end-4) 'fb_eyedat.mat'],'eye_fb');
    fname_fbSeg = [pathdat '/' name(1:end-4) 'fb_eyedat.mat'];

    data_fb = cell2mat(eye_fb.trial');
    figure,errorbar(nanmedian(data_fb),nanstd(data_fb)./sqrt(size(data_fb,1))),title(name);
    h = gcf;
    saveas(h,[pathdat '/' name(1:end-4) 'fb_eyemedian'],'fig')
    figure,plot(nanmedian(data_fb)),title(name);
 
    
    close all
    
    fprintf('\n Eyetracking data imported and segmented for ID: %d\n',VPs(vps))
    
    %% load behavioral dataset
    
    eye_rv.eyedat.stim = eye_stim;
    eye_rv.eyedat.fb = eye_fb;
    
    if days == 1
        eye_rv.MainTaskResults = expdat{vps,1}.dat.day1.maintask;
        eye_rv.MemoryTestResults=expdat{vps,1}.dat.day1.memorytest;
%         eye_rv.trl_raus = trl_raus;
    else
        eye_rv.MainTaskResults = expdat{vps,1}.dat.day2.maintask;
        eye_rv.MemoryTestResults=expdat{vps,1}.dat.day2.memorytest;
%         eye_rv.accuracy = expdat{vps,1}.dat.day2.maintask.results.accuracy(trl_raus==0);
%         eye_rv.keypress = expdat{vps,1}.dat.day2.maintask.results.keypress(trl_raus==0);
%         eye_rv.rts      = expdat{vps,1}.dat.day2.maintask.results.rt(trl_raus==0);
%         eye_rv.trlord   = expdat{vps,1}.dat.day2.maintask.results.trl(trl_raus==0,:);
%         eye_rv.trl_raus = trl_raus;
    end
    
    save([path_save 're' num2str(days) num2str(VPs(vps)) '_eye_rv200.mat'],'eye_rv')
%     
%     cd(path_dat); cd ..
%     if exist('EB_cleaned')==7
%     else
%         mkdir('EB_cleaned')
%         cd EB_cleaned
%         mkdir('EBB_pupil')
%         cd EBB_pupil
%         mkdir('pilot_3T_2sess')
%     end
%     
%     save([path_save 'pilot_3T/re' num2str(days) num2str(VPs(vps)) '_eye_rv200_fb.mat'],'eye_rv','blinkbig','blinksmall')
%     keep cutleft cutright pathdat pathdatbehav fenster smallblinksize id VPs session expdat i2
    close all    
end
 
