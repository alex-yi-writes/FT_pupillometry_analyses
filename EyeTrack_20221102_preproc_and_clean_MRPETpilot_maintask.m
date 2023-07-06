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
    VPs = [2207];
%     VPs = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 ...
%      2127 2128 2129 2130 2131 2132 2233 2234 2235 2236 2237 2238 2239 2240 2142 2143 2144 2145 2146 2147 2148 2249 2250]; % 2203: empty
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
 
for vps = 1:length(VPs)
    
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
    
    %% clean: stimulus segment
    
    cutleft = 1; cutright = 0; fenster = 30; smallblinksize = 200;
    
    fprintf('\n cleaning (stimulus segment) ID: %d\n', VPs(vps))
    
    datacall = []; dataiall = []; trls = []; samples = [];
    
    load([pathdat '/' 're' num2str(days) num2str(VPs(vps)) 'stim_eyedat.mat']);
    
    data_stim    = cell2mat(eye_stim.trial');
    samples = [samples;eye_stim.sampleinfo+1*1*10^6];%%% newly added, should be done with all future datasets
    trls    = [trls;eye_stim.cfg.trl+1*1*10^6];%%% newly added, should be done with all future datasets
    
    %  exclude data parts with closed eyes
    for tr = 1:size(data_stim,1)
        
        c = 0; emptySa = [];
        temp = data_stim(tr,:);
        x    = 1:length(temp);
        
        for sa = cutleft+1:length(temp)-(cutright+5)
            
            if temp(sa) == 0 & temp(sa-1) ~= 0
                cc = 0;c = c+1;
                
            elseif (sa == cutleft+1) & (temp(cutleft+1) == 0)% starting in blink
                c = c+1;cc = 1;
                emptySa{c}(cc) = sa;
                
            elseif temp(sa) == 0 & temp(sa-1) == 0 % continuous streaks of blinks
                cc = cc+1;
                emptySa{c}(cc) = sa;
                
            end
            
        end
        clear c cc
        
        if ~isempty(emptySa)
            
            for e =1:size(emptySa,2)
                temp(min(emptySa{e}) - cutleft:max(emptySa{e}) + cutright) = NaN; % make 200ms around blink nan
            end
            
            datac(tr,1:size(data_stim,2)) = temp(1:size(data_stim,2));
            blinkbig(tr) = sum(isnan(temp)); % calculate how many ms exclued because of blinks
            
            % if beginning missing, replace with mean
            if isnan(temp(1))
                temp(1) = nanmean(temp);
            end
            
            % interpolate data
            xi      = 1:length(temp);
            zs      = isnan(temp);
            temp(zs) = []; x(zs)=[];
            temp2   = interp1(x, temp, xi);
            datai(tr,1:size(data_stim,2)) = temp2(1:size(data_stim,2));
            
        elseif isempty(emptySa)
            
            datac(tr,:) = temp;
            datai(tr,:) = temp;
            
        end
        clear emptySa e temp temp2
        
        
        % and remove tiny blinks in interpolated data and interpolate
        % again
        
        jumps = [];
        datasel = datai(tr,:);
        cc = 0;
        
        for w = 1:fenster:length(datasel)-fenster
            
            cc  = cc+1;
            d   = datasel(w:w+fenster-1);
            jumps(cc,1) = sum(max(d)-min(d));
            jumps(cc,2) = w;
            
        end
        
        
        blinked = find(jumps(:,1) > smallblinksize);
        if ~isempty(blinked)
            
            cuttimes = unique([blinked-2:[blinked(end)+2]]); % take two windows before and afterwards
            cuttimes(find(cuttimes < 1)) = [];
            cuttimes(find(cuttimes > length(jumps))) = []; % remove parts that are outside of the range
            
            % replace separately for each blink
            blinknums = 1; cutcount = 0;
            blinkrepl{blinknums}(1:2) = jumps(cuttimes(1:2),2); % first cut
            
            for b = 3:length(cuttimes)
                
                if cuttimes(b) <= length(jumps)% not if end reached
                    
                    cutcount = cutcount +1;
                    blinkrepl{blinknums}(cutcount) = jumps(cuttimes(b),2);
                    
                    if cuttimes(b-1) - cuttimes(b) < -1 % not adjacent
                        cutcount = 0;% start new blink count
                        blinknums = blinknums +1;
                    end
                    
                end
                
            end
            clear cutcount b
            
            temp = datai(tr,:);
            for bb = 1:blinknums
                temp(min(blinkrepl{bb}):max(blinkrepl{bb})) = NaN; % replace blnks with nans
            end
            datac(tr,:)     = temp;
            blinksmall(tr)  = length(cuttimes)*fenster;
            clear bb blinknums blinkrepl
            
            % if beginning missing, replace with mean
            if isnan(temp(1))
                temp(1) = nanmean(temp);
            end
            
            % cut to correct size and interpolate data
            xi  = 1:length(temp);
            x   = 1:length(temp);
            zs  = isnan(temp);
            temp(zs) = []; x(zs)=[];
            temp2 = interp1(x, temp, xi);
            datai(tr,1:size(data_stim,2)) = temp2(1:size(data_stim,2));
            
        end
        
    end
    

    %% check
    
    figure,
    for i = 1:80
        subplot(10,8,i);plot(data_stim(i,:));hold on
        plot(datai(i,:),'r');hold off
    end
    datacall = [datacall; datac]; clear datac
    dataiall = [dataiall; datai]; clear datai
    
    %%
    
    eyeinterp = eye_stim;
    
    % make data structure new, all three runs together, and baselinecorrect
    for tr = 1: size(dataiall,1)
        
        temp    = zscore(dataiall(tr,:));
        basel   = nanmean(temp(1:200)); % 200 ms before stim onset
        eyeinterp.trial{1,tr}   = temp-basel; clear temp
        eyeinterp.time{1,tr}    = eye_stim.time{1};
        
    end
    
    eyeinterp.cfg.trl    = trls;
    eyeinterp.sampleinfo = samples;
    
    
    %% reject visual
    
    cfg.method = 'channel';
    cfg.alim = 4;
    % cfg.metric = 'var';
    cfg.latency = [-.2 2.5];
    eyeinterp = ft_rejectvisual(cfg, eyeinterp)
    
    
    cfg.method = 'trial';
    cfg.alim = 4;
    % cfg.metric = 'var';
    cfg.latency = [-.2 2.5];
    eye_rv = ft_rejectvisual(cfg, eyeinterp)
    
    
    % determine which trials out
    trlnew = eye_rv.sampleinfo(:,1);
    trlold = eye_rv.cfg.previous.previous.trl(:,1);
    trl_raus = ~ismember(trlold,trlnew);
    sum(trl_raus)
    
    %% load behavioral dataset
    
    if days == 1
        eye_rv.accuracy = expdat{vps,1}.dat.day1.maintask.results.accuracy(trl_raus==0);
        eye_rv.keypress = expdat{vps,1}.dat.day1.maintask.results.keypress(trl_raus==0);
        eye_rv.rts      = expdat{vps,1}.dat.day1.maintask.results.rt(trl_raus==0);
        eye_rv.trlord   = expdat{vps,1}.dat.day1.maintask.results.trl(trl_raus==0,:);
        eye_rv.trl_raus = trl_raus;
    else
        eye_rv.accuracy = expdat{vps,1}.dat.day2.maintask.results.accuracy(trl_raus==0);
        eye_rv.keypress = expdat{vps,1}.dat.day2.maintask.results.keypress(trl_raus==0);
        eye_rv.rts      = expdat{vps,1}.dat.day2.maintask.results.rt(trl_raus==0);
        eye_rv.trlord   = expdat{vps,1}.dat.day2.maintask.results.trl(trl_raus==0,:);
        eye_rv.trl_raus = trl_raus;
    end

    
    save([path_save 're' num2str(days) num2str(VPs(vps)) '_eye_rv200_stim.mat'],'eye_rv','blinkbig','blinksmall')
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


    %% clean: fb segment
    
    cutleft = 1; cutright = 0; fenster = 30; smallblinksize = 200;
    
    fprintf('\n cleaning (fb segment) ID: %d\n', VPs(vps))
    
    datacall = []; dataiall = []; trls = []; samples = [];
    
    load([pathdat '/' 're' num2str(days) num2str(VPs(vps)) 'fb_eyedat.mat']);
    
    data_fb    = cell2mat(eye_fb.trial');
    samples = [samples;eye_fb.sampleinfo+1*1*10^6];%%% newly added, should be done with all future datasets
    trls    = [trls;eye_fb.cfg.trl+1*1*10^6];%%% newly added, should be done with all future datasets
    
    %  exclude data parts with closed eyes
    for tr = 1:size(data_fb,1)
        
        c = 0; emptySa = [];
        temp = data_fb(tr,:);
        x    = 1:length(temp);
        
        for sa = cutleft+1:length(temp)-(cutright+5)
            
            if temp(sa) == 0 & temp(sa-1) ~= 0
                cc = 0;c = c+1;
                
            elseif (sa == cutleft+1) & (temp(cutleft+1) == 0)% starting in blink
                c = c+1;cc = 1;
                emptySa{c}(cc) = sa;
                
            elseif temp(sa) == 0 & temp(sa-1) == 0 % continuous streaks of blinks
                cc = cc+1;
                emptySa{c}(cc) = sa;
                
            end
            
        end
        clear c cc
        
        if ~isempty(emptySa)
            
            for e =1:size(emptySa,2)
                temp(min(emptySa{e}) - cutleft:max(emptySa{e}) + cutright) = NaN; % make 200ms around blink nan
            end
            
            datac(tr,1:size(data_fb,2)) = temp(1:size(data_fb,2));
            blinkbig(tr) = sum(isnan(temp)); % calculate how many ms exclued because of blinks
            
            % if beginning missing, replace with mean
            if isnan(temp(1))
                temp(1) = nanmean(temp);
            end
            
            % interpolate data
            xi      = 1:length(temp);
            zs      = isnan(temp);
            temp(zs) = []; x(zs)=[];
            temp2   = interp1(x, temp, xi);
            datai(tr,1:size(data_fb,2)) = temp2(1:size(data_fb,2));
            
        elseif isempty(emptySa)
            
            datac(tr,:) = temp;
            datai(tr,:) = temp;
            
        end
        clear emptySa e temp temp2
        
        
        % and remove tiny blinks in interpolated data and interpolate
        % again
        
        jumps = [];
        datasel = datai(tr,:);
        cc = 0;
        
        for w = 1:fenster:length(datasel)-fenster
            
            cc  = cc+1;
            d   = datasel(w:w+fenster-1);
            jumps(cc,1) = sum(max(d)-min(d));
            jumps(cc,2) = w;
            
        end
        
        
        blinked = find(jumps(:,1) > smallblinksize);
        if ~isempty(blinked)
            
            cuttimes = unique([blinked-2:[blinked(end)+2]]); % take two windows before and afterwards
            cuttimes(find(cuttimes < 1)) = [];
            cuttimes(find(cuttimes > length(jumps))) = []; % remove parts that are outside of the range
            
            % replace separately for each blink
            blinknums = 1; cutcount = 0;
            blinkrepl{blinknums}(1:2) = jumps(cuttimes(1:2),2); % first cut
            
            for b = 3:length(cuttimes)
                
                if cuttimes(b) <= length(jumps)% not if end reached
                    
                    cutcount = cutcount +1;
                    blinkrepl{blinknums}(cutcount) = jumps(cuttimes(b),2);
                    
                    if cuttimes(b-1) - cuttimes(b) < -1 % not adjacent
                        cutcount = 0;% start new blink count
                        blinknums = blinknums +1;
                    end
                    
                end
                
            end
            clear cutcount b
            
            temp = datai(tr,:);
            for bb = 1:blinknums
                temp(min(blinkrepl{bb}):max(blinkrepl{bb})) = NaN; % replace blnks with nans
            end
            datac(tr,:)     = temp;
            blinksmall(tr)  = length(cuttimes)*fenster;
            clear bb blinknums blinkrepl
            
            % if beginning missing, replace with mean
            if isnan(temp(1))
                temp(1) = nanmean(temp);
            end
            
            % cut to correct size and interpolate data
            xi  = 1:length(temp);
            x   = 1:length(temp);
            zs  = isnan(temp);
            temp(zs) = []; x(zs)=[];
            temp2 = interp1(x, temp, xi);
            datai(tr,1:size(data_fb,2)) = temp2(1:size(data_fb,2));
            
        end
        
    end
    
    %% check
    
    figure,
    for i = 1:80
        subplot(10,8,i);plot(data_stim(i,:));hold on
        plot(datai(i,:),'r');hold off
    end
    datacall = [datacall; datac]; clear datac
    dataiall = [dataiall; datai]; clear datai
    
    %%
    
    eyeinterp = eye_fb;
    
    % make data structure new, all three runs together, and baselinecorrect
    for tr = 1: size(dataiall,1)
        
        temp    = zscore(dataiall(tr,:));
        basel   = nanmean(temp(1:200)); % 200 ms before fb onset
        eyeinterp.trial{1,tr}   = temp-basel; clear temp
        eyeinterp.time{1,tr}    = eye_fb.time{1};
        
    end
    
    eyeinterp.cfg.trl    = trls;
    eyeinterp.sampleinfo = samples;
    
    
    %% reject visual
    
    cfg.method = 'channel';
    cfg.alim = 4;
    % cfg.metric = 'var';
    cfg.latency = [-.2 2.5];
    eyeinterp = ft_rejectvisual(cfg, eyeinterp)
    
    
    cfg.method = 'trial';
    cfg.alim = 4;
    % cfg.metric = 'var';
    cfg.latency = [-.2 2.5];
    eye_rv = ft_rejectvisual(cfg, eyeinterp)
    
    
    % determine which trials out
    trlnew = eye_rv.sampleinfo(:,1);
    trlold = eye_rv.cfg.previous.previous.trl(:,1);
    trl_raus = ~ismember(trlold,trlnew);
    sum(trl_raus)
    
    %% load behavioral dataset
    
    if days == 1
        eye_rv.accuracy = expdat{vps,1}.dat.day1.maintask.results.accuracy(trl_raus==0);
        eye_rv.keypress = expdat{vps,1}.dat.day1.maintask.results.keypress(trl_raus==0);
        eye_rv.rts      = expdat{vps,1}.dat.day1.maintask.results.rt(trl_raus==0);
        eye_rv.trlord   = expdat{vps,1}.dat.day1.maintask.results.trl(trl_raus==0,:);
        eye_rv.trl_raus = trl_raus;
    else
        eye_rv.accuracy = expdat{vps,1}.dat.day2.maintask.results.accuracy(trl_raus==0);
        eye_rv.keypress = expdat{vps,1}.dat.day2.maintask.results.keypress(trl_raus==0);
        eye_rv.rts      = expdat{vps,1}.dat.day2.maintask.results.rt(trl_raus==0);
        eye_rv.trlord   = expdat{vps,1}.dat.day2.maintask.results.trl(trl_raus==0,:);
        eye_rv.trl_raus = trl_raus;
    end
    

    
    save([path_save 're' num2str(days) num2str(VPs(vps)) '_eye_rv200_fb.mat'],'eye_rv','blinkbig','blinksmall')
%     keep cutleft cutright pathdat pathdatbehav fenster smallblinksize id VPs session expdat i2
    close all   
end
 
