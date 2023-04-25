% import eyetracker data for continuous pupil assessment

% how to best link pupil and fmri data?

% time point of first stimulus - MSG = 100
% then just set segmentation points every 3.04 sec , shifted by 1 sec (so
% 1-2.04
% allows to infer when first volume was let's get timelock of stim presentation from MR rergessors and link it to time of first stim presentation,
% then cut volumes accordingly


% clear
close all
% addpath(genpath('/Users/dorothea/Dropbox/matlab/toolboxes/fieldtrip-20170119'));
% ft_defaults
% addpath(genpath('/Users/dorothea/Dropbox/matlab/toolboxes/spm12/'));
% pathdat = ['/Users/dorothea/Dropbox/eyetracker_data/2opt_fmri/'];

pathdat = ['/Users/dorothea/Desktop/'];

load('/Users/dorothea/Documents/Mitarbeiter/Mareike_Mastearbeit/Skripte/ids.mat');

% not all of these will have pupil data
rangeids = ids;

TR = 3.04;

% 3109 - convert first ascii file agian .. weird )
% 3271 4108 4 doesn't work ...)
suff = 3;% % 2d or 3d

% settings for interpolation of missing data
cutleft = 200;cutright = 200;fenster = 30; smallblinksize = 200;

done = [];
for idd = 31%:length(rangeids)
%     try
        for p = 4%:6
            blinkbig = []; blinksmall = [];
%             try
                % change trigger name
                name = [num2str(rangeids(idd)) '_' num2str(p) num2str(suff) '.asc']
                eyeAdjustTrigNam_continuous([pathdat name],'MSG','INPUT');% change MSG to INPUT for trials
                
                [type] = ft_filetype([pathdat name])
                
                cfg                     = [];
                cfg.trialfun            = 'my_trialfunction_eyetracker_continuous';
                cfg.dataset             = [pathdat name]; % ascii converted eyelink filename
                cfg.headerformat        = 'eyelink_asc';
                cfg.dataformat          = 'eyelink_asc';
                cfg.trialdef.eventtype  = 'MSG';
                cfg.trialdef.eventvalue = [1 2 3]; % triggers
                cfg.trialdef.prestim    = .2; % add start trigger fixation cross before the whole thing starts
                cfg.trialdef.poststim   = TR;
                cfg.channel             = {'4'}; % channel 2 is the x-coordinate
                % channel 3 is the y-coordinate
                % channel 4 is the pupil dilation
                cfg.dataformat          = 'eyelink_asc';
                cfg.headerformat        = 'eyelink_asc';
                
                
                cfg                     = ft_definetrial(cfg);
                eye                     = ft_preprocessing(cfg);
                
                %% pupil data
                data = cell2mat(eye.trial');
                datamat = data;clear data
                % maybe put data as continuum together for interpolating, and then separate into trial again (with baseline)
                data = [];
                data = [data datamat(1,1:end)];% first trial with baseline
                for trls = 2:size(datamat,1)
                    data = [data datamat(trls,202:end)];
                end
                %
                %  exclude data parts with closed eyes and interpolate curves
                for tr = 1:size(data,1)
                    c = 0; emptySa = [];
                    temp = data(tr,:);
                    x = 1:length(temp);
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
                            temp(min(emptySa{e})-cutleft:max(emptySa{e})+cutright) = NaN;% make 200ms around blink nan
                        end
                        datac(tr,1:size(data,2)) = temp(1:size(data,2)); % datac is data with mssing data set NaN
                        blinkbig(tr)  = sum(isnan(temp));% calculate how many ms exclued because of blinks
                        % if beginning missing, replace with mean
                        if isnan(temp(1))
                            temp(1) = nanmean(temp);
                        end
                        % interpolate missing data
                        xi = 1:length(temp);
                        zs = isnan(temp);
                        temp(zs) = [];
                        x(zs)=[];
                        temp2 = interp1(x, temp, xi);
                        datai(tr,1:size(data,2)) = temp2(1:size(data,2));  % datai is data with mssing data interpolated
                    elseif isempty(emptySa)
                        datac(tr,:) = temp;
                        datai(tr,:) = temp;
                    end
                    clear emptySa e temp temp2
                    %%%%% new
                    % and remove tiny blinks in interpolated data and interpolate
                    % again
                    jumps = [];
                    datasel = datai(tr,:);
                    cc = 0;
                    for w = 1:fenster:length(datasel)-fenster % move through it in TR intervals?
                        cc = cc+1;
                        d = datasel(w:w+fenster-1);
                        jumps(cc,1) = sum(max(d)-min(d));
                        jumps(cc,2) = w;
                    end
                    blinked = find(jumps(:,1) > smallblinksize); % select starting points of 30ms time windows with blinks
                    
                    vec_fenster = 1:fenster:length(datasel)-fenster;
                    % make blink periods NaN (with 2 sec before and after) and interpolate
                    cuttimes = [];
                    if ~isempty(blinked)
                        for bb = 1:length(blinked)
                            cuttimes = [cuttimes vec_fenster(blinked(bb)):vec_fenster(blinked(bb))+fenster];
                        end
                    end
                    clear bb
                    dataout = unique(cuttimes);dataout(find(dataout==0)) = [];
                    temp = datai(tr,:);
                    temp(dataout) = NaN;
                    clear datai
                    blinksmall= length(dataout); % how many sample rejected because of small blinks?
                    % if beginning missing, replace with mean
                    if isnan(temp(1))
                        temp(1) = nanmean(temp);
                    end
                    % cut to correct size and interpolate data
                    xi = 1:length(temp);
                    x = 1:length(temp);
                    zs = isnan(temp);
                    temp(zs) = [];
                    x(zs)=[];
                    temp2 = interp1(x, temp, xi);
                    datai(tr,1:size(data,2)) = temp2(1:size(data,2));
                end
                %% %%%%%%%%%%% can be commented out, is just showing old (blue) and red (interpolated) data
                figure,
                plot(data(1:60000));hold on
                plot(datai(1:60000),'r');hold off
                %% %%%%%%%%%%% can be commented out, is just showing old (blue) and red (interpolated) data
                %% make new data structure from interpolated data for manual rejection, using the original segment sizes again
                eyenew = eye;
                % make data structure new, and baseline correct before rejection trials
                % individually
                eyenew.trial{1,1}  = datai(1:cfg.trialdef.poststim*eye.fsample+1-eye.cfg.trl(1,3)) % add prestimtime
                eyenew.time{1,tr} = eye.time{1};
                steps = cfg.trialdef.poststim*eye.fsample; % how long is one volume
                for tr = 2:size(datamat,1)
                    eyenew.trial{1,tr} = datai([(tr-1)*steps+eye.cfg.trl(1,3)]:[(tr-1)*steps+cfg.trialdef.poststim*eye.fsample]);
                    eyenew.time{1,tr} = eye.time{1};
                end
                save([pathdat name(1:end-4) '_eyedat_cont_run' num2str(p) '.mat'],'eyenew','eye','blinksmall','blinkbig');
                keep names n pathdat idd rangeids suff done tasks p TR cutleft cutright fenster smallblinksize
                done(idd,p) = rangeids(idd)
%             catch end
        end
%     catch end
    close all
end
