%% Pupillometric data preprocessing: recall phase
%   memory tests eyedata

%% work log

%   27_11_2018 created the script
%              now incorporated an optimised data structure

%% preparation

clear;clc;close all; warning('on')

% path setup
addpath(genpath('/Users/yeojin/Documents/MATLAB/fieldtrip-20220104')); % change those two
addpath(genpath('/Users/yeojin/Desktop/B_scripts/BA_preprocessing/'))

% which data will you be processing? (for now, it's all 1)
% answ        = inputdlg({'which session? 1 or 2 (for now, it''s all 1)'});
% sessions    = str2num(answ{1,1}); clear answ
% answ        = inputdlg({'which pilot?'});
% pilotNr     = str2num(answ{1,1}); clear answ
% answ        = inputdlg({'which memory test? 1=immediate or 2=delayed'});
% memtest     = str2num(answ{1,1}); clear answ
pilotNr     = 2;

path_parent = '/Users/yeojin/Desktop/E_data/';
path_dat    = [ path_parent 'EA_raw/' ]; % tailor it to your original path
path_stim   = [ path_dat 'EAA_pupil/stiminfo_all_20181004.mat' ]; % stimulus information
path_behav  = [ path_dat 'EAC_behav/MRPETpilot/' ];
path_save   = [ path_parent 'EB_cleaned/EBB_pupil/MRPETpilot/'];

% create file list and VP ID list
% VPs = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2111 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221];
% tag = [1 2; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
% d1m = [1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 0; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
% d2m = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
% IDs
VPs  = [2202 2203 2204 2205 2206 2207 2208 2109 2110 2112 2113 2114 2115 2116 2217 2218 2219 2220 2221 2222 2223 2224 2125 2126 2127 2128 2129 2130 2131 2132 2233 2234 2235 2236 2237 2238 2239 2240 2142 2143 2144 2145 2146 2147 2148 2249 2250];
days = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
d1m  = [1 2; 1 0; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2]; % 1=immediate 2=delayed
d2m  = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

% if sessions == 1
%     VPs = [2202 2203 2204 2207 2208];
% else
%     VPs = [2201 2202 2205 2206];
% end

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
            fname_beh{i1,d}     = [ num2str(VPs(i1)) '_' num2str(d) '.mat' ];
            expdat{i1,d} = load([path_behav fname_beh{i1,d}]);
        end
    end
    
end

% script specific variables
pathdat       = fullfile(path_dat,'EAB_pupil','MRPETpilot','ASC');
pathdatbehav  = fullfile(path_behav,fname_beh);

fprintf('\n Preparation done \n')
%% import eyetracker data

close all

for vps = 1:length(VPs)
    
    for d = 1:2
        if tag(vps,d) == 0
            fprintf('\n skipping day %d for ID %d\n', d, VPs(vps))
        else
            fprintf('\n importing and segmenting ID and day %d: %d\n', d, VPs(vps))
            for m = 1:2
                
                if eval(['d' num2str(d) 'm(vps,m) == 1'])
                    fprintf('\n***immediate***\n')
                    name = fname_raweye_immediate{vps,d}
                elseif eval(['d' num2str(d) 'm(vps,m) == 2'])
                    fprintf('\n***delayed***\n')
                    name = fname_raweye_delayed{vps,d}
                else
                    name = {};
                    if m == 1; errorname = 'immediate'; else errorname = 'delayed';end
                    fprintf(['no data for day %d for ID %d %s test\n'], d, VPs(vps),errorname)
                end
                
                try
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
                    
                    i1 = 0; kk1 = 0; clear startline endline
                    for k=1:numel(C{1,1})
                        tmp1 = regexp(C{1,1}(k),'CAMERA_CONFIG'); % find where header ends
                        if ~isemptycell(tmp1)
                            i1 = i1+1;
                            startline(i1) = k+2; % mark the start of lines deleted
                        end
                        tmp2 = regexp(C{1,1}(k),'RECCFG CR '); % find where the recording starts
                        %         tmp2 = regexp(C{1,1}(k),'RECCFG CR 1000 2 1 R'); % find where the recording starts
                        if ~isemptycell(tmp2)
                            kk1 = kk1+1;
                            endline(kk1) = k-1; % mark the end of lines deleted
                        end
                        
                    end
                    fprintf('\n importing done \n')
                    
                    newtext = [];
                    if ~isempty(startline) && ~isempty(endline)
                        
                        % okay somehow this bit doesn't work, and i'm too tired now... i'll
                        % just write some workarounds
                        %         for k1 = startline(1):endline
                        %             C{1,1}{k1,1} = []; % now delete those lines containing calibration info
                        %         end
                        %         for k1 = startline(1):endline
                        %             if isemptycell(C{1,1}(k1))
                        %                 C{1,1}(k1)= [];
                        %             else
                        %             end% now delete those lines containing calibration info
                        %         end
                        
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
                    
                    
                    eyeAdjustTrigNam([pathdat '/' name],'MSG','INPUT');
                    [type] = ft_filetype([pathdat '/' name])
                    
                    % stimulus segment
                    cfg                     = [];
                    
                    if m == 1
                        cfg.trialfun            = 'my_trialfunction_eyetracker_memorytest_3T';
                        cfg.trialdef.eventvalue = [1 2 3]; % triggers
                    elseif m == 2
                        cfg.trialfun            = 'my_trialfunction_eyetracker_memorytest_3T_delayed';
                        cfg.trialdef.eventvalue = [1 2 3]; % triggers
                    end
                    cfg.dataset             = [pathdat '/' name]; % ascii converted eyelink filename
                    cfg.headerformat        = 'eyelink_asc';
                    cfg.dataformat          = 'eyelink_asc';
                    cfg.trialdef.eventtype  = 'MSG';
                    
                    cfg.trialdef.prestim    = 0.2; % add start trigger fixation cross before the whole thing starts
                    cfg.trialdef.poststim   = 2.5;
                    cfg.channel             = {'4'}; % channel 2 is the x-coordinate
                    % channel 3 is the y-coordinate
                    % channel 4 is the pupil dilation
                    cfg.dataformat          = 'eyelink_asc';
                    cfg.headerformat        = 'eyelink_asc';
                    
                    cfg                     = ft_definetrial(cfg);
                    eye_mem                 = ft_preprocessing(cfg);
                    
                    
                    %     % resample if necessary
                    %     sampleinfo = eye_mem.sampleinfo;
                    %     if eye_mem.fsample ~= 1000
                    %         cfg = [];
                    %         cfg.resamplefs = 1000;
                    %         cfg.detrend = 'no';
                    %         eye_mem = ft_resampledata(cfg, eye_mem);
                    %         eye_mem.cfg = eye_mem.cfg.previous;
                    %         eye_mem.sampleinfo = sampleinfo;
                    %     end
                    
                    
                    
                    % pupil data
                    
                    save([pathdat '/' name(1:end-4) 'mem_eyedat.mat'],'eye_mem');
                    fname_memSeg = [pathdat '/' name(1:end-4) 'mem_eyedat.mat'];
                    try
                        data_mem = cell2mat(eye_mem.trial');
                    catch
                        tmp4 = eye_mem.trial;
                        for imc = 1:length(tmp4)
                            if length(tmp4{1,imc}) < 2701
                                tmp4{1,imc} = horzcat(tmp4{1,imc}(1,:),tmp4{1,imc}(end));
                                %                 tmp4{1,imc}(end) = [];
                            else
                            end
                        end
                        data_mem = cell2mat(tmp4)';
                        eye_mem.trial = tmp4; clear tmp4
                    end
                    
                    %     data_mem = cell2mat(eye_mem.trial');
                    figure,errorbar(nanmedian(data_mem),nanstd(data_mem)./sqrt(size(data_mem,1))),title(name);
                    h = gcf;
                    saveas(h,[pathdat '/' name(1:end-4) 'mem_eyemedian'],'fig')
                    figure,plot(nanmedian(data_mem)),title(name);
                    
                    
                    close all
                    
                    fprintf('\n Eyetracking data imported and segmented for ID: %d\n',VPs(vps))
                    
                    %% clean
                    
                    cutleft = 1; cutright = 0; fenster = 30; smallblinksize = 200;
                    
                    fprintf('\n cleaning ID: %d\n', VPs(vps))
                    
                    datacall = []; dataiall = []; trls = []; samples = [];
                    
                    if m == 1
                        load([pathdat '/' 'ri' num2str(d) num2str(VPs(vps)) 'mem_eyedat.mat']);
                    elseif m == 2
                        load([pathdat '/' 'rd' num2str(d) num2str(VPs(vps)) 'mem_eyedat.mat']);
                    end
                    
                    try
                        data_mem = cell2mat(eye_mem.trial');
                    catch
                        tmp4 = eye_mem.trial;
                        for imc = 1:length(tmp4)
                            if length(tmp4{1,imc}) < 2701
                                tmp4{1,imc} = horzcat(tmp4{1,imc}(1,:),tmp4{1,imc}(end));
                                %                 tmp4{1,imc}(end) = [];
                            else
                            end
                        end
                        data_mem = cell2mat(tmp4)';
                        eye_mem.trial = tmp4; clear tmp4
                    end
                    
                    %     data_mem    = cell2mat(eye_mem.trial');
                    %     samples = [samples;eye.sampleinfo+(p-1)*1*10^6];%%% newly added, should be done with all future datasets
                    %     trls    = [trls;eye.cfg.trl+(p-1)*1*10^6];%%% newly added, should be done with all future datasets
                    samples = [samples;eye_mem.sampleinfo+1*1*10^6];%%% newly added, should be done with all future datasets
                    trls    = [trls;eye_mem.cfg.trl+1*1*10^6];%%% newly added, should be done with all future datasets
                    
                    %  exclude data parts with closed eyes
                    for tr = 1:size(data_mem,1)
                        
                        c = 0; emptySa = [];
                        temp = data_mem(tr,:);
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
                            
                            datac(tr,1:size(data_mem,2)) = temp(1:size(data_mem,2));
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
                            datai(tr,1:size(data_mem,2)) = temp2(1:size(data_mem,2));
                            
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
                            d1   = datasel(w:w+fenster-1);
                            jumps(cc,1) = sum(max(d1)-min(d1));
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
                            datai(tr,1:size(data_mem,2)) = temp2(1:size(data_mem,2));
                            
                        end
                        
                    end
                    
                    
                    %% check
                    
%                     figure,
%                     for i = 1:80
%                         subplot(10,8,i);plot(data_mem(i,:));hold on
%                         plot(datai(i,:),'r');hold off
%                     end
                    datacall = [datacall; datac]; clear datac
                    dataiall = [dataiall; datai]; clear datai
                    
                    %%
                    
                    eyeinterp = eye_mem;
                    
                    % make data structure new, all three runs together, and baselinecorrect
                    for tr = 1: size(dataiall,1)
                        
                        temp    = zscore(dataiall(tr,:));
                        basel   = nanmean(temp(1:200)); % 200 ms before stim onset
                        eyeinterp.trial{1,tr}   = temp-basel; clear temp
                        eyeinterp.time{1,tr}    = eye_mem.time{1};
                        
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
                    
                    % resample if necessary
                    sampleinfo = eye_rv.sampleinfo;
                    if eye_rv.fsample ~= 1000
                        cfg = [];
                        cfg.resamplefs = 1000;
                        cfg.detrend = 'no';
                        eye_rv = ft_resampledata(cfg, eye_rv);
                        eye_rv.cfg = eye_rv.cfg.previous;
                        eye_rv.sampleinfo = sampleinfo;
                    end
                    
                    %% load behavioral dataset
                    
                    cd(path_dat); cd ..
                    if exist('EB_cleaned')==7
                    else
                        mkdir('EB_cleaned')
                        cd EB_cleaned
                        mkdir('EBB_pupil')
                        cd EBB_pupil
                        mkdir('pilot_3T')
                    end
                    
                    if d == 1
                        if m == 1
                            eye_rv.immediate.memresults = expdat{vps,1}.dat.day1.memorytest.immediate.results.all(trl_raus==0,:);
                            eye_rv.immediate.memstimvec = expdat{vps,1}.dat.day1.memorytest.immediate.config.stim.stimlist_all(trl_raus==0,:);
                            eye_rv.immediate.memlabels  = expdat{vps,1}.dat.day1.memorytest.labels;
                            eye_rv.immediate.trl_raus   = trl_raus;
                            save([path_save 'pilot_3T/ri' num2str(d) num2str(VPs(vps)) '_eye_rv200_stim.mat'],'eye_rv','blinkbig','blinksmall')
                        elseif m == 2
                            eye_rv.delayed.memresults = expdat{vps,1}.dat.day1.memorytest.delayed.results.all(trl_raus==0,:);
                            eye_rv.delayed.memstimvec = expdat{vps,1}.dat.day1.memorytest.delayed.config.stim.stimlist_all(trl_raus==0,:);
                            eye_rv.delayed.memlabels  = expdat{vps,1}.dat.day1.memorytest.labels;
                            eye_rv.delayed.trl_raus   = trl_raus;
                            try
                                save([path_save 'pilot_3T/rd' num2str(d) num2str(VPs(vps)) '_eye_rv200_stim.mat'],'eye_rv','blinkbig','blinksmall')
                            catch
                                save([path_save 'pilot_3T/rd' num2str(d) num2str(VPs(vps)) '_eye_rv200_stim.mat'],'eye_rv','blinked','blinksmall')
                                
                            end
                        end
                    elseif d == 2
                        if m == 1
                            eye_rv.immediate.memresults = expdat{vps,2}.dat.day2.memorytest.immediate.results.all(trl_raus==0,:);
                            eye_rv.immediate.memstimvec = expdat{vps,2}.dat.day2.memorytest.immediate.config.stim.stimlist_all(trl_raus==0,:);
                            eye_rv.immediate.memlabels  = expdat{vps,2}.dat.day2.memorytest.labels;
                            eye_rv.immediate.trl_raus   = trl_raus;
                            save([path_save 'pilot_3T/ri' num2str(d) num2str(VPs(vps)) '_eye_rv200_stim.mat'],'eye_rv','blinkbig','blinksmall')
                        elseif m == 2
                            eye_rv.delayed.memresults = expdat{vps,2}.dat.day2.memorytest.delayed.results.all(trl_raus==0,:);
                            eye_rv.delayed.memstimvec = expdat{vps,2}.dat.day2.memorytest.delayed.config.stim.stimlist_all(trl_raus==0,:);
                            eye_rv.delayed.memlabels  = expdat{vps,2}.dat.day2.memorytest.labels;
                            eye_rv.delayed.trl_raus   = trl_raus;
                            try
                                save([path_save 'pilot_3T/rd' num2str(d) num2str(VPs(vps)) '_eye_rv200_stim.mat'],'eye_rv','blinkbig','blinksmall')
                            catch
                                save([path_save 'pilot_3T/rd' num2str(d) num2str(VPs(vps)) '_eye_rv200_stim.mat'],'eye_rv','blinked','blinksmall')
                                
                            end
                        end
                    end
                    
                    keep cutleft cutright pathdat pathdatbehav fenster smallblinksize id VPs vps d expdat tag m d1m d2m fname_raweye_main fname_raweye_immediate fname_raweye_delayed fname_beh path_save path_dat
                    close all
                catch
                    if isempty(name)
                        warning('no data for this particular date and test')
                    else
                        error('errors!')
                    end
                end
            end % close memory tests loop
        end % close session conditional
        
        
    end % close session loop
end % close participants loop
diary off
