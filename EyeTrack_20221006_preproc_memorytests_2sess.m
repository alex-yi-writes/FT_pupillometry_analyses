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
% days = [1 2; 1 2; 1 2; 1 2; 1 2; 0 2; 1 2; 1 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
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
        if days(i1,d) == 0
            fname_raweye_main{i1,d}  = []; % rewardtask-encoding (main task)
            fname_raweye_immediate{i1,d}  = []; % rewardtask-immediate recall
            fname_raweye_delayed{i1,d}  = []; % rewardtask-delayed recall
            fname_beh{i1,d} = [];
            expdat{i1,d} = [];
        elseif days(i1,d) == 1
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
        elseif days(i1,d) == 2
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
pathdat       = fullfile(path_dat,'EAB_pupil','MRPETpilot','ASC','immediate');
pathdatbehav  = fullfile(path_behav,fname_beh);

fprintf('\n Preparation done \n')
%% import eyetracker data

close all

for vps = [36]%1:length(VPs)

    for d = 2%1:2
        if days(vps,d) == 0
            fprintf('\n skipping day %d for ID %d\n', d, VPs(vps))
        else
            fprintf('\n importing and segmenting ID and day %d: %d\n', d, VPs(vps))
            for m = 1%:2

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

%                 try
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

%                 catch
%                 warning('no data')
%                 end
            end % close memory tests loop
        end % close session conditional


    end % close session loop
end % close participants loop
diary off
