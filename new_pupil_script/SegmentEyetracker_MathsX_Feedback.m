% import and segment eyetracker data
clear
close all
clc

% for e = 1:9:length(event)
%     event(e).value = 1; %cue (named t_cue, confusingly)
%     event(e+1).value = 2;% fixation
%     event(e+2).value = 3;% task
%     event(e+3).value = 4; % solution
%     event(e+4).value = 5; % fixation
%     event(e+5).value = 6; % difficulty question
%     event(e+6).value = 7; % fixation 
%     event(e+7).value = 8; % feedback
%     event(e+8).value = 9; % fixation (ISI)
% end

% addpath(genpath('\\imed\fme\kliniken\iknd\rakshit\Eigene Dateien\Dorothea\Pupil_clining\fieldtrip-master\fieldtrip-master'));
% pathdat = ['/Volumes/LIDO_DH/LIDO/DATA/'];
%pathdat =['/Users/dorothea/OneDrive - University College London/LIDO/DATA_LIDO/problematic_pupil_data/EM_task/'];
% pathdat = ['//imed/fme/kliniken/iknd/rakshit/Eigene Dateien/Dorothea/Pupil_clining/EM_raw/'];

addpath(genpath('/Users/dorothea/Dropbox/matlab/toolboxes/fieldtrip-20170119'));
addpath('/Users/dorothea/Dropbox/work/Lehre_Innsbruck/Veranstaltungen_SS22/Sem_Forschung_II/MathsX/Eyetracker_Skripte/');
pathdat = '/Users/dorothea/Dropbox/work/Lehre_Innsbruck/Veranstaltungen_SS22/Sem_Forschung_II/MathsX/Daten/';


ids = [25:34]; % sort out 2 and 3

done = [];
%% break

for n = 1:length(ids)
    for run = 1:2
        %     try
        name = [num2str(ids(n)) '/math' num2str(ids(n))  num2str(run) '.asc'];
%         eyeAdjustTrigNam_DH([pathdat name],'MSG','INPUT'); % only has to run  once, is overwriting original text file (replaces MSG with INPUT in .asc file)
        [type] = ft_filetype([pathdat name]);
        
        cfg                     = [];
        cfg.trialfun            = 'trialfunction_eyetracker_mathsx';
        cfg.dataset             = [pathdat name]; % ascii converted eyelink filename
        cfg.headerformat        = 'eyelink_asc';
        cfg.dataformat          = 'eyelink_asc';
        cfg.trialdef.eventtype  = 'MSG';
        cfg.trialdef.eventvalue = [1 2 3]; % triggers
        cfg.trialdef.prestim    = .2; % add start trigger fixation cross before the whole thing starts
        cfg.trialdef.poststim   = 3; % feedback (.5 seconds into fixation cross afterwards
        cfg.triggernumber       = 7;
        cfg.channel             = {'4'}; % channel 2 is the x-coordinate
        % channel 3 is the y-coordinate
        % channel 4 is the pupil dilation
        cfg.dataformat          = 'eyelink_asc';
        cfg.headerformat        = 'eyelink_asc';
        
        cfg                     = ft_definetrial(cfg);
        eye                     = ft_preprocessing(cfg);
        
        %% save segmented data (not interpolated, baseline corrected etc.)
        
        save([pathdat name(1:end-4) '_eyedat_feedback.mat'],'eye'); % save original raw data (eye) and those where 0s are replaced with missing values
        
        done(n) = ids(n)
        
        %     catch end
        
        keep done pathdat ids n run
    end
end
close all
