% clean eyetracker data

clear
close all

pathdat = ['/Users/dorothea/Dropbox/eyetracker_data/2opt_fmri/'];
mkdir([pathdat 'dataclean/']);

% settings for cleaning of data

TR = 3.04;
rejectlatency = [-.2 TR];
suff = 3;% % 2d or 3d

load('/Users/dorothea/Documents/Mitarbeiter/Mareike_Mastearbeit/Skripte/ids.mat');

% not all of these will have pupil data
rangeids = ids;

for n = 1:length(ids)
    for p = 1:6
        %         try
        name = [num2str(ids(n)) '_' num2str(p) num2str(suff) '.asc']
        load([pathdat name(1:end-4) '_eyedat_cont_run' num2str(p) '.mat']);
               
        %% reject visual
        cfg.method = 'channel'; % first overview
        cfg.metric = 'zvalue';
        cfg.latency = rejectlatency;
        eyeclean1 = ft_rejectvisual(cfg, eyenew)
        
        
        cfg.method = 'trial'; % then individual trials
        cfg.metric = 'zvalue';
        cfg.latency = rejectlatency;
        eye_rv = ft_rejectvisual(cfg, eyeclean1)
        
        % determine which trials out
        trlnew = eye_rv.sampleinfo(:,1);
        trlold = eye_rv.cfg.previous.previous.trl(:,1);
        trl_raus = ~ismember(trlold,trlnew);
        sum(trl_raus)
        
        save([pathdat 'dataclean/' name(1:end-4) '_eyedat_cont_run' num2str(p) '_eye_rv200.mat'],'eye_rv','blinkbig','blinksmall');
        keep cutleft cutright pathdat cuttimes fenster smallblinksize n ids rejectlatency rangeids suff
        %         catch end
    end
    
    close all
end