function [trl,eventlist] = my_trialfunction_eyetracker_artefacts(cfg)

hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

raus = [];

trlnums = [event.value];
rausrows = find(ismember(trlnums,raus));
event(rausrows) = [];

trl = [];eventlist = [];
% adjust read events
for e = 1:4:length(event)
    event(e).value = 1;   % fix1
    event(e+1).value = 2; % inst
    event(e+2).value = 3; % fix2
    event(e+3).value = 4; % stim
end

stimtype = 4; % first trigger is fixation cross, second target

count = 0;
for i=1:length(event)
    if strcmp(event(i).type, 'INPUT')
        % it is a trigger, see whether it has the right value
        if ismember(event(i).value, stimtype)
            count = count +1;
            % define beginning, end, and offset of trial
            begsample = event(i).sample - cfg.trialdef.prestim*hdr.Fs; % beginning of trial
            endsample = event(i).sample + cfg.trialdef.poststim*hdr.Fs; % end of trial
            offset = -cfg.trialdef.prestim*hdr.Fs; % offset (time before stimulus onset)
            % store individual sample information of every trial
            trl(end+1, :) = round([begsample endsample offset]);
            eventlist{end+1} = event(i).value;
        end
    end
end
