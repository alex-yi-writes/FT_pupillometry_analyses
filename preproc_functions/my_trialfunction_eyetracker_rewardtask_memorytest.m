function [trl,eventlist] = my_trialfunction_eyetracker_rewardtask_memorytest(cfg)

hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

raus = [];
trlnums = [event.value];
rausrows = find(ismember(trlnums,raus));
event(rausrows) = [];

trl = [];eventlist = [];

% adjust read events
for e = 1:3:length(event)
    event(e).value = 1;   % fix1
    event(e+1).value = 2; % target
    event(e+2).value = 3; % response
end

stimtype = 2; % first trigger is fixation cross, second target

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
