function [trl,eventlist] = trialfunction_eyetracker_mathsx(cfg)

hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

trl = [];eventlist = [];
% adjust read events, 3 events, fixation, target and question, we want
% target
for e = 1:9:length(event)
    event(e).value = 1; %cue 
    event(e+1).value = 2;% fixation
    event(e+2).value = 3;% task
    event(e+3).value = 4; % solution
    event(e+4).value = 5; % fixation
    event(e+5).value = 6; % difficulty question
    event(e+6).value = 7; % fixation 
    event(e+7).value = 8; % feedback
    event(e+8).value = 9; % fixation (ISI)
end

stimtype = cfg.triggernumber; % first trigger is fixation cross, second target

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
