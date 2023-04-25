function [trl,eventlist] = my_trialfunction_eyetracker_continuous(cfg)

hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% get time of begin of task (MSG -100, which is start of first volume and
% start segmenting 1sec after that, also mark end of recording (MSG - 200)

trl = [];eventlist = [];

last_trigger = event(end).sample;
first_trigger = event(1).sample;
segmentduration = cfg.trialdef.poststim*hdr.Fs;

% how many volumes do we get in?
duration_total = last_trigger - first_trigger;
duration_total = duration_total - 1000; % one second offset
num_volumes = fix(duration_total/segmentduration);

% check that first trigger is MSG = 100 (value) 0
if event(1).value > 0
    disp('Error: trigger for beginning of recording is missing')
elseif event(1).value == 0 % (MSG = 100)
    % now segment time windows according to volume duraiton, but shift by 1 sec
    for c=1:num_volumes-1 % omit last volume
        % define beginning, end, and offset of trial
        begsample = (first_trigger + (c-1)*segmentduration + 1000) - cfg.trialdef.prestim*hdr.Fs; % beginning of trial
        endsample = (first_trigger + c*segmentduration + 1000); % end of trial
        offset = -cfg.trialdef.prestim*hdr.Fs; % offset (time before stimulus onset)
        % store individual sample information of every trial
        trl(end+1, :) = round([begsample endsample offset]);
        eventlist{end+1} = 1;
    end
end