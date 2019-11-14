function [trl,eventlist] = my_trialfunction_eyetracker_OB_3(cfg)

hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

%take out 61,70,122,171 trials (only fixation cross : ID 1001) --> adjust this later!
% raus = [81,162];
% raus = [40 93 100 110 131 136 148 160 168 173 182 198 202 216 225 228 236 238 243 251 255];
raus = [];

trlnums = [event.value];
rausrows = find(ismember(trlnums,raus));
event(rausrows) = [];

trl = [];eventlist = [];
% adjust read events
for e = 1:2:length(event)
    event(e).value = 1;   % stim
    event(e+1).value = 2; % resp
%     event(e+2).value = 3; % fix1
%     event(e+3).value = 4; % fix2
%     event(e+4).value = 5; % feedback
end

stimtype = 1; % first trigger is stim, second response

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
