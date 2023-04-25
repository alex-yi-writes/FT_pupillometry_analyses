% function outfile = eyeAdjustTrigNam_DH(infile,orgTrigNam,newTrigNam)
%
% triggers in eyelink file are indicated by 'MSG'. However, fieldtrip
% requires it to be called 'INPUT'. This function adjusts these prefixes
% and writes a new file which should be readable by fieldtrip.
%
% @infile: file of original eyelink ascii-file
% @orgTrigNam: prefix of original triggers (default: 'MSG')
% @newTrigNam: prefix of new triggers (default: 'INPUT': readable by
% fieldtrip)
% @outfile: name of output-file
%
% TH, 10.14
%
function outfile = eyeAdjustTrigNam_DH(infile,orgTrigNam,newTrigNam)

if nargin < 3; newTrigNam = 'INPUT'; end % new
if nargin < 2; orgTrigNam = 'MSG'; end % old

% read file
fid=fopen(infile,'r','n','US-ASCII');
C = textscan(fid,'%s','delimiter','\n');
fclose(fid);

startstring = ['0  100']; % leave those in
endstring = ['0  200']; %

% first replace everything with MSG
cc = 0;
for i = 1:length(C{1})
    [in,o] = regexp(C{1}{i},newTrigNam);
    if ~isempty(in)
        C{1}{i} = [orgTrigNam C{1}{i}(o+1:end)];
        cc = cc+1;
    end
end
% then replace orgTrigNam (start at line 30)
cc = 0; gonow = 0;
for i = 30:length(C{1})
    [in,o] = regexp(C{1}{i},orgTrigNam);
    if ~isempty(in)
        if strfind(C{1}{i},startstring) > 0
            gonow = 1;
        elseif strfind(C{1}{i},endstring) > 0
        elseif ~isempty(in) && gonow ==1
            C{1}{i} = [newTrigNam C{1}{i}(o+1:end)];%% addded a tab between INPUT and samplepoint
            cc = cc+1;
        end
    end
end


% outfile

[path,name,ext] = fileparts(infile);
outfile = [path '/' name  ext];

% write to file
fid=fopen(outfile,'wt','n','US-ASCII');
for i = 1:length(C{1})
    fprintf(fid,C{1}{i});
    fprintf(fid,'\n');
end
fclose(fid)
end