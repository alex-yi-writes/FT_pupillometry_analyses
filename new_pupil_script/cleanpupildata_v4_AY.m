%% clean pupil script

% mareike.luudwig@med.ovgu.de Jan 2022
% adjusted for Alex 04-08-2022

% consists of 4 parts: automatic cleaning(1), visual check(2), define in and out trials(3), zscore data(4)
clear; close all

root = ['/Users/yeojin/Desktop/E_data/EA_raw/EAB_pupil/MRPETpilot/ASC/test/'];
pathdat = [root 'pupil_seg/'];
% savepath =[root 'pupil_clea\'];

%% parameteres

% test 
ids = [22207];
name_stim = 're';
%% parameters for cleaning

trls = []; samples = [];
cutleft = 25; cutright = 25;        % data based decision

%%
for sub = 1:length(ids)

    load([pathdat name_stim num2str(ids(sub)) 'stim_eyedat.mat']);
%     
%     for t = 1:length(eye_rv.eyedat.trial)
%         puptrial = cell2mat(eye_rv.eyedat.trial(1,73)');
%         figure,plot(puptrial)
%     end


        data = cell2mat(eye_rv.eyedat.trial');
        samples = [eye_rv.eyedat.sampleinfo];
        trls = [eye_rv.eyedat.cfg.trl];


        for tr = 1:size(data,1)
            
            %storing trial number
            tt{tr,1} = tr; 
            
            %storing data in temp 
            temp = data(tr,:);

            %Step1: Smoothing with moving window average using hanning window 11ms: Mathot et al., 2013
            L = 15;
            win = hann(L);
            filtered = conv(temp, win, 'same');
            smoothedsignal = movmean(filtered, 15); % filtered moving winow over sliding factor k
            %figure,plot(smoothedsignal)

            %Step2: Generate a velocity profile -----------------------------------------------------------------------------
            v=diff(smoothedsignal);

            %Step3: Define individual threshold -----------------------------------------------------------------------------
            mean_t = mean(v);
            std_t = std(v);
            manualthreshold = abs(mean_t - std_t);

            %Step4: get start,stop and big blink(zero data) of v ------------------------------------------------------------
            % 1 is index 
            % 2 is the value 
            
            %start: where velocity is below -threshold
            start(:,1) = find(v <= -manualthreshold)';
            start(:,2) = v(v <= -manualthreshold)';

            %stop: where velocity is above manual threshold
            stop(:,1) = find(v>=manualthreshold)';
            stop(:,2) = v(v>=manualthreshold)';

            %big blink: where velocity is zero 
            bb(:,1) = find(v == 0)';
            bb(:,2) = v(v == 0)';


            %Step 5: Finding first & last points; cutting left & right ------------------------------------------------------
            
            %for start: 
            STRT = diff([0,diff(start(:,1)')==1,0]);
            START(:,2) = start((STRT>0),1);                  % second column  = first ascending = 1
            START(:,3) = start((STRT<0),1);                  % third column   = last ascending  = -1  % take data points between 1 and - 1
            
            for i = 1:size(START,1)      
               % --- cutleft --- %
               if  START(i,2) > cutleft %20                        
                    START(i,1) = START(i,2) - cutleft;       % first colomn   = - cutleft
               else
                    START(i,1) = START(i,2);
               end
               % --- cutright --- %
               if START(i,3) <=  2676 %2681                        
                    START(i,4) = START(i,3) + cutright;      % fourth colomn   = + cutright
               else
                    START(i,4) = START(i,3);
               end
            end

            %for stop: 
            STP = diff([0,diff(stop(:,1)')==1,0]);
            STOP(:,2) = stop((STP>0),1);                     %first
            STOP(:,3) = stop((STP<0),1);                     %last
            
            for i = 1:size(STOP,1)                           %cutleft
               if  STOP(i,2) > cutleft %20  
                    STOP(i,1) = STOP(i,2) - cutleft;   
               else
                    STOP(i,1) = STOP(i,2);
               end
               if STOP(i,3) <=   2676 %size(1:length(    size(eyeinterp.trial{1,1})-cutright)   ) %2681 
                    STOP(i,4) = STOP(i,3) + cutright;        %cutright
                else
                    STOP(i,4) = STOP(i,3);
               end
            end


            %for big blinks: 
            BB = diff([0,diff(bb(:,1)')==1,0]);
            BIGBLINK(:,2) = bb((BB>0),1);                    %first
            BIGBLINK(:,3) = bb((BB<0),1);                    %last
            
             for i = 1:size(BIGBLINK,1)                      %cutleft
               if  BIGBLINK(i,2) > cutleft %20 
                    BIGBLINK(i,1) = BIGBLINK(i,2) - cutleft;
               else
                    BIGBLINK(i,1) = BIGBLINK(i,2);
               end
               if BIGBLINK(i,3) <=  2676 %2681                      %cutright
                   BIGBLINK(i,4) = BIGBLINK(i,3) + cutright;
               else
                   BIGBLINK(i,4) = BIGBLINK(i,3);
               end
             end
             
             % get data points of small and big blinks
             if ~isempty (START)
                 for i = 1:size(START,1)
                     a{i,1} = length(START(i,1):START(i,4));
                     a{i,2} = [START(i,1), START(i,4)];
                 end
             else
                 a{i,1} = 0;
                 a{i,2} = 0;
                
             end
             %sum sample points start
             tt{tr,2} = sum(cell2mat(a(:,1)));
             
             
             if ~isempty (STOP)
                 for i = 1:size(STOP,1)
                     a{i,3}  = length(STOP(i,1):STOP(i,4));
                     a{i,4} = [STOP(i,1), STOP(i,4)];
                 end
             else
                 a{i,3} = 0;
                 a{i,4} = 0;
             end
             %sum sample points stop
             tt{tr,3} = sum(cell2mat(a(:,3)));
             
             if ~isempty (BIGBLINK)
                 for i = 1:size(BIGBLINK,1)
                     a{i,5} = length(BIGBLINK(i,1):BIGBLINK(i,4));
                     a{i,6} = [BIGBLINK(i,1), BIGBLINK(i,4)];
                 end
             else
                 a{i,5} = 0;
                 a{i,6} = 0;
             end
             %sum sample points stop
             tt{tr,4} = sum(cell2mat(a(:,5)));
             
             % find sum of all blinks
             tt{tr,5} = sum(tt{tr,2} + tt{tr,3} + tt{tr,4});
             
             
             %Step 6: Set the values to NaN; ---------------------------------------------------------------------------------
             
             %for start
             for i=1:size(START,1)
                     v(START(i,1):START(i,4))=NaN;
                     temp(START(i,1):START(i,4))=NaN;
             end
             
             %for stop
             for i=1:size(STOP,1)
                     v(STOP(i,1):STOP(i,4))=NaN;
                     temp(STOP(i,1):STOP(i,4))=NaN;
             end
             
             %for big blinks
             for i=1:size(BIGBLINK,1)
                     v(BIGBLINK(i,1):BIGBLINK(i,4))= NaN;
                     temp(BIGBLINK(i,1):BIGBLINK(i,4))=NaN;
             end

             %Interpolation & Plotting -----------------------------------------------------------------------------------------
                
             % if first signal is NaN take first var
             if isnan(temp(1))
                 var = temp(~isnan(temp));
                 temp(1) = var(1);
             end
             
             if isnan(temp(end))
                 var = temp(~isnan(temp));
                 temp(end) = var(end); %replace NaN with last legitimate val
             end
             
             
             
             % get zero trials
             if sum(isnan(temp)) > 2160
                 datai(tr,1:size(data,2)) = temp;
                 
             else

                 %cut to correct size and interpolate missing data
                 xxi = 1:size(data,2);
                 xx = 1:size(data,2);
                 zs2 = find(isnan(temp));
                 temp1 = temp;
                 temp1(zs2)=[];
                 xx(zs2)=[];
                 temp3 = interp1(xx, temp1, xxi); %recon
                 data3(tr,1:size(data,2)) = temp3(1:size(data,2));
%                 
%                  
             end
             
             figure;
             plot(data(tr,:));hold on
             plot(data3(tr,:),'r');hold off
             
             %Labelling trial as good or bad based on criteria
             if (2701-(tt{tr,5})) < 1891 % 70%    2701 = 100% data, 2160 = 80% data (percentBadCriteria * size(datai(tr,:)))
                 tt{tr,6} = 1 ;%'bad';
             else
                 tt{tr,6} = 0 ;%'good';
             end
            

            clearvars -except tt pathdat sub ids eye_rv eye_rv.eyedat trls samples name_stim data datai data3 cutleft cutright percentBadCriteria
        end

%              
%         figure,
%         for i = 1:60
%             subplot(12,5,i);
%             plot(data(i,:));hold on
%             plot(datai(i,:),'r');hold off
%         end
%  
        eye_rv.eyedat.tt = tt;
        
        % ------------------------- %
        %  make data structure new  %
        % ------------------------- %
       
        eyeinterp = eye_rv.eyedat;
       
        for i = 1:length(eyeinterp.tt)
            
            temp = data3(i,:);
            eyeinterp.trial{1,i} = temp;
            eyeinterp.time{1,i} = eye_rv.eyedat.time{1};
            
        end

        eyeinterp.cfg.trl = trls;
        eyeinterp.sampleinfo = samples;
   
%         pathsave = ['C:\Users\FriedaLubre\Documents\AlexData\pupil_clea\'];
 
%         save([pathsave name_stim num2str(ids(sub)) '_eyedat_cleaned.mat'],'eyeinterp', 'eye_rv');
     
end

%% 2) optional: visual check

% here you can use the fieldtrip version, however I had to adjust this
% because my segmented data needed to be subdivided again, so fieldtrip was
% not working with my inputs, so I was looking for an alternative, which
% works the same:

% What I do is: I check all trials again and set 1 for rejection in tt(:,7)
% What you can also do is, to just check the remaining trials after the autom
% cleaning: eye_rv.eyedat.tt(:,6) contains 1 for removed and 0 for
% remaining trials

for sub = 1:length(ids)
 
        load([pathsave name_stim num2str(ids(sub)) '_eyedat_cleaned.mat'])
   
        for ii = 1:length(eyeinterp.trial)
            plot(cell2mat(eyeinterp.trial(1,ii))); title([num2str(ii)]);
            [x,y,button] = ginput(1);
            if button == 97 % A - for accept
                indices(ii) = 0;
                disp('Input Accepted!');
            elseif button == 114 % R - for reject
                indices(ii) = 1;
                disp('Input Rejected!');
            else
                disp('Button not recongnized!')
            end
        end
        manurej = [indices'];
        eyeinterp.tt(:,7) = num2cell([indices']);
        
        clear manurej indices
        
        pathsave = ['C:\Users\FriedaLubre\Documents\AlexData\pupil_clea\'];
 
        save([pathsave name_stim num2str(ids(sub)) '_eyedat_cleaned_vis.mat'],'eyeinterp', 'eye_rv');
     
end

%% 3) define final data for in and out trials
%  4) zscore

for sub = 1:length(ids)
 
        load([pathsave name_stim num2str(ids(sub)) '_eyedat_cleaned.mat'])
        
        for i = 1:length(eyeinterp.tt)
            
            % auto == rej && manual == rej
            if cell2mat(eyeinterp.tt(i,6)) == 1 && cell2mat(eyeinterp.tt(i,7)) == 1
                eyeinterp.tt(i,8) = num2cell(1);
                
            % auto == rej     
            elseif  cell2mat(eyeinterp.tt(i,6)) == 1 && cell2mat(eyeinterp.tt(i,7)) == 0
                eyeinterp.tt(i,8) = num2cell(1);
            
            % manual == rej
            elseif  cell2mat(eyeinterp.tt(i,6)) == 0 && cell2mat(eyeinterp.tt(i,7)) == 1
                eyeinterp.tt(i,8) = num2cell(1);
            
            % no rejection
            else
                eyeinterp.tt(i,8) = num2cell(0);
            end

        end
        
        % exclude bad trials
        for i = 1:length(eyeinterp.trial)
            if cell2mat(eyeinterp.tt(i,end)) == 1
                trial_out{i,sub} =  cell2mat(eyeinterp.tt(i));
            else
                trial_in{i,sub} = cell2mat(eyeinterp.tt(i));
            end
        end
        
        % sort good trials
        Trials_In = find(~cellfun(@isempty,trial_in));
        
        for i = 1:length(Trials_In)
            P1_Adat(i) = eyeinterp.trial(Trials_In(i));
            P1_Atime(i) = eyeinterp.time(Trials_In(i));
        end
        
        P1_dat = cell2mat(P1_Adat');
        
        P1_time = P1_Atime{1}(:);
        
        
        % change format
        eval(['P1_size' '=size(P1_dat,1);' ]);
        
        eval(['P1_time' '=P1_time;' ]); %=eyeinterp.time{1}(:);'
        
        clear trial_in trial_in Trials_In P1_Adat P1_Atime P1_dat P1_time
       
         
        %Finding z-score
        eval(['P1_dat' '= zscore(P1_dat'  ',1,2);']);
        eval(['storebslnoise' '= nanmean(nanstd(P1_dat' '(:,1:200)));']); % mean of std within 200ms baseline per person, assess whether baseline noise is different
        
        
        %baseline correction (0 to 200 ms time window)
        for trl = 1:size(eval(['P1_dat']),1)
            eval(['P1_dat' '(' num2str(trl) ',:) = P1_dat' '(' num2str(trl) ',:)-nanmean(P1_dat' '(' num2str(trl) ',1:200));']);
        end
        
end

% show final data
for i = 1:length(P1_Adat)
    figure,plot(cell2mat(P1_Adat(1,73)))
end


