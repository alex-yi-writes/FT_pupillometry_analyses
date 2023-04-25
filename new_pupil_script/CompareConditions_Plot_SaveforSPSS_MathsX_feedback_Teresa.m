% compare conditions on eyetracker data, plot results and save for analyses
% with SPSS - shortened training version

clear
close all


% addpath('/Users/dorothea/Dropbox/matlab/toolboxes/fieldtrip-20170119'); % add your path to fieldtrip folderpathdat = ['\\imed\fme\kliniken\iknd\rakshit\Eigene Dateien\Dorothea\Pupil_clining\EM_raw\'];
% addpath('/Users/dorothea/Dropbox/work/Lehre_Innsbruck/Veranstaltungen_SS22/Sem_Forschung_II/FACES/Eyetracker_Skripte/');
pathdat = '/Users/dorothea/Dropbox/work/Lehre_Innsbruck/Veranstaltungen_SS22/Sem_Forschung_II/MathsX/Daten/';


pathsave = '/Users/dorothea/Dropbox/work/Lehre_Innsbruck/Veranstaltungen_SS22/Sem_Forschung_II/MathsX/Daten/figures/';mkdir(pathsave)

baseline = 1:200;
stats_window = [1000:2000]; % in which time window mean of data for analyses?

ids = [2:22 24:34]; % list of ids


for id = 1:length(ids)
    
    load([pathdat num2str(ids(id)) '/cleaned/' num2str(ids(id)) '_eye_cleaned_feedback.mat'],'eye_rv');
    
    dat = cell2mat(eye_rv.trial');
    
    trlind = find(eye_rv.trlraus == 0); % which pupil trials are still in (for matching behavioral data)
    
    % baseline correct
    bsldat = nanmean(dat(:,baseline)');
    dat = dat-repmat(bsldat',1,size(dat,2));
    storebslnoise(id) = nanmean(nanstd(dat(:,baseline)'));% mean of std within baseline per person for covariate lateron
    
    % pic definitions 
    %% correct or incorrect feedback, only when response in time 
    condnames = {'correct';'incorrect'};
    indcorr = find(eye_rv.corrresp(trlind) == 1);% correct
    indwrong = find(eye_rv.corrresp(trlind) == 0 & eye_rv.rts_task(trlind) > 0.3);% incorrect & no missing / two quick response
    dat1(id,:) = nanmean(dat(indcorr,:));
    dat2(id,:) = nanmean(dat(indwrong,:));
    storenum(id,1) = ids(id); % store trialnum per condition
    storenum(id,2) = length(indcorr);
    storenum(id,3) = length(indwrong);
         
    
    %% prepare for export to SPSS: get mean time window 1-2sec after stimulus
    windowexp = stats_window+max(baseline);
    storeexp(id,1) = nanmean(dat1(id,windowexp));
    storeexp(id,2) = nanmean(dat2(id,windowexp));
    
    % plot individual data
    subplot(7,6,id);
    time = eye_rv.time{1};
    plot(time,[dat1(id,:)' dat2(id,:)']);title(['ID ' num2str(ids(id))]);
    if ids(id) == max(ids)
        legend( condnames{1}, condnames{2});
    end
end
h = gcf;
set(h,'Position',[0 100 1400 800]);
set(gcf,'PaperPositionMode','auto')
saveas(h,[pathsave 'indivs_pupil_teresa_' condnames{1} '_' condnames{2}],'fig')


%% plot summary across ids
range = 1:length(time);
mean1 = nanmean(dat1);mean2 = nanmean(dat2);
se1 = nanstd(dat1)./sqrt(length(ids));se2 = nanstd(dat2)./sqrt(length(ids));

figure,
shadedErrorBar(time(range),mean1(range),se1(range),'or','transparent');hold on % cond1
shadedErrorBar(time(range),mean2(range),se2(range),'ogr','transparent');hold off % cond2
set(gca,'FontSize',20);
set(gcf,'PaperPositionMode','auto');ylabel('Pupil diameter (zscore)','FontSize',20);
xlabel('time after target onset (ms)','FontSize',20);title([ condnames{1} ' red, ' condnames{2} ' green']);
axis([-.2 2.5 -1 3]);
h = gcf; 
saveas(h,[pathsave 'summary_pupil_teresa_' condnames{1} '_' condnames{2}],'fig')
print(h,'-djpeg','-r300',[pathsave 'summary_pupil_teresa_' condnames{1} '_' condnames{2}]);


%% dataexport for SPSS analyses (or R)
% including variance in baseline time window, number of trials per condition
dataexport = [storenum storebslnoise' storeexp];
% close all
pathspss = [pathdat '/spss/'];mkdir(pathspss);

h = fopen([pathspss  condnames{1} '_' condnames{2} 'feedb_dataexp_' num2str(windowexp(1)-max(baseline)) '_' num2str(windowexp(end)-max(baseline)) '_Teresa.txt'],'wt');
fprintf(h,[['id\t '] ['numtr'  condnames{1} '\t '] [' numtr'  condnames{2} '\t '] 'blsnoise\t ' ['mean'  condnames{1} '\t '] ['mean'  condnames{2} '\n']]);
for k=1:size(dataexport,1)
    fprintf(h,'%4.0f\t %2.0f\t %2.0f\t %1.4f\t %1.4f\t %1.4f\n',...
        dataexport(k,1),dataexport(k,2),dataexport(k,3),dataexport(k,4),dataexport(k,5),dataexport(k,6));
end
fclose(h);
clear k h

%% prelim plots
% pick which variables to plot
d = dataexport;
range = 5:6; 
figure,barweb(nanmean(d(:,range))',nanstd(d(:,range))'./sqrt(size(d,1)));
legend(condnames{range-4});
