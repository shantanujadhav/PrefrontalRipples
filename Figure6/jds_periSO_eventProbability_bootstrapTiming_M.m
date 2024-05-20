function jds_periSO_eventProbability_bootstrapTiming_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates and plots the fold change occurence of ripples and spindles
%relative to slow oscillation troughs (up states).
%Also, resampling of occurence matrix and comparison of peak bin to 0 with a
%ttest to determine significance of peaks of fold change occurence result.
%Timing of events relative to SO.
%%------------------------------------------------------------------------
binsize = 20;
day = 1;

pret=2000; postt=2000; %% Larger Window
peakbins = find(abs(-pret:binsize:postt)<=250);
bwin = [-3000 -2000];
g1 = gaussian(3, 3);

dataCompileSep1 = [];
dataCompileSep2 = [];
dataCompileSep3 = [];

dataCompileSepBck1 = [];
dataCompileSepBck2 = [];
dataCompileSepBck3 = [];

dataCompile1 = [];
dataCompile2 = [];
dataCompile3 = [];

dataCompileTime1 = [];
dataCompileTime2 = [];
dataCompileTime3 = [];

dataCompileBck1 = [];
dataCompileBck2 = [];
dataCompileBck3 = [];

epTime1 = [];
epTime2 = [];
epTime3 = [];

alignTo = 'SO'; %SO or delta
peakAlign = 0;
iricrit = 0;

evList = {'rippletime_leadlag','rippletime_leadlag','ctxspindletime_SWS'};
CA1lead = 1;
totalTrigs = 0;
epochs = 1:2:17;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    load(sprintf('%s%sslowoscdeltatimesSep_SWS%02d.mat', dir,animalprefix,day)) %load delta and so
    
    for evs = 1:length(evList)
        for ep = 1:length(epochs)
            
            epoch = epochs(ep);
            evTimeEp = [];
            if peakAlign == 1
                if strcmp(alignTo,'SO')
                    alignTimes = slowosc{day}{epoch}.SOpeaktime;
                else
                    alignTimes = slowosc{day}{epoch}.Deltapeaktime;
                end
            else
                if strcmp(alignTo,'SO')
                    alignTimes = slowosc{day}{epoch}.SOtroughtime;
                else
                    alignTimes = slowosc{day}{epoch}.Deltatroughtime;
                end
            end

            if (evs == 1) && (strcmp(evList{evs},'rippletime_leadlag'))
                load(sprintf('%s%s%s%02d.mat',dir,animalprefix,evList{evs},day));
                if CA1lead
                    if ~isempty(ripplecoupling{day}{epoch}.hpleadriptimes)
                        currEvents = ripplecoupling{day}{epoch}.hpleadriptimes(:,1);
                    else
                        currEvents = [];
                    end
                else
                    if ~isempty(ripplecoupling{day}{epoch}.hplagriptimes)
                        currEvents = ripplecoupling{day}{epoch}.hplagriptimes(:,1);
                    else
                        currEvents = [];
                    end
                end
            elseif (evs == 2) && (strcmp(evList{evs},'rippletime_leadlag')) %get ctx ripple events
                if CA1lead
                    if ~isempty(ripplecoupling{day}{epoch}.ctxlagriptimes)
                        currEvents = ripplecoupling{day}{epoch}.ctxlagriptimes(:,1);
                    else
                        currEvents = [];
                    end
                else
                    if ~isempty(ripplecoupling{day}{epoch}.ctxleadriptimes)
                        currEvents = ripplecoupling{day}{epoch}.ctxleadriptimes(:,1);
                    else
                        currEvents = [];
                    end
                end
            else
                load(sprintf('%s%s%s%02d.mat',dir,animalprefix,evList{evs},day)); %load ripple, ctxripple, spindle
                if evs == 1
                    currEvents = ripple{day}{epoch}.starttime;
                elseif evs == 2
                    currEvents = ctxripple{day}{epoch}.starttime;
                elseif evs == 3
                    currEvents = ctxspindle{day}{epoch}.starttime;
                end
            end
            totalTrigs = totalTrigs + length(alignTimes);
            evTmp = [];
            evTmpBck = [];
            alignTimes = alignTimes*1000; %to ms
            if isempty(alignTimes)
                continue
            end
            if iricrit == 1
                iri = diff(alignTimes);
                keepidx = [1;find(iri>=1000)+1];
                alignTimes = alignTimes(keepidx);
            end
            if isempty(currEvents)
                continue
            end
            currEvents = currEvents*1000;
            allEvs = [];
            for i=1:length(alignTimes)
                % Align to event trough or peak
                % ------------------------------------
                currTrig = alignTimes(i);

                currEvs =  currEvents(find((currEvents>=(currTrig-pret)) & (currEvents<=(currTrig+postt))));
                currEvs = currEvs-(currTrig); %time relative to trigger
                allEvs = [allEvs; currEvs];
                evTimeEp = [evTimeEp; currEvs];
                if isempty(currEvs)
                    histEvs = zeros(1,length(-pret:binsize:postt));
                    histEvs = histEvs(:).';
                else
                    histEvs = histc(currEvs,-pret:binsize:postt);
                    histEvs = histEvs(:).';
                end
                evTmp = [evTmp; histEvs];

                %background window
                currEvsBck =  currEvents(find((currEvents>=(currTrig+bwin(1))) & (currEvents<=(currTrig+bwin(2)))));
                currEvsBck = currEvsBck-(currTrig); %time relative to trigger
                if isempty(currEvsBck)
                    histEvsBck = zeros(1,length(bwin(1):binsize:bwin(2)));
                    histEvsBck = histEvsBck(:).';
                else
                    histEvsBck = histc(currEvsBck,bwin(1):binsize:bwin(2));
                    histEvsBck = histEvsBck(:).';
                end
                evTmpBck = [evTmpBck; histEvsBck];

            end
            if evs == 1
                dataCompile1 = [dataCompile1; evTmp];
                dataCompileBck1 = [dataCompileBck1; evTmpBck];
                dataCompileTime1 = [dataCompileTime1; allEvs];
                epTime1 = [epTime1; nanmean(evTimeEp)];
                if length(currEvents) > 10
                    dataCompileSep1 = [dataCompileSep1; ...
                        smoothvect(sum(evTmp)./length(evTmp(:,1)),g1)];
                    dataCompileSepBck1 = [dataCompileSepBck1;...
                        mean(sum(evTmpBck)./length(evTmpBck(:,1)))];
                end
            elseif evs == 2
                dataCompile2 = [dataCompile2; evTmp];
                dataCompileBck2= [dataCompileBck2; evTmpBck];
                dataCompileTime2 = [dataCompileTime2; allEvs];
                epTime2 = [epTime2; nanmean(evTimeEp)];
                if length(currEvents) > 10
                    dataCompileSep2 = [dataCompileSep2; ...
                        smoothvect(sum(evTmp)./length(evTmp(:,1)),g1)];
                    dataCompileSepBck2 = [dataCompileSepBck2;...
                        mean(sum(evTmpBck)./length(evTmpBck(:,1)))];
                end
            elseif evs == 3
                dataCompile3 = [dataCompile3; evTmp];
                dataCompileBck3 = [dataCompileBck3; evTmpBck];
                dataCompileTime3 = [dataCompileTime3; allEvs];
                epTime3 = [epTime3; nanmean(evTimeEp)];
                if length(currEvents) > 10
                    dataCompileSep3 = [dataCompileSep3; ...
                        smoothvect(sum(evTmp)./length(evTmp(:,1)),g1)];
                    dataCompileSepBck3 = [dataCompileSepBck3;...
                        mean(sum(evTmpBck)./length(evTmpBck(:,1)))];
                end
            end
        end
    end
end
%% Resampling
numsamps = 100;
window = 250;
numruns = 100;
timingDist1 = [];
timingDist2 = [];
timingDist3 = [];
resamp1 = [];
resamp2 = [];
resamp3 = [];
for r = 1:numruns
    randidx1 = randperm(length(dataCompile1(:,1)));
    randidx2 = randperm(length(dataCompile2(:,1)));
    randidx3 = randperm(length(dataCompile3(:,1)));
    dataCompile1perm = dataCompile1(randidx1,:);
    dataCompile2perm = dataCompile2(randidx2,:);
    dataCompile3perm = dataCompile3(randidx3,:);
    dataBck1perm = dataCompileBck1(randidx1,:);
    dataBck2perm = dataCompileBck2(randidx2,:);
    dataBck3perm = dataCompileBck3(randidx3,:);
    cnt = 1;
    cnt2 = 100;
    for i = 1
        %fold change
        probdata1 = (nansum(dataCompile1perm(cnt:cnt2,:)))/numsamps;
        probdata2 = (nansum(dataCompile2perm(cnt:cnt2,:)))/numsamps;
        probdata3 = (nansum(dataCompile3perm(cnt:cnt2,:)))/numsamps;

        probbck1 = mean((nansum(dataBck1perm(cnt:cnt2,:)))/numsamps);
        probbck2 = mean((nansum(dataBck2perm(cnt:cnt2,:)))/numsamps);
        probbck3 = mean((nansum(dataBck3perm(cnt:cnt2,:)))/numsamps);

        tmpSamps1 = (probdata1 - probbck1)./probbck1;
        tmpSamps1 = smoothvect(tmpSamps1,g1);

        tmpSamps2 = (probdata2 - probbck2)./probbck2;
        tmpSamps2 = smoothvect(tmpSamps2,g1);

        tmpSamps3 = (probdata3 - probbck3)./probbck3;
        tmpSamps3 = smoothvect(tmpSamps3,g1);

        [M I] = max(tmpSamps1(peakbins)); %look at bins around SO trough
        timingDist1 = [timingDist1; I];
        resamp1 = [resamp1; tmpSamps1];
        [M I] = max(tmpSamps2(peakbins));
        timingDist2 = [timingDist2; I];
        resamp2 = [resamp2; tmpSamps2];
        [M I] = max(tmpSamps3(peakbins));
        timingDist3 = [timingDist3; I];
        resamp3 = [resamp3; tmpSamps3];
        cnt = cnt + numsamps;
        cnt2 = cnt2 + numsamps;
    end
end
peakbinsresamp1 = resamp1(:,peakbins);
peakbinsresamp2 = resamp2(:,peakbins);
peakbinsresamp3 = resamp3(:,peakbins);
plot(mean(peakbinsresamp1))
hold on
plot(mean(peakbinsresamp2))
plot(mean(peakbinsresamp3))

figure; hold on %plot bounded line
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(-250:binsize:250-1,nanmean(peakbinsresamp1),'-k','LineWidth',1)
boundedline(-250:binsize:250-1,nanmean(peakbinsresamp1),nanstd(peakbinsresamp1)./sqrt(length(peakbinsresamp1(:,1))),'-k');
pl2 = plot(-250:binsize:250-1,nanmean(peakbinsresamp2),'-r','LineWidth',1)
boundedline(-250:binsize:250-1,nanmean(peakbinsresamp2),nanstd(peakbinsresamp2)./sqrt(length(peakbinsresamp2(:,1))),'-r');
pl3 = plot(-250:binsize:250-1,nanmean(peakbinsresamp3),'-b','LineWidth',1)
boundedline(-250:binsize:250-1,nanmean(peakbinsresamp3),nanstd(peakbinsresamp3)./sqrt(length(peakbinsresamp3(:,1))),'-b');
set(gcf, 'renderer', 'painters')
ylabel('Probability (Fold change)')
if peakAlign == 0
    title([alignTo '-' 'trough'])
    xlabel('Time from SO trough (ms)')
else
    title([alignTo '-' 'peak'])
    xlabel('Time from SO peak (ms)')
end
xlim([-250 250])
[M1 I1] = max(sum(peakbinsresamp1))
[M2 I2] = max(sum(peakbinsresamp2))
[M3 I3] = max(sum(peakbinsresamp3))

[H1,P1] = ttest(peakbinsresamp1(:,I1))
[H2,P2] = ttest(peakbinsresamp2(:,I2))
[H3,P3] = ttest(peakbinsresamp3(:,I3))

figure; hold on %plot bounded line
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(-pret:binsize:postt,nanmean(resamp1),'-k','LineWidth',1)
boundedline(-pret:binsize:postt,nanmean(resamp1),nanstd(resamp1)./sqrt(length(resamp1(:,1))),'-k');
pl2 = plot(-pret:binsize:postt,nanmean(resamp2),'-r','LineWidth',1)
boundedline(-pret:binsize:postt,nanmean(resamp2),nanstd(resamp2)./sqrt(length(resamp2(:,1))),'-r');
pl3 = plot(-pret:binsize:postt,nanmean(resamp3),'-b','LineWidth',1)
boundedline(-pret:binsize:postt,nanmean(resamp3),nanstd(resamp3)./sqrt(length(resamp3(:,1))),'-b');
xlim([peakbins(1) peakbins(end)])
set(gcf, 'renderer', 'painters')
ylabel('Probability (Fold change)')
if peakAlign == 0
    title([alignTo '-' 'trough'])
    xlabel('Time from SO trough (ms)')
else
    title([alignTo '-' 'peak'])
    xlabel('Time from SO peak (ms)')
end

figure; hold on
plotcdf(timingDist1)
plotcdf(timingDist2)
plotcdf(timingDist3)
%%

figure; hold on

plot(-pret:binsize:postt,smoothvect(sum(dataCompile1)./length(dataCompile1(:,1)),g1),'LineWidth',2)
plot(-pret:binsize:postt,smoothvect(sum(dataCompile2)./length(dataCompile2(:,1)),g1),'LineWidth',2)
plot(-pret:binsize:postt,smoothvect(sum(dataCompile3)./length(dataCompile3(:,1)),g1),'LineWidth',2)

ylabel('Event probability')
xlabel('Time from trigger onset')
if peakAlign == 0
    title([alignTo '-' 'trough'])
else
    title([alignTo '-' 'peak'])
end
xlim([-1000 1000])
set(gcf, 'renderer', 'painters')

%% Fold change all data
figure; hold on

dataTmp = dataCompile1;
dataTmpBck = dataCompileBck1;
meanProb = mean(sum(dataTmpBck)./length(dataTmpBck(:,1))); %mean in back window
dataProb = sum(dataTmp)./length(dataTmp(:,1));
foldChange = (dataProb - meanProb)./meanProb;
plot(-pret:binsize:postt,smoothvect(foldChange,g1),'LineWidth',2)

dataTmp = dataCompile2;
dataTmpBck = dataCompileBck2;
meanProb = mean(sum(dataTmpBck)./length(dataTmpBck(:,1))); %mean in back window
dataProb = sum(dataTmp)./length(dataTmp(:,1));
foldChange = (dataProb - meanProb)./meanProb;
plot(-pret:binsize:postt,smoothvect(foldChange,g1),'LineWidth',2)

dataTmp = dataCompile3;
dataTmpBck = dataCompileBck3;
meanProb = mean(sum(dataTmpBck)./length(dataTmpBck(:,1))); %mean in back window
dataProb = sum(dataTmp)./length(dataTmp(:,1));
foldChange = (dataProb - meanProb)./meanProb;
plot(-pret:binsize:postt,smoothvect(foldChange,g1),'LineWidth',2)
ylabel('Event probability (fold change)')
xlabel('Time from trigger onset')
if peakAlign == 0
    title([alignTo '-' 'trough'])
else
    title([alignTo '-' 'peak'])
end
xlim([-1000 1000])
set(gcf, 'renderer', 'painters')

% legend(evList)
%% All event timings
figure
%Xms window around trig
windowData1 = dataCompileTime1(find(dataCompileTime1>-window & dataCompileTime1 < window))/1000;
windowData2 = dataCompileTime2(find(dataCompileTime2>-window & dataCompileTime2 < window))/1000;
windowData3 = dataCompileTime3(find(dataCompileTime3>-window & dataCompileTime3 < window))/1000;
datacombinedEvents = [windowData1; windowData2; windowData3];
g1 = repmat({'Ev1'},length(windowData1),1);
g2 = repmat({'Ev2'},length(windowData2),1);
g3 = repmat({'Ev3'},length(windowData3),1);
g = [g1;g2;g3];

figure;
h = boxplot(datacombinedEvents,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['Timing of events relative to trigger'])
set(gcf, 'renderer', 'painters')

figure
hold on
plotcdf(windowData1);
plotcdf(windowData2);
plotcdf(windowData3);
set(gcf, 'renderer', 'painters')

%% Epoch separated probs
figure; hold on %plot bounded line
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(-pret:binsize:postt,mean(dataCompileSep1),'-k','LineWidth',1)
boundedline(-pret:binsize:postt,mean(dataCompileSep1),std(dataCompileSep1)./sqrt(length(dataCompileSep1(:,1))),'-k');
pl2 = plot(-pret:binsize:postt,mean(dataCompileSep2),'-r','LineWidth',1)
boundedline(-pret:binsize:postt,mean(dataCompileSep2),std(dataCompileSep2)./sqrt(length(dataCompileSep2(:,1))),'-r');
pl3 = plot(-pret:binsize:postt,mean(dataCompileSep3),'-b','LineWidth',1)
boundedline(-pret:binsize:postt,mean(dataCompileSep3),std(dataCompileSep3)./sqrt(length(dataCompileSep3(:,1))),'-b');
xlim([-1000 1000])
set(gcf, 'renderer', 'painters')

%% Fold change timing separated by epochs
maxidx1 = [];
maxidx2 = [];
maxidx3 = [];
foldChng1 = (dataCompileSep1-dataCompileSepBck1)./dataCompileSepBck1;
for i = 1:length(foldChng1(:,1))
    tmp = foldChng1(i,:);
    [M I] = max(tmp(peakbins))
    maxidx1 = [maxidx1; I];
end
foldChng2 = (dataCompileSep2-dataCompileSepBck2)./dataCompileSepBck2;
for i = 1:length(foldChng2(:,1))
    tmp = foldChng2(i,:);
    [M I] = max(tmp(peakbins))
    maxidx2 = [maxidx2; I];
end
foldChng3 = (dataCompileSep3-dataCompileSepBck3)./dataCompileSepBck3;
foldChng3 = foldChng3(2:end,:);
for i = 1:length(foldChng3(:,1))
    tmp = foldChng3(i,:);
    [M I] = max(tmp(peakbins))
    maxidx3 = [maxidx3; I];
end

figure
plotcdf(maxidx1);
hold on
plotcdf(maxidx2);
plotcdf(maxidx3);
set(gcf, 'renderer', 'painters')

figure; hold on %plot bounded line
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(-pret:binsize:postt,nanmean(foldChng1),'-k','LineWidth',1)
boundedline(-pret:binsize:postt,nanmean(foldChng1),nanstd(foldChng1)./sqrt(length(foldChng1(:,1))),'-k');
pl2 = plot(-pret:binsize:postt,nanmean(foldChng2),'-r','LineWidth',1)
boundedline(-pret:binsize:postt,nanmean(foldChng2),nanstd(foldChng2)./sqrt(length(foldChng2(:,1))),'-r');
pl3 = plot(-pret:binsize:postt,nanmean(foldChng3),'-b','LineWidth',1)
boundedline(-pret:binsize:postt,nanmean(foldChng3),nanstd(foldChng3)./sqrt(length(foldChng3(:,1))),'-b');
xlim([-250 250])
set(gcf, 'renderer', 'painters')
ylabel('Probability (Fold change)')
if peakAlign == 0
    title([alignTo '-' 'trough'])
    xlabel('Time from SO trough (ms)')
else
    title([alignTo '-' 'peak'])
    xlabel('Time from SO peak (ms)')
end
