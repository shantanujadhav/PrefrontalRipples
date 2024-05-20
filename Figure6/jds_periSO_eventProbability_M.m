function jds_periSO_eventProbability_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates and plots the fold change occurence of ripples and spindles
%relative to slow oscillation troughs (up states)
%%------------------------------------------------------------------------
binsize = 20;
day = 1;

pret=2000; postt=2000; 
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
peakAlign = 0; %align to peak (1) or trough (0)
iricrit = 0;

% evList = {'rippletime_coordSWS','ctxrippletime_coordSWS','ctxspindletime_SWS'}; %coordinated ripples
% evList = {'rippletime_noncoordSWS','ctxrippletime_noncoordSWS','ctxspindletime_SWS'}; %noncoordinated ripples
evList = {'rippletime_leadlag_embedded','rippletime_leadlag_embedded','ctxspindletime_SWS'};
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

            if (evs == 1) && (strcmp(evList{evs},'rippletime_leadlag_embedded'))
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
            elseif (evs == 2) && (strcmp(evList{evs},'rippletime_leadlag_embedded')) %get ctx ripple events
                if CA1lead
                    if ~isempty(ripplecoupling{day}{epoch}.ctxlagriptimes)
                        currEvents = ripplecoupling{day}{epoch}.ctxlagriptimes(:,1);
                    else
                        currEvents = [];
                    end
                else
                    if ~isempty(ripplecoupling{day}{epoch}.ctxlagriptimes)
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
                        smoothvect(sum(evTmp,1)./length(evTmp(:,1)),g1)];
                    dataCompileSepBck1 = [dataCompileSepBck1;...
                        mean(sum(evTmpBck,1)./length(evTmpBck(:,1)))];
                end
            elseif evs == 2
                dataCompile2 = [dataCompile2; evTmp];
                dataCompileBck2= [dataCompileBck2; evTmpBck];
                dataCompileTime2 = [dataCompileTime2; allEvs];
                epTime2 = [epTime2; nanmean(evTimeEp)];
                if length(currEvents) > 10
                    dataCompileSep2 = [dataCompileSep2; ...
                        smoothvect(sum(evTmp,1)./length(evTmp(:,1)),g1)];
                    dataCompileSepBck2 = [dataCompileSepBck2;...
                        mean(sum(evTmpBck,1)./length(evTmpBck(:,1)))];
                end
            elseif evs == 3
                dataCompile3 = [dataCompile3; evTmp];
                dataCompileBck3 = [dataCompileBck3; evTmpBck];
                dataCompileTime3 = [dataCompileTime3; allEvs];
                epTime3 = [epTime3; nanmean(evTimeEp)];
                if length(currEvents) > 10
                    dataCompileSep3 = [dataCompileSep3; ...
                        smoothvect(sum(evTmp,1)./length(evTmp(:,1)),g1)];
                    dataCompileSepBck3 = [dataCompileSepBck3;...
                        mean(sum(evTmpBck,1)./length(evTmpBck(:,1)))];
                end
            end
        end
    end
end
%Plot the fold change occurence
figure; hold on
dataTmp = dataCompile1;
dataTmpBck = dataCompileBck1;
meanProb = mean(sum(dataTmpBck)./length(dataTmpBck(:,1))); %mean prob in back window over all bins
dataProb = sum(dataTmp)./length(dataTmp(:,1)); %probability in each bin
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

figure
window = 250;
%Xms window around trig
windowData1 = dataCompileTime1(find(dataCompileTime1>-window & dataCompileTime1 < window))/1000;
windowData2 = dataCompileTime2(find(dataCompileTime2>-window & dataCompileTime2 < window))/1000;
windowData3 = dataCompileTime3(find(dataCompileTime3>-window & dataCompileTime3 < window))/1000;
datacombinedEvents = [windowData1; windowData2; windowData3];
g1 = repmat({'Ev1'},length(windowData1),1);
g2 = repmat({'Ev2'},length(windowData2),1);
g3 = repmat({'Ev3'},length(windowData3),1);
g = [g1;g2;g3];

%Boxplot of timings
figure;
h = boxplot(datacombinedEvents,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['Timing of events relative to trigger'])
set(gcf, 'renderer', 'painters')

%CDF
figure; hold on
plotcdf(windowData1);
plotcdf(windowData2);
plotcdf(windowData3);
set(gcf, 'renderer', 'painters')

figure; hold on 
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
figure; hold on
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
