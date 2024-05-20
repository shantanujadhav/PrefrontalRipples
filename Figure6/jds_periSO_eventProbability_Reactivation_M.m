function jds_periSO_eventProbability_Reactivation_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plots CA1 and PFC assembly reactivation aligned to SO troughs and
%calculates the difference in peak timing distributions
%%------------------------------------------------------------------------
binsize = 20;
day = 1;
reactbins = 100;

pret=2000; postt=2000; 
peakbins = find(-pret:binsize:postt<=100 & -pret:binsize:postt>=-100);
rwin = [-2000 2000];
bwin = [-3000 -2000];
nstd = round(60/binsize);
g1 = gaussian(nstd, nstd*5);

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
CA1react = [];
PFCreact = [];
CA1reactShuf = [];
PFCreactShuf = [];

alignTo = 'SO'; %SO or delta
peakAlign = 0; %align to peak (1) or trough (0)

evList = {'rippletime_leadlag','rippletime_leadlag','ctxspindletime_SWS'};
CA1lead = 0; %analyze leading events?
totalTrigs = 0;
epochs = 3:2:17;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    load(sprintf('%s%srippletime_leadlag%02d.mat', dir,animalprefix,day)) %load ripples
    load(sprintf('%s%sslowoscdeltatimesSep_SWS%02d.mat', dir,animalprefix,day)) %load delta and so
    reactHp = load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    reactCtx = load(sprintf('%s%sPFC_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    
    for evs = 1:length(evList)
        for ep = 1:length(epochs)
            epoch = epochs(ep);
            assemblytmpHp = reactHp.RtimeStrength{epoch}.reactivationStrength;
            assemblytmpCtx = reactCtx.RtimeStrength{epoch}.reactivationStrength;
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

            if isempty(currEvents)
                continue
            end
            if CA1lead == 1
                ripCoord = ripplecoupling{day}{epoch}.hpleadriptimes;
            else
                ripCoord = ripplecoupling{day}{epoch}.hplagriptimes;
            end
            currEvents = currEvents*1000;
            allEvs = [];
            rTvec = assemblytmpHp{1}(:,1)*1000;
            if ~isempty(ripCoord)
                ripCoord = ripCoord(:,1)*1000;
                react_idx = [];
                for t = 1:length(alignTimes)
                    tDiff = abs(alignTimes(t) - ripCoord);
                    if ~isempty(find(ripCoord > (alignTimes(t)-250) & ripCoord < (alignTimes(t)+250)))
                        idxtmp = lookup(alignTimes(t), rTvec);
                        react_idx = [react_idx; idxtmp];
                    end
                end
            end

            if evs == 1 %only align reactivation strength once
                %ALIGN BOTH CA1 and PFC
                for c = 1:length(assemblytmpHp)
                    CA1reactTmp = [];
                    CA1reactTmpShuf = [];
                    strengthstmp = assemblytmpHp{c}(:,2);
                    strengthstmpShuf = circshift(strengthstmp, randi(length(strengthstmp)));
                    for r = 1:length(react_idx)
                        if ((react_idx(r) + reactbins) < length(strengthstmp)) && ((react_idx(r) - reactbins) > 1)
                        tmp = strengthstmp((react_idx(r) - reactbins):(react_idx(r) + reactbins)); %get vector of reactivation strenths for specified time period
                        CA1reactTmp = [CA1reactTmp; tmp'];
                        end
                    end
                    CA1react = [CA1react; smoothvect(zscore(mean(CA1reactTmp,1)),g1)];
                    for s = 1
                        strengthstmpShuf = circshift(strengthstmp, randi(length(strengthstmp)));
                        for r = 1:length(react_idx)
                            if ((react_idx(r) + reactbins) < length(strengthstmpShuf)) && ((react_idx(r) - reactbins) > 1)
                                tmp = strengthstmpShuf((react_idx(r) - reactbins):(react_idx(r) + reactbins)); %get vector of reactivation strenths for specified time period
                                CA1reactTmpShuf = [CA1reactTmpShuf; tmp'];
                            end
                        end
                        CA1reactShuf = [CA1reactShuf; smoothvect(zscore(mean(CA1reactTmpShuf,1)),g1)];
                    end
                end
                
                for d = 1:length(assemblytmpCtx)
                    PFCreactTmp = [];
                    strengthstmp = assemblytmpCtx{d}(:,2);
                    PFCreactTmpShuf = [];
                    for r = 1:length(react_idx)
                        if ((react_idx(r) + reactbins) < length(strengthstmp)) && ((react_idx(r) - reactbins) > 1)
                        tmp = strengthstmp((react_idx(r) - reactbins):(react_idx(r) + reactbins)); %get vector of reactivation strenths for specified time period
                        PFCreactTmp = [PFCreactTmp; tmp'];
                        end
                    end
                    PFCreact = [PFCreact; smoothvect(zscore(mean(PFCreactTmp,1)),g1)];
                    for s = 1
                        strengthstmpShuf = circshift(strengthstmp, randi(length(strengthstmp)));
                        for r = 1:length(react_idx)
                            if ((react_idx(r) + reactbins) < length(strengthstmpShuf)) && ((react_idx(r) - reactbins) > 1)
                                tmp = strengthstmpShuf((react_idx(r) - reactbins):(react_idx(r) + reactbins)); %get vector of reactivation strenths for specified time period
                                PFCreactTmpShuf = [PFCreactTmpShuf; tmp'];
                            end
                        end
                        PFCreactShuf = [PFCreactShuf; smoothvect(zscore(mean(PFCreactTmpShuf,1)),g1)];
                    end
                end
            end
            for i=1:length(alignTimes)
                currTrig = alignTimes(i);
                % Align to event trough or peak
                % ------------------------------------
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

figure; hold on
yyaxis left
dataTmp = dataCompile1;
dataTmpBck = dataCompileBck1;
meanProb = mean(sum(dataTmpBck)./length(dataTmpBck(:,1))); %mean in back window
dataProb = sum(dataTmp)./length(dataTmp(:,1));
foldChange1 = (dataProb - meanProb)./meanProb;
plot(-pret:binsize:postt,smoothvect(foldChange1,g1),'LineWidth',2)

dataTmp = dataCompile2;
dataTmpBck = dataCompileBck2;
meanProb = mean(sum(dataTmpBck)./length(dataTmpBck(:,1))); 
dataProb = sum(dataTmp)./length(dataTmp(:,1));
foldChange2 = (dataProb - meanProb)./meanProb;
plot(-pret:binsize:postt,smoothvect(foldChange2,g1),'LineWidth',2)

dataTmp = dataCompile3;
dataTmpBck = dataCompileBck3;
meanProb = mean(sum(dataTmpBck)./length(dataTmpBck(:,1)));
dataProb = sum(dataTmp)./length(dataTmp(:,1));
foldChange3 = (dataProb - meanProb)./meanProb;
plot(-pret:binsize:postt,smoothvect(foldChange3,g1),'LineWidth',2)
ylabel('Event probability (fold change)')
xlabel('Time from trigger onset')
if peakAlign == 0
    title([alignTo '-' 'trough'])
else
    title([alignTo '-' 'peak'])
end
set(gcf, 'renderer', 'painters')
yyaxis right

figure; hold on
pl1 = plot([-pret:binsize:postt],mean(CA1react),'-k','LineWidth',1);
boundedline([-pret:binsize:postt],mean(CA1react),std(CA1react)./sqrt(length(CA1react(:,1))),'-m');
pl2 = plot([-pret:binsize:postt],mean(PFCreact),'-r','LineWidth',1);
boundedline([-pret:binsize:postt],mean(PFCreact),std(PFCreact)./sqrt(length(PFCreact(:,1))),'-r');

pl3 = plot([-pret:binsize:postt],mean(CA1reactShuf),'-k','LineWidth',1);
boundedline([-pret:binsize:postt],mean(CA1reactShuf),std(CA1reactShuf)./sqrt(length(CA1reactShuf(:,1))),'--m');
pl4 = plot([-pret:binsize:postt],mean(PFCreactShuf),'-r','LineWidth',1);
boundedline([-pret:binsize:postt],mean(PFCreactShuf),std(PFCreactShuf)./sqrt(length(PFCreactShuf(:,1))),'--r');

[M pkBinPFC] = max(mean(PFCreact));
[M pkBinCA1] = max(mean(CA1react));
[M pkBinCA1rip] = max(smoothvect(foldChange1,g1));
[M pkBinPFCrip] = max(smoothvect(foldChange2,g1));
[M pkBinPFCspin] = max(smoothvect(foldChange3,g1));

%sig test vs shuffle
[p1 h] = ranksum(PFCreact(:,pkBinPFC),PFCreactShuf(:,pkBinPFC))
[p2 h] = ranksum(CA1react(:,pkBinCA1),CA1reactShuf(:,pkBinCA1))

%sig test for PFC vs CA1 timing
[M I] = max(CA1react(:,peakbins)');
[S SI] = sort(M,'descend');
figure
imagesc(-pret:binsize:postt,1:length(CA1react(:,1)),CA1react(SI,:))
colormap(inferno)
title('CA1 Assemblies')
caxis([-1 5])
set(gcf, 'renderer', 'painters')

[M I2] = max(PFCreact(:,peakbins)');
[S SI] = sort(M,'descend');
figure
imagesc(-pret:binsize:postt,1:length(PFCreact(:,1)),PFCreact(SI,:))
colormap(inferno)
title('PFC Assemblies')
caxis([-1 5])
set(gcf, 'renderer', 'painters')

[p3 h] = ranksum(I,I2); % peak timing
keyboard
