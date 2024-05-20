function jds_periSO_eventProbability_ReactivationQuartiles_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plots peri-SO CA1 reactivation strength based on quartiles. Plots 1st and
%4th quartiles, and heatmap of sorted assemblies. Accompanies figure
%showing relationship between peri-SO CA1 reactivation strength and
%independent PFC ripple triggered suppression.
%%------------------------------------------------------------------------
binsize = 20;
day = 1;
reactbins = 100;


pret=2000; postt=2000; 
peakbins = find(-pret:binsize:postt<=250 & -pret:binsize:postt>=0);

g1 = gaussian(3, 3);

CA1react = [];
CA1reactStrength = [];
PFCreact = [];
CA1reactShuf = [];
PFCreactShuf = [];

alignTo = 'SO'; %SO or delta
peakAlign = 0; %align to peak (0) or trough (1)

% evList = {'rippletime_coordSWS','ctxrippletime_coordSWS','ctxspindletime_SWS'}; %coordinated ripples
% evList = {'rippletime_noncoordSWS','ctxrippletime_noncoordSWS','ctxspindletime_SWS'}; %noncoordinated ripples
evList = {'rippletime_leadlag','rippletime_leadlag','ctxspindletime_SWS'};
CA1lead = 0;
totalTrigs = 0;
epochs = 3:2:17;
Q1h = [];
Q2h = [];
Q3h = [];
Q4h = [];
Q1c = [];
Q2c = [];
Q3c = [];
Q4c = [];
Q1xcorr = [];
Q2xcorr = [];
Q3xcorr = [];
Q4xcorr = [];

Q1xcorrOther = [];
Q2xcorrOther = [];
Q3xcorrOther = [];
Q4xcorrOther = [];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    load(sprintf('%s%srippletime_leadlag%02d.mat', dir,animalprefix,day)) %load delta and so
    load(sprintf('%s%sslowoscdeltatimesSep_SWS%02d.mat', dir,animalprefix,day)) %load delta and so
    reactHp = load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    reactCtx = load(sprintf('%s%sPFC_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    strengthVecHp = [];
    strengthVecCtx = [];
    cntCA1 = 1;
    cntPFC = 1;
    hpAssemNum = [];
    ctxAssemNum = [];
    hpRtimes = [];
    ctxRtimes = [];
    epVecHp = [];
    epVecCtx = [];
    alignTimesEp = [];
    for evs = 1:length(evList)
        for ep = 1:length(epochs)
            epoch = epochs(ep);
            assemblytmpHp = reactHp.RtimeStrength{epoch}.reactivationStrength;
            assemblytmpCtx = reactCtx.RtimeStrength{epoch}.reactivationStrength;
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
            end

            totalTrigs = totalTrigs + length(alignTimes);
            evTmp = [];
            evTmpBck = [];
            alignTimes = alignTimes*1000; %to ms
            if length(alignTimes) < 10
                continue
            end

            if isempty(currEvents)
                continue
            end

            ripCoord = ripplecoupling{day}{epoch}.hplagriptimes;
            currEvents = currEvents*1000;
            allEvs = [];
            rTvec = assemblytmpHp{1}(:,1)*1000;
            %find ripples in vicinity of SO
            if ~isempty(ripCoord)
                ripCoord = ripCoord(:,1)*1000;
                react_idx = [];
                alignIdx = [];
                for t = 1:length(alignTimes)
                    tDiff = abs(alignTimes(t) - ripCoord);
                    if ~isempty(find(ripCoord > (alignTimes(t)-250) & ripCoord < (alignTimes(t)+250)))
                        idxtmp = lookup(alignTimes(t), rTvec);
                        react_idx = [react_idx; idxtmp];
                        alignIdx = [alignIdx; t];
                    end
                end
                alignTimesEp{epoch} = [alignTimes(alignIdx)-1000 alignTimes(alignIdx)+1000];
            end

            if evs == 1 %only align reactivation strength once
                %ALIGN BOTH CA1 and PFC
                if length(alignTimes(alignIdx)) > 10
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
                        rIdx = find(assemblytmpHp{c}(:,2) > 5);
                        hpRtimes{cntCA1} = assemblytmpHp{c}(rIdx,1);
                        aTmp = smoothvect(zscore(mean(CA1reactTmp,1)),g1);
                        CA1react = [CA1react; smoothvect(zscore(mean(CA1reactTmp,1)),g2)];
                        CA1reactStrength = [CA1reactStrength; mean(aTmp(peakbins))];
                        hpAssemNum(cntCA1) = cntCA1;
                        cntCA1 = cntCA1 + 1;
                        epVecHp = [epVecHp; epoch];
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
                        rIdx = find(assemblytmpCtx{d}(:,2) > 5);
                        ctxRtimes{cntPFC} = assemblytmpCtx{d}(rIdx,1);
                        aTmp = smoothvect(zscore(mean(PFCreactTmp,1)),g1);
                        PFCreact = [PFCreact; smoothvect(zscore(mean(PFCreactTmp,1)),g1)];
                        strengthVecCtx = [strengthVecCtx; mean(aTmp(peakbins))];
                        ctxAssemNum(cntPFC) = cntPFC;
                        cntPFC = cntPFC + 1;
                        epVecCtx = [epVecCtx; epoch];
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
            end
        end
    end

    quartile_sepHp = floor(length(strengthVecHp(:,1))/4);
    [HpQuar Idx] = sortrows(strengthVecHp,1);
    epVecHp = epVecHp(Idx);
    valsHp = [];
    cnt = 1;
    for s = 1:4
        if s < 4
            tmp = HpQuar(cnt:quartile_sepHp*s);
            tmp(:,2) = s;
            tmp(:,3) = Idx(cnt:quartile_sepHp*s);
            tmp(:,4) = epVecHp(cnt:quartile_sepHp*s);
            if s == 1
                Q1h = [Q1h; CA1react(Idx(cnt:quartile_sepHp*s),:)];
            elseif s == 2
                Q2h = [Q2h; CA1react(Idx(cnt:quartile_sepHp*s),:)];
            elseif s == 3
                Q3h = [Q3h; CA1react(Idx(cnt:quartile_sepHp*s),:)];
            end
        else
            tmp = HpQuar(cnt:end);
            tmp(:,2) = s;
            tmp(:,3) = Idx(cnt:end);
            tmp(:,4) = epVecHp(cnt:end);
            Q4h = [Q4h; CA1react(Idx(cnt:end),:)];
        end
        valsHp = [valsHp; tmp];
        cnt = cnt + quartile_sepHp;
        clear tmp
    end

    quartile_sepCtx = floor(length(strengthVecCtx(:,1))/4);
    [CtxQuar Idx] = sortrows(strengthVecCtx,1);
    epVecCtx = epVecCtx(Idx);
    valsCtx = [];
    cnt = 1;
    for s = 1:4
        if s < 4
            tmp = CtxQuar(cnt:quartile_sepCtx*s);
            tmp(:,2) = s;
            tmp(:,3) = Idx(cnt:quartile_sepCtx*s);
            tmp(:,4) = epVecCtx(cnt:quartile_sepCtx*s);
            if s == 1
                Q1c = [Q1c; PFCreact(Idx(cnt:quartile_sepCtx*s),:)];
            elseif s == 2
                Q2c = [Q2c; PFCreact(Idx(cnt:quartile_sepCtx*s),:)];
            elseif s == 3
                Q3c = [Q3c; PFCreact(Idx(cnt:quartile_sepCtx*s),:)];
            end
        else
            tmp = CtxQuar(cnt:end);
            tmp(:,2) = s;
            tmp(:,3) = Idx(cnt:end);
            tmp(:,4) = epVecCtx(cnt:end);
            Q4c = [Q4c; PFCreact(Idx(cnt:end),:)];
        end
        valsCtx = [valsCtx; tmp];
        cnt = cnt + quartile_sepCtx;
        clear tmp
    end
end

figure; hold on
pl1 = plot([-1000:10:1000-1],mean(Q1h),'-m','LineWidth',1);
boundedline([-1000:10:1000-1],mean(Q1h),std(Q1h)./sqrt(length(Q1h(:,1))),'-m');
pl2 = plot([-1000:10:1000-1],mean(Q4h),'-m','LineWidth',1);
boundedline([-1000:10:1000-1],mean(Q4h),std(Q4h)./sqrt(length(Q4h(:,1))),'-m');
xlim([-1000 1000])

[B I] = sort(CA1reactStrength);

figure; hold on
%sorted by reactivation strength in peak bins
imagesc(CA1react(I,:))

keyboard
