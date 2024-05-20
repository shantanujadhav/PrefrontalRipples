function jds_periSO_eventProbability_ReactivationSuppression_M(animalprefixlist)

%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates the relationship between peri-SO CA1 reactivation and
%Independent PFC ripple triggered CA1 assembly suppression
%%------------------------------------------------------------------------

binsize = 20;
day = 1;
reactbins = 100;
reactbins2 = 400; %bins for PFC ripl aligned (+-8 sec win)

pret=2000; postt=2000; %% Larger Window
peakbins = find(-pret:binsize:postt>=0 & -pret:binsize:postt<=200);
peakbins2 = find(abs(-reactbins2:reactbins2)<=10);

CA1react = [];
CA1reactSupp = [];
CA1reactShuf = [];
totalCA1lagrips = 0;
totalSOcoordrips = 0;
alignTo = 'SO'; %SO or delta
peakAlign = 1; %align to peak (1) or trough (0)

totalTrigs = 0;
epochs = 3:2:17; %no epoch 1 because no reactivation

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    load(sprintf('%s%sctxrippletime_noncoordSWS%02d.mat', dir,animalprefix,day))
    load(sprintf('%s%srippletime_leadlag%02d.mat', dir,animalprefix,day)) %load delta and so
    load(sprintf('%s%sslowoscdeltatimesSep_SWS%02d.mat', dir,animalprefix,day)) %load delta and so
    reactHp = load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    reactCtx = load(sprintf('%s%sPFC_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    for ep = 1:length(epochs)
        epoch = epochs(ep);
        ctxrips = ctxripple{day}{epoch}.starttime;
        if length(ctxrips) < 10
            continue;
        end
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

        totalTrigs = totalTrigs + length(alignTimes);
        alignTimes = alignTimes*1000; %to ms

        ripCoord = ripplecoupling{day}{epoch}.hplagriptimes;
        totalCA1lagrips = totalCA1lagrips + length(ripCoord);

        if isempty(alignTimes)
            continue
        end
        allEvs = [];
        rTvec = assemblytmpHp{1}(:,1)*1000;
        react_idx = [];
        if ~isempty(ripCoord)
            ripCoord = ripCoord(:,1)*1000;
            for t = 1:length(alignTimes)
                tDiff = abs(alignTimes(t) - ripCoord);
                if ~isempty(find(ripCoord > (alignTimes(t)-250) &...
                        ripCoord < (alignTimes(t)+250)))
                    idxtmp = lookup(alignTimes(t), rTvec);
                    react_idx = [react_idx; idxtmp];
                end
            end
        end
        totalSOcoordrips = totalSOcoordrips + length(react_idx);
        if isempty(react_idx)
            continue
        end
        ctxrip_idx = [];
        for t = 1:length(ctxrips)
            idxtmp = lookup(ctxrips(t), assemblytmpHp{1}(:,1));
            ctxrip_idx = [ctxrip_idx; idxtmp];
        end
        
        % CA1 aligned to SO trough
        for c = 1:length(assemblytmpHp)
            CA1reactTmp = [];
            CA1reactTmpShuf = [];
            strengthstmp = assemblytmpHp{c}(:,2);
            strengthstmpShuf = circshift(strengthstmp, randi(length(strengthstmp)));
            for r = 1:length(react_idx)
                if ((react_idx(r) + reactbins) < length(strengthstmp)) &&...
                        ((react_idx(r) - reactbins) > 1)
                    %get vector of reactivation strenths for specified time period
                    tmp = strengthstmp((react_idx(r) - reactbins):(react_idx(r) + reactbins)); 
                    CA1reactTmp = [CA1reactTmp; tmp'];
                end
            end
            react_z = zscore(mean(CA1reactTmp,1));
            react_z = mean(react_z(peakbins));
            CA1react = [CA1react; react_z];
            % Shuffle each assembly once
            for s = 1
                strengthstmpShuf = circshift(strengthstmp, randi(length(strengthstmp)));
                for r = 1:length(react_idx)
                    if ((react_idx(r) + reactbins) < length(strengthstmpShuf)) &&...
                            ((react_idx(r) - reactbins) > 1)
                        tmp = strengthstmpShuf((react_idx(r) - reactbins):(react_idx(r) + reactbins));
                        CA1reactTmpShuf = [CA1reactTmpShuf; tmp'];
                    end
                end
                react_z = zscore(mean(CA1reactTmpShuf,1));
                react_z = mean(react_z(peakbins));
                CA1reactShuf = [CA1reactShuf; react_z];
            end
            % CA1 aligned to Ind. PFC ripples to get suppression
            atmp = [];
            for r = 1:length(ctxrip_idx)
                if ((ctxrip_idx(r) + reactbins2) < length(strengthstmp)) &&...
                        ((ctxrip_idx(r) - reactbins2) > 1)
                    tmp = strengthstmp((ctxrip_idx(r) - reactbins2):(ctxrip_idx(r) + reactbins2)); %get vector of reactivation strenths for specified time period
                    atmp = [atmp; tmp'];
                end
            end
            react_z = zscore(mean(atmp));
            react_z = mean(react_z(peakbins2));
            CA1reactSupp = [CA1reactSupp; react_z];
        end
    end
end
%% CA1
[r p] = corrcoef(CA1react,CA1reactSupp)
[shufr shufp] = corrcoef(CA1reactShuf,CA1reactSupp)
combinedData = [CA1react CA1reactSupp];
combinedShuf = [CA1reactShuf CA1reactSupp];

quartile_sep = floor(length(combinedData(:,1))/4);
suppressionquar = sortrows(combinedData,1);

vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = suppressionquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = suppressionquar(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});

X = [1:4];
prank = ranksum(vals{1}(:,2),vals{4}(:,2));

figure; hold on; errorbar(X, data_means(:,2), data_sems(:,2),'b','LineWidth',3);
xlim([0.5 4.5])

ylabel(['Suppression - pval Q1vQ4 - p=' num2str(prank)])

suppressionquarshuf = sortrows(combinedShuf,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = suppressionquarshuf(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = suppressionquarshuf(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});
errorbar(X, data_means(:,2), data_sems(:,2),'k','LineWidth',3);

ylabel('Suppression')
xlabel('SO aligned reactivation quartile')

legend({['Data-p = ' num2str(p(1,2)) 'r=' num2str(r(1,2))] ['Shuf-p = ' ...
    num2str(shufp(1,2)) 'r=' num2str(shufr(1,2))]})
title('CA1 SO reactivation PFCrip Suppression - Quartile')
xlabel([num2str(totalSOcoordrips) ' out of ' num2str(totalCA1lagrips) ' - ' ...
    num2str(totalSOcoordrips/totalCA1lagrips)])

xticks([1:4])
set(gcf, 'renderer', 'painters')

keyboard;