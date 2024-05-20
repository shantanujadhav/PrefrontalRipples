function jds_updown_replayProbabilityCoordNoncoord_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates and plots the probability of replay during up and down states
%of the slow oscillation for different SWR types
%%------------------------------------------------------------------------
binsize = 100;
day = 1;

pret=2000; postt=2000; %% Larger Window

rwin = [-2000 2000];
bwin = [-3000 -2000];
g1 = gaussian(3, 3);

dataCompileUp = [];
dataCompileDown = [];

alignTo = 'SO'; %SO or delta
peakAlign = 0; %align to peak (0) or trough (1)

upreplayprop = [];
downreplayprop = [];
epochs = 3:2:17;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    load(sprintf('%s%sslowoscdeltatimesSep_SWS%02d.mat', dir,animalprefix,day)) %load delta and so

    for ep = 1:length(epochs)

        epoch = epochs(ep);

        load(sprintf('%s%sPCCreplaydecode_alldata_50msRipNoIriCrit_sequenceScores0%d_%02d.mat',dir,animalprefix,day,epoch));
        pMats = PCC_events.pMats;
        upTimes = [slowosc{day}{epoch}.SOd2u slowosc{day}{epoch}.SOu2d]; %down to up - to - up to down transition
        downTimes = [slowosc{day}{epoch}.SOstarttime slowosc{day}{epoch}.SOd2u]; %SO start to down to up transition
        replayTimes = [];
        for event = 1:length(pMats)
            tmp_pmat = pMats{event};
            traj = find(min(tmp_pmat.traj_pvals)==tmp_pmat.traj_pvals);
            if length(traj) > 1
                continue
            end
            if min(tmp_pmat.traj_pvals) < 0.05
                if pMats{event}.rippletype == 1
                    eventtime = pMats{event}.eventtime;
                    replayTimes = [replayTimes; eventtime(1)];%(((eventtime(2) - eventtime(1))/2) + eventtime(1))]; %centers
                end
            end
        end
        if length(replayTimes) > 10 %at least 10 significant replay events
            upBins = periodAssign(replayTimes, upTimes);
            downBins = periodAssign(replayTimes, downTimes);
            upreplayprop = [upreplayprop; length(find(upBins~=0))/length(replayTimes)];
            downreplayprop = [downreplayprop; length(find(downBins~=0))/length(replayTimes)];
        end
    end
end

[p1 h1] = ranksum(upreplayprop,downreplayprop)
datacombinedReplay = [upreplayprop; downreplayprop];
g1 = repmat({'SO-Up replay'},length(upreplayprop),1);
g2 = repmat({'SO-Down replay'},length(downreplayprop),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedReplay,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['CA1 replay probability during SO-p = ' num2str(p1)])
set(gcf, 'renderer', 'painters')
keyboard
