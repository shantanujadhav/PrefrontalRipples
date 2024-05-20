function jds_sequenceRobustness_withawake_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates sequence degradation for significant replay events in wake and
%sleep (independent and coordinated in sleep). This is a measure of
%sequence robustness or fidelity.
%%------------------------------------------------------------------------
day = 1;

rZDegradedInd = [];
rZDegradedCoord = [];
rZDegradedWake = [];

all_rZ_real_sleep = [];
all_rZ_real_wake = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    epochs = [3:2:17];
    runepochs = [2:2:16];

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    for ep=1:length(epochs)

        epoch = epochs(ep);

        load(sprintf('%s%sPCCreplaydecode_alldata_50msRipNoIriCrit_sequenceScores0%d_%02d.mat',...
            dir,animalprefix,day,epoch));
        sleepPcc = PCC_events; clear PCC_events
        
        sleepEvents = sleepPcc.alleventsAllcell;

        sleepCellidx = sleepPcc.cellidx;
        sleep_pMats = sleepPcc.pMats;

        %Filter for low pvalue
        events_sig_sleep = [];
        for event = 1:length(sleep_pMats)
            tmp_pmat = sleep_pMats{event};
            traj = find(min(tmp_pmat.traj_pvals)==tmp_pmat.traj_pvals);
            if min(tmp_pmat.traj_pvals) < 0.05
                idx1 = find(sleepEvents(:,7) == event);
                idx2 = find(sleepEvents(:,2) == traj);
                idx3 = intersect(idx1, idx2);
                events_sig_sleep = [events_sig_sleep; sleepEvents(idx3,:)];
            end
        end

        uniqueEvents = unique(events_sig_sleep(:,7));
        for ev = 1:length(uniqueEvents)
            idxtmp = find(events_sig_sleep(:,7) == uniqueEvents(ev));
            tmp = events_sig_sleep(idxtmp(1),5);
            all_rZ_real_sleep = [all_rZ_real_sleep; tmp];
        end

        if ~isempty(events_sig_sleep)
            for i = 1:length(sleepCellidx(:,1))
                currCell = sleepCellidx(i,[1 2]);

                idx1 = find(events_sig_sleep(:,3) == currCell(1));
                idx2 = find(events_sig_sleep(:,4) == currCell(2));
                idx3 = intersect(idx1, idx2);
                seqScore = events_sig_sleep(idx3,6);
                seqScoreReal = events_sig_sleep(idx3,5);

                eventidx = events_sig_sleep(idx3,7);
                notnanidx = find(~isnan(seqScore));
                seqScore = seqScore(notnanidx);
                eventidx = eventidx(notnanidx);

                for r = 1:length(eventidx)
                    rtype = sleepPcc.pMats{eventidx(r)}.rippletype;
                    if rtype == 1
                        rZDegradedInd = [rZDegradedInd; seqScore(r)];
                    elseif rtype == 2
                        rZDegradedCoord = [rZDegradedCoord; seqScore(r)];
                    end
                end
            end
        end
    end
    %% WAKE
    for runep = 1:length(runepochs)
        repoch = runepochs(runep);
        load(sprintf('%s%sPCCreplaydecode_run_alldataNew_50msRipNoIriCritTbin15_sequenceScores0%d_%02d.mat',...
            dir,animalprefix,day,repoch));
        wakePcc = PCC_events; clear PCC_events

        wakeEvents = wakePcc.alleventsAllcell;

        wakeCellidx = wakePcc.cellidx;
        wake_pMats = wakePcc.pMats;
        events_sig_wake = [];
        for event = 1:length(wake_pMats)
            tmp_pmat = wake_pMats{event};
            traj = find(min(tmp_pmat.traj_pvals)==tmp_pmat.traj_pvals);
            if min(tmp_pmat.traj_pvals) < 0.05
                idx1 = find(wakeEvents(:,7) == event);
                idx2 = find(wakeEvents(:,2) == traj);
                idx3 = intersect(idx1, idx2);
                events_sig_wake = [events_sig_wake; wakeEvents(idx3,:)];
            end
        end

        uniqueEvents = unique(events_sig_wake(:,7));
        for ev = 1:length(uniqueEvents)
            idxtmp = find(events_sig_wake(:,7) == uniqueEvents(ev));
            tmp = events_sig_wake(idxtmp(1),5);
            all_rZ_real_wake = [all_rZ_real_wake; tmp];
        end

        if ~isempty(events_sig_wake)
            for i = 1:length(wakeCellidx(:,1))
                currCell = wakeCellidx(i,[1 2]);

                idx1 = find(events_sig_wake(:,3) == currCell(1));
                idx2 = find(events_sig_wake(:,4) == currCell(2));
                idx3 = intersect(idx1, idx2);
                seqScore = events_sig_wake(idx3,6);
                seqScoreReal = events_sig_wake(idx3,5);

                eventidx = events_sig_wake(idx3,7);
                notnanidx = find(~isnan(seqScore));
                seqScore = seqScore(notnanidx);
                eventidx = eventidx(notnanidx);

                for r = 1:length(eventidx)
                    rZDegradedWake = [rZDegradedWake; seqScore(r)];
                end
            end
        end
    end
end
%%
n1 = [];
n1(1:length(rZDegradedInd)) = 1;
n1 = string(n1)';

n2 = [];
n2(1:length(rZDegradedCoord)) = 2;
n2 = string(n2)';

n3 = [];
n3(1:length(rZDegradedWake)) = 3;
n3 = string(n3)';

group_data = [rZDegradedInd; rZDegradedCoord; rZDegradedWake];
grouping = [n1;n2;n3];

[p,tbl,stats] = kruskalwallis(group_data,grouping);
[c,m,h,gnames] = multcompare(stats,"CriticalValueType","bonferroni");
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

datameans = [nanmean(rZDegradedInd) nanmean(rZDegradedCoord) nanmean(rZDegradedWake)];
datasems = [nanstd(rZDegradedInd)./sqrt(length(rZDegradedInd))...
    nanstd(rZDegradedCoord)./sqrt(length(rZDegradedCoord))...
    nanstd(rZDegradedWake)./sqrt(length(rZDegradedWake))];

figure
bar([1:3],datameans,'k')
hold on
er = errorbar([1:3],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Sequence Score (rZ)')
xticklabels({'Independent','Coordinated','Wake'}); xtickangle(45)
title('Sequence degradation')

keyboard




