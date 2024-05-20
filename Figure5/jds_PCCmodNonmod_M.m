function jds_PCCmodNonmod_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Compare per cell comtribution (PCC) to replay for mod and nonmod CA1 cells
%%------------------------------------------------------------------------
day = 1;

PCCexc = [];
PCCinh = [];
PCCnonmod = [];
sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    
    if isequal(animalprefix, 'JS34')
        epochs = [3:2:17];
    elseif (isequal(animalprefix, 'JS17')) || (isequal(animalprefix, 'KL8'))
        epochs = [3:2:17];
    else
        epochs = [1:2:17];
    end
    
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    
    load(sprintf('%s%sCA1ctxripallmodNewWin_epsIncludeHigh0%d.mat',dir,animalprefix,day));
    modcellsall = epochModulation.cellidx;
    modMat = epochModulation.modMat;
    
    for ep=1:length(epochs)
        
        epoch = epochs(ep);

        ep2 = find(sleeps(:,2) == epoch);
        
        load(sprintf('%s%sPCCreplaydecode_alldata_50msRipNoIriCrit_sequenceScores0%d_%02d.mat',dir,animalprefix,day,epoch));
        
        events = PCC_events.alleventsAllcell;
        cellidx = modcellsall;
        modcells = [cellidx modMat(:,ep2)];
        pMats = PCC_events.pMats;

        %Filter for low pvalue
        events_sig = [];
        for event = 1:length(pMats)
            tmp_pmat = pMats{event};
%             if tmp_pmat.rippletype == 1
                traj = find(min(tmp_pmat.traj_pvals)==tmp_pmat.traj_pvals);
                if min(tmp_pmat.traj_pvals) < 0.05
                    idx1 = find(events(:,7) == event);
                    idx2 = find(events(:,2) == traj);
                    idx3 = intersect(idx1, idx2);
                    events_sig = [events_sig; events(idx3,:)];
                end
%             end
        end
%         events_sig = events;
        
        if ~isempty(events_sig)
            for i = 1:length(cellidx(:,1))
                currCell = cellidx(i,[1 2]);
                [b cellLoc] = ismember(currCell, modcells(:,[1,2]), 'rows', 'legacy');
                if (cellLoc ~= 0)
                    cellmod = modcells(cellLoc,3);
                    
                    idx1 = find(events_sig(:,3) == currCell(1));
                    idx2 = find(events_sig(:,4) == currCell(2));
                    idx3 = intersect(idx1, idx2);
                    PCC = nanmean(events_sig(idx3,1));
                    
                    if (length(idx3) > 10) && (~isinf(PCC))
                        if cellmod == 1
                            PCCexc = [PCCexc; PCC];
                        elseif cellmod == -1
                            PCCinh = [PCCinh; PCC];
                        elseif cellmod == 0
                            PCCnonmod = [PCCnonmod; PCC];
                        end
                    end
                    
                elseif cellLoc == 0
                    cellmod = 0;
                    idx1 = find(events_sig(:,3) == currCell(1));
                    idx2 = find(events_sig(:,4) == currCell(2));
                    idx3 = intersect(idx1, idx2);
                    PCC = nanmean(events_sig(idx3,1));

                    if (length(idx3) > 10) && (~isinf(PCC))
                        PCCnonmod = [PCCnonmod; PCC];
                    end
                end
            end
        end
    end
end

PCCmod = [PCCinh;PCCexc];
datameans = [nanmean(PCCexc) nanmean(PCCinh)];
datasems = [(nanstd(PCCexc)/sqrt(length(find(~isnan(PCCexc)))))...
    (nanstd(PCCinh)/sqrt(length(find(~isnan(PCCinh)))))];
datamedians = [nanmedian(PCCexc) nanmedian(PCCinh)];

figure
bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('PCC')
title('Per Cell Contribution (PCC) to Replay Sequenceness')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

datameans = [nanmean(PCCexc) nanmean(PCCinh)];
datasems = [(nanstd(PCCexc)/sqrt(length(find(~isnan(PCCexc)))))...
    (nanstd(PCCinh)/sqrt(length(find(~isnan(PCCinh)))))];


[p1 h1] = ranksum(PCCmod,PCCnonmod)
datameans = [nanmean(PCCmod) nanmean(PCCnonmod)];
datasems = [(nanstd(PCCmod)/sqrt(length(find(~isnan(PCCmod)))))...
    (nanstd(PCCnonmod)/sqrt(length(find(~isnan(PCCnonmod)))))];
datamedians = [nanmedian(PCCmod) nanmedian(PCCnonmod)];

figure
bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('PCC')
title(['Per Cell Contribution (PCC) to Replay Sequenceness p=' num2str(p1)])
xticklabels({'CA1mod','CA1nonmod'}); xtickangle(45)
keyboard