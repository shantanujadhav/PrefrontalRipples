function jds_PCCrippleReplayDifference_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates and plots the PCC difference for mod and non mod cells.
%Compares coordinated and independent SWR replay events.
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
        events_sig_ind = [];
        for event = 1:length(pMats)
            tmp_pmat = pMats{event};
            if tmp_pmat.rippletype == 1
                traj = find(min(tmp_pmat.traj_pvals)==tmp_pmat.traj_pvals);
                if min(tmp_pmat.traj_pvals) < 0.05
                    idx1 = find(events(:,7) == event);
                    idx2coord = find(events(:,2) == traj);
                    idx3coord = intersect(idx1, idx2coord);
                    events_sig_ind = [events_sig_ind; events(idx3coord,:)];
                end
            end
        end

        events_sig_coord = [];
        for event = 1:length(pMats)
            tmp_pmat = pMats{event};
            if tmp_pmat.rippletype == 2
                traj = find(min(tmp_pmat.traj_pvals)==tmp_pmat.traj_pvals);
                if min(tmp_pmat.traj_pvals) < 0.05
                    idx1 = find(events(:,7) == event);
                    idx2coord = find(events(:,2) == traj);
                    idx3coord = intersect(idx1, idx2coord);
                    events_sig_coord = [events_sig_coord; events(idx3coord,:)];
                end
            end
        end

        if (~isempty(events_sig_ind)) && ~isempty(events_sig_coord)
            for i = 1:length(cellidx(:,1))
                currCell = cellidx(i,[1 2]);
                [b cellLoc] = ismember(currCell, modcells(:,[1,2]), 'rows', 'legacy');
                if (cellLoc ~= 0)
                    cellmod = modcells(cellLoc,3);
                    
                    idx1ind = find(events_sig_ind(:,3) == currCell(1));
                    idx2ind = find(events_sig_ind(:,4) == currCell(2));
                    idx3ind = intersect(idx1ind, idx2ind);
                    PCCind = nanmean(events_sig_ind(idx3ind,1));

                    idx1coord = find(events_sig_coord(:,3) == currCell(1));
                    idx2coord = find(events_sig_coord(:,4) == currCell(2));
                    idx3coord = intersect(idx1coord, idx2coord);
                    PCCcoord = nanmean(events_sig_coord(idx3coord,1));
                    
                    if (length(idx3coord) > 5) && (length(idx3ind) > 5) ...
                            && (~isinf(PCCind)) && (~isinf(PCCcoord))
                        if cellmod == 1
                            PCCexc = [PCCexc; PCCcoord-PCCind];
                        elseif cellmod == -1
                            PCCinh = [PCCinh; PCCcoord-PCCind];
                        elseif cellmod == 0
                            PCCnonmod = [PCCnonmod; PCCcoord-PCCind];
                        end
                    end
                    
                elseif cellLoc == 0
                    cellmod = 0;
                    idx1ind = find(events_sig_ind(:,3) == currCell(1));
                    idx2ind = find(events_sig_ind(:,4) == currCell(2));
                    idx3ind = intersect(idx1ind, idx2ind);
                    PCCind = nanmean(events_sig_ind(idx3ind,1));

                    idx1coord = find(events_sig_coord(:,3) == currCell(1));
                    idx2coord = find(events_sig_coord(:,4) == currCell(2));
                    idx3coord = intersect(idx1coord, idx2coord);
                    PCCcoord = nanmean(events_sig_coord(idx3coord,1));
                    if (length(idx3coord) > 5) && (length(idx3ind) > 5) ...
                            && (~isinf(PCCind)) && (~isinf(PCCcoord))
                        PCCnonmod = [PCCnonmod; PCCcoord-PCCind];
                    end
                end
            end
        end
    end
end

PCCmod = [PCCexc;PCCinh];
datameans = [nanmean(PCCmod) nanmean(PCCnonmod)];
datasems = [(nanstd(PCCmod)/sqrt(length(find(~isnan(PCCmod)))))...
    (nanstd(PCCnonmod)/sqrt(length(find(~isnan(PCCnonmod)))))];
[p1 h1] = ranksum(PCCmod,PCCnonmod)
figure
bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('PCC')
title(['PCC difference (Coordinated - Independent) p=' num2str(p1)])
xticklabels({'CA1mod','CA1nonmod'}); xtickangle(45)

keyboard

