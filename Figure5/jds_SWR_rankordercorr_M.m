function jds_SWR_rankordercorr_M(animalprefixlist,day,eps,cellcountthresh)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates and plots the fold change occurence of ripples and spindles
%relative to slow oscillation troughs (up states)
%%------------------------------------------------------------------------

independentRho = [];
coordinatedRho = [];
independentRhoEp = [];
coordinatedRhoEp = [];
plotranks = 0;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    %%
    %-----match neurons across epochs-----%

    %%
    %-----create the event matrix during SWRs-----%
    spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
    load(sprintf('%s/%srippletime_noncoordSWS%02d.mat', dir, animalprefix, day));
    nc_rip = ripple; clear ripple
    load(sprintf('%s/%srippletime_coordSWS%02d.mat', dir, animalprefix, day));
    load(sprintf('%s/%sspikes%02d.mat', dir, animalprefix, day));
    c_rip = ripple; clear ripple
    for ep = eps
        [ctxidx, hpidx] = jds_getallepcells_includeall(dir, animalprefix, day, ep, []);
        hpnum = length(hpidx(:,1));

        nc_riptimes = [nc_rip{day}{ep}.starttime nc_rip{day}{ep}.endtime];
        nc_riptimes(:,3) = 1;

        c_riptimes = [c_rip{day}{ep}.starttime c_rip{day}{ep}.endtime];
        c_riptimes(:,3) = 2;

        if (length(nc_riptimes) > 10) && (length(c_riptimes) > 10)

            %combine riptimes
            riptimes = sortrows([nc_riptimes; c_riptimes],1);

            rip_starttime = 1000*riptimes(:,1);  % in ms
            rip_endtime = 1000*riptimes(:,2);  % in ms

            riplength = rip_endtime - rip_starttime;
            keepidx2 = find(riplength >= 50);% use the ripple last for more than 50 ms
            riptimes = riptimes(keepidx2,:);

            ripnum = size(riptimes,1);
            tmp_c = [];
            tmp_nc = [];
            %%
            if ripnum > 1
                celldata = [];
                spikecounts = [];
                for cellcount = 1:hpnum %get spikes for each cell
                    index = [day,ep,hpidx(cellcount,:)] ;
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                    else
                        spiketimes = [];
                    end
                    spikebins = periodAssign(spiketimes, riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
                    if ~isempty(spiketimes)
                        validspikes = find(spikebins);
                        spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
                        spikebins = spikebins(validspikes);
                        tmpcelldata = [spiketimes spikebins];
                    end
                    if ~isempty(spiketimes)
                        tmpcelldata(:,3) = cellcount; %keep count of how many cells active during rip event
                    else
                        tmpcelldata = [0 0 cellcount];
                    end
                    celldata = [celldata; tmpcelldata];
                    spikecount = zeros(1,size(riptimes,1));
                    for i = 1:length(spikebins)
                        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                    end
                    spikecounts = [spikecounts; spikecount]; %concatenating num spikes per cell, per event
                end
                cellcounts = sum((spikecounts > 0));
                eventindex = find(cellcounts >= cellcountthresh); %Find event indices that have more than cellthresh # of cells
                if ~isempty(eventindex)
                    for event = 1:length(eventindex)
                        rippletype = riptimes(eventindex(event),3);
                        cellsi = celldata(find(celldata(:,2)==eventindex(event)),3);
                        spktimes = celldata(find(celldata(:,2)==eventindex(event)),1);

                        [cellsi,ia] = unique(cellsi,'first');
                        spktimes = spktimes(ia);
                        [~, I] = sort(spktimes);
                        %                     [~,sortorder] = sort(ia);
                        I(:,2) = 1:length(I(:,1));
                        rankorder = [cellsi(I(:,1)) normalize(I(:,2),'range')];
                        templateranks = [];
                        for evs = 1:length(eventindex)
                            nanVec = nan(hpnum,1);
                            rippletype2 = riptimes(eventindex(evs),3);

                            if (evs ~= event) && (rippletype == rippletype2)
                                cellsi2 = celldata(find(celldata(:,2)==eventindex(evs)),3);
                                spktimes2 = celldata(find(celldata(:,2)==eventindex(evs)),1);
                                [cellsi2,ia2] = unique(cellsi2,'first');
                                spktimes2 = spktimes2(ia2);
                                [~, I2] = sort(spktimes2);
                                %                             [~,sortorder2] = sort(ia2);
                                I2(:,2) = 1:length(I2(:,1));
                                rankorder2 = [cellsi2(I2(:,1)) normalize(I2(:,2),'range')];
                                nanVec(rankorder2(:,1)) = rankorder2(:,2);
                                tmpRank = nanVec;
                                templateranks = [templateranks tmpRank];
                            end
                        end
                        templatemeanranks = nanmean(templateranks,2);
                        ranksems = (nanstd(templateranks')./sqrt(length(templateranks(1,:))))';
                        currcells = templatemeanranks(rankorder(:,1));

                        if plotranks == 1
                            [~, sortidx] = sort(templatemeanranks);
                            sortmean = templatemeanranks(sortidx);
                            sortsem = ranksems(sortidx);
                            errorbar([1:length(sortmean)],...
                                sortmean',sortsem','r','LineStyle','none','Marker','o','MarkerFaceColor','auto');
                            title('rank order template')
                            ylabel('Normalized rank')
                            xlabel('Cell #')
                            set(gcf, 'renderer', 'painters')
                            keyboard;
                        end

%                         rho = corr(currcells,rankorder(:,2),'Type','Spearman','rows','complete');
                        rho = corr(currcells,rankorder(:,2),'rows','complete');
                        

                        if rippletype == 1
                            independentRho = [independentRho; rho];
                            tmp_nc = [tmp_nc; rho];
                        elseif rippletype == 2
                            coordinatedRho = [coordinatedRho; rho];
                            tmp_c = [tmp_c; rho];
                        end
                    end
                end
            end
            independentRhoEp = [independentRhoEp; mean(tmp_nc)];
            coordinatedRhoEp = [coordinatedRhoEp; mean(tmp_c)];
        end
    end
end
keyboard
[p h] = ranksum(independentRho,coordinatedRho);
[p2 h2] = signrank(independentRhoEp,coordinatedRhoEp);
figure
bar([mean(independentRho) mean(coordinatedRho)],'k')
hold on
errorbar([1:2],[mean(independentRho) mean(coordinatedRho)],...
    [(std(independentRho)./sqrt(length(independentRho))) ...
    (std(coordinatedRho)./sqrt(length(coordinatedRho)))],'k','LineStyle','none')
ylim([0.28 0.315])
yticks([0.28:0.01:0.31])
xticklabels({'Ind','Coord'})
ylabel('Rank order correlation (r)')
title(['Rank order correlation p=' num2str(p)])
set(gcf, 'renderer', 'painters')

keyboard