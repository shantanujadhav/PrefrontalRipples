function jds_CA1_SWRdimensionality_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates the dimensionality of independent and coordinated SWR events
%and compares
%%------------------------------------------------------------------------

day = 1; 
epochs = [1:2:17];
dim_coord = [];
dim_noncoord = [];
for a = 1:length(animalprefixlist)

    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
    % get ripple time
    load(sprintf('%s%srippletime_coordSWS0%d.mat',dir,animalprefix,day));
    cripple = ripple; clear ripple
    load(sprintf('%s%srippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    ncripple = ripple; clear ripple
    for e = 1:length(epochs)
        ep = epochs(e);

        crip = cripple{day}{ep};
        ncrip = ncripple{day}{ep};

        criptimes = [crip.starttime crip.endtime];
        numcoordrips = length(criptimes(:,1));
        ncriptimes = [ncrip.starttime ncrip.endtime];
        numnoncoordrips = length(ncriptimes(:,1));

        %at least 100 ripples of each type in an epoch to be considered
        if (length(criptimes(:,1)) > 100) && (length(ncriptimes(:,1)) > 100)
            for samps = 1:100
                randsamp_nc = randi(numnoncoordrips,[1 50]);
                randsamp_c = randi(numcoordrips,[1 50]);
                ncriptimes_samp = ncriptimes(randsamp_nc,:);
                criptimes_samp = criptimes(randsamp_c,:);
                if samps == 1
                    %include interneurons and pyramidal cells
                    [ctxidx, hpidx] = jds_getallepcells_includeall(dir, animalprefix, day, ep, []);
                end
                hpnum = length(hpidx(:,1));
                if hpnum < 20
                    continue
                end

                celldata = [];
                CA1matrix_nc = [];
                CA1matrix_c = [];
                for cellcount = 1:hpnum %get spikes for each cell
                    index = [day,ep,hpidx(cellcount,:)] ;
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                    else
                        spiketimes = [];
                    end
                    spikebins = periodAssign(spiketimes, ncriptimes_samp(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
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
                    spikecount = zeros(1,size(ncriptimes_samp,1));
                    for i = 1:length(spikebins)
                        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                    end
                    CA1matrix_nc = [CA1matrix_nc spikecount']; %concatenating num spikes per cell, per event
                end

                for cellcount = 1:hpnum %get spikes for each cell
                    index = [day,ep,hpidx(cellcount,:)] ;
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                    else
                        spiketimes = [];
                    end
                    spikebins = periodAssign(spiketimes, criptimes_samp(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
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
                    spikecount = zeros(1,size(criptimes_samp,1));
                    for i = 1:length(spikebins)
                        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                    end
                    CA1matrix_c = [CA1matrix_c spikecount']; %concatenating num spikes per cell, per event
                end
                [COEFFnc, SCOREnc, LATENTnc, TSQUAREDnc, EXPLAINEDnc] = pca(CA1matrix_nc);
                eqn_nc = cumsum(EXPLAINEDnc) > 90;
                sol_nc = min(find(eqn_nc==1));
                dim_nc = sol_nc/hpnum;
                [COEFFc, SCOREc, LATENTc, TSQUAREDc, EXPLAINEDc] = pca(CA1matrix_c);
                eqn_c = cumsum(EXPLAINEDc) > 90;
                sol_c = min(find(eqn_c==1));
                dim_c = sol_c/hpnum;
                dim_noncoord = [dim_noncoord; dim_nc];
                dim_coord = [dim_coord; dim_c];
            end
        end
    end
end
[p h] = ranksum(dim_noncoord,dim_coord)
figure; hold on
bar([mean(dim_noncoord) mean(dim_coord)],'k')
errorbar([1:2],[mean(dim_noncoord) mean(dim_coord)],...
    [(std(dim_noncoord)./sqrt(length(dim_noncoord))) ...
    (std(dim_coord)./sqrt(length(dim_coord)))],'k','LineStyle','none')
ylim([0.2 0.245])
yticks([0.2:0.01:0.24])
ylabel('Scaled dimensionality')
xticks([1 2])
xticklabels({'Ind','Coord'})
title(['SWR dimensionality p = ' num2str(p)])
set(gcf, 'renderer', 'painters')

keyboard;