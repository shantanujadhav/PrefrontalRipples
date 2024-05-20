function jds_rippletrig_spindlepwr_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates ripple triggered spindle power for different ripples
%%------------------------------------------------------------------------
day = 1;
epochs = [1:2:17];
allenv_coord = [];
allenv_noncoord = [];
for a = 1:length(animalprefixlist)

    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    daystring = ['0',num2str(day)];

    plotindiv = 0;

    load(sprintf('%s%sctxrippletime_noncoordSWS0%d.mat',dir,animalprefix,day));% get ripple time
    noncoord_ctx = ctxripple; clear ctxripple
    load(sprintf('%s%sctxrippletime_coordSWS0%d.mat',dir,animalprefix,day));
    coord_ctx = ctxripple;
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));

    tets = tetinfo{1}{1};

    pfctets = [];
    for t = 1:length(tets)
        tmp = tets{t};
        if isfield(tmp, 'descrip')
            if isequal(tmp.descrip, 'ctxriptet') %use ctx ripple tetrodes
                pfctets = [pfctets; t];
            end
        end
    end
    spintets = pfctets;

    for ep=1:length(epochs)
        epoch = epochs(ep);

        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end


        coord = coord_ctx{day}{epoch};
        noncoord = noncoord_ctx{day}{epoch};

        curreegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring,'-','01']; %use tet1 to extract starttime of epoch
        load(curreegfile);

        ep_starttime = eegref{1, 1}{epoch}{1}.starttime;

        times = geteegtimes(eegref{day}{epoch}{1});

        coordriplist = [coord.starttime coord.endtime];%.starttime relative to starttime of epoch

        noncoordriplist = [noncoord.starttime noncoord.endtime];%.starttime

        if ~isempty(coordriplist)
            allenv = [];
            for ii = 1:length(spintets)
                allenvtmp = [];
                tetrode = spintets(ii);
                if tetrode <10
                    tetstring = ['0',num2str(tetrode)];
                else
                    tetstring = num2str(tetrode);
                end
                currspinfile = [dir,'/EEG/',animalprefix,'spindlegnd', daystring,'-',epochstring,'-',tetstring];
                load(currspinfile);
                env = zscore(double(spindlegnd{day}{epoch}{tetrode}.data(:,3)));
                for i = 1:length(coordriplist(:,1))
                    riptimetmp = [coordriplist(i,1) coordriplist(i,2)];
                    tmp_idx = find(times >= riptimetmp(1) & times <= riptimetmp(2));
                    envtmp = max(env(tmp_idx));
                    allenvtmp = [allenvtmp envtmp];
                end
                allenv = [allenv; allenvtmp];
            end
            allenv_coord = [allenv_coord; mean(allenv)']; %mean over all tetrodes for each event
        end

        if ~isempty(noncoordriplist)
            allenv = [];
            for ll = 1:length(spintets)
                allenvtmp = [];
                tetrode = spintets(ll);
                if tetrode <10
                    tetstring = ['0',num2str(tetrode)];
                else
                    tetstring = num2str(tetrode);
                end
                currspinfile = [dir,'/EEG/',animalprefix,'spindlegnd', daystring,'-',epochstring,'-',tetstring];
                load(currspinfile);
                env = zscore(double(spindlegnd{day}{epoch}{tetrode}.data(:,3)));
                for l = 1:length(noncoordriplist(:,1))
                    riptimetmp = [noncoordriplist(l,1) noncoordriplist(l,2)];
                    tmp_idx = find(times >= riptimetmp(1) & times <= riptimetmp(2));
                    envtmp = max(env(tmp_idx));
                    allenvtmp = [allenvtmp envtmp];
                end
                allenv = [allenv; allenvtmp];
            end
            allenv_noncoord = [allenv_noncoord; mean(allenv)'];
        end
        clear coordriplist noncoordriplist
    end
end
meancoord = mean(allenv_coord);
meannoncoord = mean(allenv_noncoord);
semcoord = std(allenv_coord)./sqrt(length(allenv_coord));
semnoncoord = std(allenv_noncoord)./sqrt(length(allenv_noncoord));
bar([meancoord meannoncoord]);
hold on
errorbar([meancoord meannoncoord], [semcoord semnoncoord],'LineStyle','none')
xticklabels({'Coord PFC Rips','NonCoord PFC Rips'})
ylabel('Spindle Power (zscore)')
keyboard

