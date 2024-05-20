function jds_extractSlowOscillations_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Extracts slow oscillations and categorize as SO (high amplitude events) or
%other/delta (lower amplitude events)
%%------------------------------------------------------------------------

day = 1;
epochs = [1:2:17];
daystring = sprintf('%02d',day);
savedata = 1;
%load SO filter
load('/Users/justinshin/Desktop/Code/usrlocal/filtering/sofilter.mat');

compilewavesSWSSlow = [];
compilewavesSWSDelta = [];
deltaISI = [];
slowISI = [];
deltaPeaks = [];
slowPeaks = [];
slowDownMua = [];
slowUpMua = [];
deltaTroughs = [];
slowTroughs = [];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    load(sprintf('%s%s_spikematrix_ev_singleepochallcellints10_%02d.mat',dir,animalprefix,day));

    tets = tetinfo{1}{epochs(1)};

    ctxtets = []; %get all ctxriptet tetrodes
    for t = 1:length(tets)
        tmp = tets{t};
        if isfield(tmp, 'descrip')
            if isequal(tmp.descrip, 'ctxriptet')
                ctxtets = [ctxtets; t];
            end
        end
    end

    for e = 1:length(epochs)
        compilewavesAllSlow = [];
        compilewavesAllDelta = [];
        lfpdataall = [];
        epoch = epochs(e);

        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end

        for i = 1:length(ctxtets)
            ctxtet = ctxtets(i);

            if (ctxtet<10)
                ctxtetstring = ['0',num2str(ctxtet)];
            else
                ctxtetstring = num2str(ctxtet);
            end

            curreegfile = [dir,'/EEG/',animalprefix,'eeg', daystring,'-',epochstring,'-',ctxtetstring];
            load(curreegfile);

            lfpdata = eeg{day}{epoch}{ctxtet}.data*-1; %re-flip lfp to match original - negative is up state

            lfpdataall = [lfpdataall; lfpdata'];
        end
       
        [bhigh,ahigh] = butter(2,0.1/(1500/2),'high');
        [blow,alow] = butter(5,4/(1500/2),'low');
%         [B A] = butter(2,[0.1/(1500/2) 4/(1500/2)]);
        ampData = filtfilt(bhigh,ahigh,mean(lfpdataall));
        ampData = filtfilt(blow,alow,ampData);
%         ampData = filtfilt(B,A,mean(lfpdataall));
        hdata = hilbert(ampData);
        env = abs(hdata);
        phase = angle(hdata);
        phasedata = int16(phase*10000);
        times = geteegtimes(eeg{day}{epoch}{ctxtet}); % construct time array

        %Zscore amplitude data and find epochs where peak is >=-2
        z_amp = zscore(ampData);
        zci = find(diff(sign(z_amp))); %find the zero crossing indices
        startidx = zci(1:end-1); %look at pretty much every interval
        endidx = zci(2:end);
        indices = [startidx' endidx'];

        ctxData = mean(observation_matrix{epoch}.ctxdata);
        ctxDataZ = zscore(mean(observation_matrix{epoch}.ctxdata));
        prctileMUAhigh = prctile(mean(observation_matrix{epoch}.ctxdata),70);
        prctileMUAlow = prctile(mean(observation_matrix{epoch}.ctxdata),30);
        muaTimes = observation_matrix{epoch}.timeeeg(1):0.01:observation_matrix{epoch}.timeeeg(end);

        overthresh_idx = [];
        trough_idx = [];
        peak_idx = [];
        all_peaks = [];
        all_troughs = [];
        subset_peaks = [];
        subset_troughs = [];
        stendidx = [];
        d2u = [];
        u2d = [];
        muaDown = [];
        muaUp = [];
        compilewaves = [];

        for l = 1:length(indices(:,1))-2 %look at two intervals at a time
            idx = indices(l,:);
            if z_amp(idx(1) + 1) < 0 %only candidate if first crossing is neg to pos
                continue
            end
            idx2 = indices(l+2,:);
            stendidxtmp = [idx(1) idx2(2)];
            ampvector = z_amp(idx(1):idx2(2)); %index 1:4
            [max_z maxidx] = max(z_amp(idx(1):idx(2))); %preceeding up state trough
            [min_z minidx] = min(z_amp(idx(2):idx2(1))); %preceeding down state peak
            [max_zPost maxidxPost] = max(z_amp(idx2(1):idx2(2))); %post down state peak
            downLength = (idx(2) - idx(1))/1500; %down time in ms
            upLength = (idx2(2) - idx2(1))/1500; %up time in ms
            all_peaks = [all_peaks; max_z]; %to get pctile
            all_troughs = [all_troughs; min_z]; %to get pctile
            p2t = minidx + (length(z_amp(idx(1):idx(2))) - maxidx); %number of data points
            p2t = p2t/1500; %convert to ms
            alignIdx = (length(z_amp(idx(1):idx(2))) + minidx + idx(1));
            downPeak = maxidx + idx(1);
            muaIdx = [lookup(times(idx(1)), muaTimes) lookup(times(idx(2)), muaTimes)];
            muaIdx2 = [lookup(times(idx2(1)), muaTimes) lookup(times(idx2(2)), muaTimes)];

            if (downLength > 0.05) && (downLength < 0.5) && (upLength > 0.1) && (upLength < 1) 
                if ((alignIdx-3000) > 0) && ((alignIdx+3000) < length(times))
                    muaDown = [muaDown; mean(ctxDataZ(muaIdx(1):muaIdx(2)))];
                    muaUp = [muaUp; mean(ctxDataZ(muaIdx2(1):muaIdx2(2)))];
                    compilewaves = [compilewaves; z_amp((alignIdx-3000):(alignIdx+3000))]; %2 sec around trough (upstate)
                    stendidx = [stendidx; stendidxtmp];
                    d2u = [d2u; idx(2)]; %zero crossing for down to up
                    u2d = [u2d; idx2(1)]; %zero crossing for up to down
                    subset_peaks = [subset_peaks; max_z];
                    subset_troughs = [subset_troughs; min_z];
                    trough_idx = [trough_idx; alignIdx];
                    peak_idx = [peak_idx; downPeak];
                end
            end
        end
        threshPeaks = prctile(all_peaks, 85); %up
        threshTroughs = prctile(all_troughs, 40); %down
        peakLogicSO = subset_peaks > threshPeaks;
        peakLogicDelta = subset_peaks < threshPeaks;
        troughLogic = subset_troughs < threshTroughs;
        slowLogic = peakLogicSO & troughLogic;
        deltalLogic = peakLogicDelta & troughLogic;
        compilewavesAllSlow = [compilewavesAllSlow; compilewaves(slowLogic,:)];
        compilewavesAllDelta = [compilewavesAllDelta; compilewaves(deltalLogic,:)];
        stendidxSlow = stendidx(slowLogic,:);
        stendidxDelta = stendidx(deltalLogic,:);
        u2dSlowTimes = times(u2d(slowLogic))';
        u2dDeltaTimes = times(u2d(deltalLogic))';
        d2uSlowTimes = times(d2u(slowLogic))';
        d2uDeltaTimes = times(d2u(deltalLogic))';

        slowTimes = times(stendidxSlow);
        deltaTimes = times(stendidxDelta);
        deltaTroughTimes = times(trough_idx(deltalLogic))';
        slowTroughTimes = times(trough_idx(slowLogic))';

        deltaTroughStrength = subset_troughs(deltalLogic);
        slowTroughStrength = subset_troughs(slowLogic);

        deltaPeakTimes = times(peak_idx(deltalLogic))';
        slowPeakTimes = times(peak_idx(slowLogic))';

        deltaPeakStrength = subset_peaks(deltalLogic);
        slowPeakStrength = subset_peaks(slowLogic);

        muaSOdown = muaDown(slowLogic);
        muaSOup = muaUp(slowLogic);

        muaDeltaDown = muaDown(deltalLogic);
        muaDeltaUp = muaUp(deltalLogic);

        %Constrain by SWS
        load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get sws time

        swstime = sws{day}{epoch};
        if ~isempty(swstime.starttime) 
            if length(swstime.starttime) < 2
                midpt = (swstime.endtime - swstime.starttime)/2;
                swslist = [[swstime.starttime(1) swstime.starttime(1)+midpt];...
                    [swstime.starttime(1)+midpt swstime.endtime(1)]];
            else
                swslist = [swstime.starttime swstime.endtime];
            end
            if length(deltaTroughTimes) > 1
                deltaSWS = logical(periodAssign(deltaTroughTimes,swslist)); %peak times that fall in SWS

                compilewavesSWSDelta = [compilewavesSWSDelta; compilewavesAllDelta(deltaSWS, 1500:4500)];

                compilewavesAllDelta = compilewavesAllDelta(deltaSWS, 1500:4500);

                slowosc{day}{epoch}.Deltastarttime = deltaTimes(deltaSWS,1);
                slowosc{day}{epoch}.Deltaendtime = deltaTimes(deltaSWS,2);
                slowosc{day}{epoch}.Deltatroughtime = deltaTroughTimes(deltaSWS);
                slowosc{day}{epoch}.Deltapeaktime = deltaPeakTimes(deltaSWS);
                slowosc{day}{epoch}.Deltatroughsmagnitude = deltaTroughStrength(deltaSWS);
                slowosc{day}{epoch}.Deltapeakmagnitude = deltaPeakStrength(deltaSWS);
                slowosc{day}{epoch}.Deltau2d = u2dDeltaTimes(deltaSWS);
                slowosc{day}{epoch}.Deltad2u = d2uDeltaTimes(deltaSWS);
                slowosc{day}{epoch}.DeltaWavs = compilewavesAllDelta;
                slowosc{day}{epoch}.DeltamuaUp = muaDeltaUp(deltaSWS);
                slowosc{day}{epoch}.DeltamuaDown = muaDeltaDown(deltaSWS);

                deltaPeaks = [deltaPeaks; deltaPeakStrength(deltaSWS)];
                deltaTroughs = [deltaTroughs; deltaTroughStrength(deltaSWS)];
                deltaISI = [deltaISI; diff(deltaTimes(deltaSWS,1))];
            else
                slowosc{day}{epoch}.Deltastarttime = [];
                slowosc{day}{epoch}.Deltaendtime = [];
                slowosc{day}{epoch}.Deltatroughtime = [];
                slowosc{day}{epoch}.Deltapeaktime = [];
                slowosc{day}{epoch}.Deltatroughsmagnitude = [];
                slowosc{day}{epoch}.Deltapeakmagnitude = [];
                slowosc{day}{epoch}.Deltau2d = [];
                slowosc{day}{epoch}.Deltad2u = [];
                slowosc{day}{epoch}.DeltaWavs = [];
                slowosc{day}{epoch}.DeltamuaUp = [];
                slowosc{day}{epoch}.DeltamuaDown = [];
            end
            if length(slowTroughTimes) > 1
                slowSWS = logical(periodAssign(slowTroughTimes,swslist));
                compilewavesSWSSlow = [compilewavesSWSSlow; compilewavesAllSlow(slowSWS, 1500:4500)];
                compilewavesAllSlow = compilewavesAllSlow(slowSWS, 1500:4500);

                slowosc{day}{epoch}.SOstarttime = slowTimes(slowSWS,1);
                slowosc{day}{epoch}.SOendtime = slowTimes(slowSWS,2);
                slowosc{day}{epoch}.SOtroughtime = slowTroughTimes(slowSWS);
                slowosc{day}{epoch}.SOpeaktime = slowPeakTimes(slowSWS);
                slowosc{day}{epoch}.SOtroughsmagnitude = slowTroughStrength(slowSWS);
                slowosc{day}{epoch}.SOpeakmagnitude = slowPeakStrength(slowSWS);
                slowosc{day}{epoch}.SOu2d = u2dSlowTimes(slowSWS);
                slowosc{day}{epoch}.SOd2u = d2uSlowTimes(slowSWS);
                slowosc{day}{epoch}.SOWavs = compilewavesAllSlow;
                slowosc{day}{epoch}.SOmuaUp = muaSOup(slowSWS);
                slowosc{day}{epoch}.SOmuaDown = muaSOdown(slowSWS);
                slowosc{day}{epoch}.phasedata = phasedata;
                slowosc{day}{epoch}.ctxtets = ctxtets;
                slowosc{day}{epoch}.tvec = times;
                slowosc{day}{epoch}.descrip = 'SO and delta waves extracted from mean of ctxriptets, p2t 150-500ms';
                slowosc{day}{epoch}.descrip2 = 'Extracted using method from SO/delta dissociation paper Cell 2019';

                slowPeaks = [slowPeaks; slowPeakStrength(slowSWS)];
                slowTroughs = [slowTroughs; slowTroughStrength(slowSWS)];
                slowISI = [slowISI; diff(slowTimes(slowSWS,1))];
                slowDownMua = [slowDownMua; muaSOdown(slowSWS)];
                slowUpMua = [slowUpMua; muaSOup(slowSWS)];
            else
                slowosc{day}{epoch}.SOstarttime = [];
                slowosc{day}{epoch}.SOendtime = [];
                slowosc{day}{epoch}.SOtroughtime = [];
                slowosc{day}{epoch}.SOpeaktime = [];
                slowosc{day}{epoch}.SOtroughsmagnitude = [];
                slowosc{day}{epoch}.SOpeakmagnitude = [];
                slowosc{day}{epoch}.SOu2d =[];
                slowosc{day}{epoch}.SOd2u = [];
                slowosc{day}{epoch}.SOWavs = [];
                slowosc{day}{epoch}.SOmuaUp = [];
                slowosc{day}{epoch}.SOmuaDown = [];
                slowosc{day}{epoch}.phasedata = phasedata;
                slowosc{day}{epoch}.ctxtets = ctxtets;
                slowosc{day}{epoch}.tvec = times;
                slowosc{day}{epoch}.descrip = 'SO and delta waves extracted from mean of ctxriptets, p2t 150-500ms';
                slowosc{day}{epoch}.descrip2 = 'Extracted using method from SO/delta dissociation paper Cell 2019';
            end
        else
            slowosc{day}{epoch}.SOstarttime = [];
            slowosc{day}{epoch}.SOendtime = [];
            slowosc{day}{epoch}.Deltastarttime = [];
            slowosc{day}{epoch}.Deltaendtime = [];
        end
        clear swslist
    end

    if savedata == 1
        save(sprintf('%s%sslowoscdeltatimesSepSO_SWS%02d.mat', dir,animalprefix,day), 'slowosc');
    end
    clear slowosc
end
bins = 1500;
figure; hold on
pl1 = plot([-bins:bins],mean(compilewavesSWSSlow),'-b','LineWidth',1)
boundedline([-bins:bins],mean(compilewavesSWSSlow),...
    std(compilewavesSWSSlow)./sqrt(size(compilewavesSWSSlow,1)),'-b');
pl2 = plot([-bins:bins],mean(compilewavesSWSDelta),'-m','LineWidth',1)
boundedline([-bins:bins],mean(compilewavesSWSDelta),...
    std(compilewavesSWSDelta)./sqrt(size(compilewavesSWSDelta,1)),'-m');

legend([pl1 pl2],{'Slow oscillations','Delta'})
xlim([-750 750])

figure; hold on
histogram(slowISI);
histogram(deltaISI)
legend({'Slow oscillations','Delta'})

combinedData = [[slowPeaks slowTroughs];[deltaPeaks deltaTroughs]];
[idx,C] = kmeans(combinedData,2);

figure; hold on
plot(combinedData(idx==1,1),combinedData(idx==1,2),'r.','MarkerSize',12)
hold on
plot(combinedData(idx==2,1),combinedData(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
    'MarkerSize',15,'LineWidth',3)

figure; hold on
scatter(slowPeaks,slowTroughs,'b')
scatter(deltaPeaks,deltaTroughs,'r')

keyboard