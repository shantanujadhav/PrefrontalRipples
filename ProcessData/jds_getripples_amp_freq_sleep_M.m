%---------------------------------------------------------------%
%  This is the main script for measuring amplitude and frequency%
%  of each SWR event                                            %
%  -- Wenbo Tang (Sep 15, 2019) 
%     Justin Shin (Adapted for Sleep analysis)
%---------------------------------------------------------------%
clc;
close all;
clear all;

%%
animalprefix = ('ZT2'); % animal prefix
animaltestday = 1;% animal experimental day
eps = 1:2:17;% epochs
%%
%---- set parameters ----%
minstd = 3; % min std above mean for ripple detection
minrip = 1; % min number of tetrode detected ripples
minenergy = 0; % min energy threshold
velfilter = 4; %velocity <= 4cm for ripple detection
matcheegtime = 0; % match EEG time?
%%
% set animal directory
animaldir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

eegdir = [animaldir,'EEG/'];% EEG directory
%%
tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo'); % get tetrode info
%%
% day loop

for day = animaltestday    
    
    for id = 1:length(eps)
        d = 1;%day
        e = eps(id);%epoch
        disp(['Animal: ',animalprefix,' Epoch:',num2str(e)])% display current animal, day and epoch
        load(sprintf('%s%ssws0%d.mat',animaldir,animalprefix,d));

        % calculate using riptet only
        tetfilter = 'isequal($descrip, ''riptet'')';% use riptet only (change to ctxriptet for PFC)
        tetlist =  evaluatefilter(tetinfo{d}{e}, tetfilter);
        tetlist = unique(tetlist(:,1))';
        
        % get the mean and std for each frequency channel for z-scoring
        % note that only run it for the first time
        baselinespecgram_forref(animalprefix, d, e, tetlist, 'fpass',[0 400])% high frequency range, 0-400 Hz, for all tets
        sws_ep = sws{d}{e};
        sws_start = sws_ep.starttime;
        sws_end = sws_ep.endtime;
        swstime_all = [sws_start, sws_end];
        
        %Choose ctx ripples or CA1 ripples
        ripples = loaddatastruct(animaldir, animalprefix, 'ripples', d);% load ripple info OR load ctxripples file
%         ripples = loaddatastruct(animaldir, animalprefix, 'ctxripples', d);
        
        r = ripples{d}{e}{tetlist(1)};
        % time range, 10ms bin
        times = r.timerange(1):0.001:r.timerange(end);
        %reset
        nrip = zeros(size(times)); 
        nstd=[];
        ripplestd = zeros(size(times));
        
        % tetrode loop
        for t = 1:length(tetlist)
             tmprip = ripples{d}{e}{tetlist(t)};
             % get the indeces for the ripples with energy above minenergy
             % and maxthresh above minstd
             rvalid = find((tmprip.energy >= minenergy) & (tmprip.maxthresh >= minstd) & (isExcluded(tmprip.midtime, swstime_all)));
             rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
             tmpripplestd = [tmprip.maxthresh(rvalid) tmprip.maxthresh(rvalid)];
             % create another parallel vector with bordering times for zeros
             nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
             rtimes = reshape(rtimes', length(rtimes(:)), 1);
             rtimes(:,2) = 1;
             tmpriplestd = [rtimes(:,1) tmpripplestd(:)];
             nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
                 length(nrtimes(:)), 1) ; r.timerange(2)];
             nrtimes(:,2) = 0;
             % create a new list with all of the times in it
             tlist = sortrows([rtimes ; nrtimes]);
             [junk, ind] = unique(tlist(:,1));
             tlist = tlist(ind,:);

             stdlist = sortrows([tmpriplestd ; nrtimes]);
             stdlist =stdlist(ind,:);
             nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
             nstd(t,:) = interp1(stdlist(:,1), stdlist(:,2), times, 'nearest');  % carry forward amplitude of ripple
        end
        
        %find the start and end borders of each ripple
        inripple = (nrip >= minrip);
        startrippleind = find(diff(inripple) == 1)+1;
        endrippleind = find(diff(inripple) == -1)+1;
        ripplestdout = [];
        
        if ((length(startrippleind) > 1) && (length(endrippleind) > 1))
            if (endrippleind(1) < startrippleind(1))
                endrippleind = endrippleind(2:end);
            end
            if (endrippleind(end) < startrippleind(end))
                startrippleind = startrippleind(1:end-1);
            end
            startripple = times(startrippleind);
            endripple = times(endrippleind);
            %----- measure amplitude of each ripple-----%
            % Get amplitude of "global" ripple: maximum across tetrodes
            [max_nstd,tetid] = max(nstd,[],1);
            ampripple = max_nstd(startrippleind);
            riptet = tetid(startrippleind);

            out = [startripple(:) endripple(:) ampripple(:)]; % amplitude of ripple
            riptet = riptet((out(:,2)-out(:,1))>.050);
            out = out( ((out(:,2)-out(:,1))>.050),:);  % ripples separated by 50 ms
            
            %----- measure frequncy of each ripple-----%
            freqripple = nan(length(riptet),1);
            for r = 1:length(riptet)
                riptime = out(r,1:2);
                %calculate frequency of each ripple event
                freqripple(r) = ripple_frequency_fun(animalprefix, d, e, tetlist(riptet(r)), riptime);
            end
            out = [out,freqripple];
        else
            out = [];
            ripplestdout = [];
        end
        rippleampfreq{d}{e} = out;% save result
        clear out; clear freqripple
    end
    save(sprintf('%s/%sripple_amp_freq_SWS_%02d.mat', animaldir, animalprefix, d), 'rippleampfreq');% save files
    clear rippleampfreq; 
end