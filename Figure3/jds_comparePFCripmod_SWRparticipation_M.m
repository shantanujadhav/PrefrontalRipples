%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates CA1 mod cell firing rates, SWR gain, and prct ripples active.
%%------------------------------------------------------------------------
clear all;
close all;
savedirX = '/Volumes/JUSTIN/SingleDay/ProcessedData/';

animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8'};
day = 1;
epochs = [1:2:17];
excluderips = 1;
exc_SWRfiring = [];
inh_SWRfiring = [];
nonsig_SWRfiring = [];
inh_FR = [];
exc_FR = [];
inh_FRgain = [];
exc_FRgain = [];
nonsig_FRgain = [];
nonsig_FR = [];
inh_FRsws = [];
exc_FRsws = [];
nonsig_FRsws = [];
inh_FRrip = [];
exc_FRrip = [];
nonsig_FRrip = [];
for a = 1:length(animalprefixlist)
    animalprefix = char(animalprefixlist{a});
    animdir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    load(sprintf('%s%srippletime_coordSWS0%d.mat',animdir,animalprefix,day));% get ripple time
    load(sprintf('%s%sspikes0%d.mat',animdir,animalprefix,day));
    load(sprintf('%s%sswsALL0%d.mat',animdir,animalprefix,day));
    load(sprintf('%s%sCA1ctxripmodsig_epsExcludeHigh0%d.mat',animdir,animalprefix,day));
       
    sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
    
    for e = 1:length(epochs)
        epoch = epochs(e);
        
        ep2 = find(sleeps(:,2) == epoch);
        modcells = epochModulation.cellidx;
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = -1;
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = 1;
        
        allcells = [inhcells; exccells];
        
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end
        
        riptimes = [ripple{day}{epoch}.starttime ripple{day}{epoch}.endtime];
        swslist = [sws{day}{epoch}.starttime sws{day}{epoch}.endtime];
        swsdur = sws{day}{epoch}.total_duration;
        curreegfile = [animdir,'/EEG/',animalprefix,'eeg', '01','-',epochstring,'-','02']; %use any tetrode
        load(curreegfile);
        time1 = geteegtimes(eeg{day}{epoch}{2}) ; % construct time array
        
        [~,swsvec] = wb_list2vec(swslist,time1);
        
        [~,ripvec] = wb_list2vec(riptimes,time1);
        %% If want to exclude ripples from SWS FR calculation
        if excluderips == 1
        swsvec_new = swsvec & ~ripvec;
        
        swslist = vec2list(swsvec_new,time1);
        end
        %%
        swsdur = sum(swslist(:,2) - swslist(:,1));
        ripdur = sum(riptimes(:,2) - riptimes(:,1));
        for cellcount = 1:length(allcells(:,1))
            cellmod = allcells(cellcount,3);
          
            if ~isempty(riptimes)
                index = [day,epoch,allcells(cellcount,[1:2])];
                
                if (length(riptimes(:,1)) > 1) && (length(swslist(:,1)) > 1)
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                        fRate = spikes{index(1)}{index(2)}{index(3)}{index(4)}.meanrate;
                        spikebins = periodAssign(spiketimes, riptimes);
                        
                        swsbins = periodAssign(spiketimes, swslist);
                        inSwsSpk = length(find(swsbins ~= 0));
                        fRateSws = inSwsSpk/swsdur;
                        
                        numSpkRip = length(find(spikebins ~= 0));
                        ripFR = numSpkRip/ripdur;
                        
                        inRipSpk = length(find(unique(spikebins) ~= 0));
                        pctRipActivity = (inRipSpk/(length(riptimes(:,1))));
                        if cellmod == 1
                            exc_SWRfiring = [exc_SWRfiring; pctRipActivity];
                            exc_FR = [exc_FR; fRate];
                            exc_FRsws = [exc_FRsws; fRateSws];
                            exc_FRrip = [exc_FRrip; ripFR];
                            exc_FRgain = [exc_FRgain; (ripFR/fRateSws)];
                        elseif cellmod == -1
                            inh_SWRfiring = [inh_SWRfiring; pctRipActivity];
                            inh_FR = [inh_FR; fRate];
                            inh_FRsws = [inh_FRsws; fRateSws];
                            inh_FRrip = [inh_FRrip; ripFR];
                            inh_FRgain = [inh_FRgain; (ripFR/fRateSws)];
                        end
                    end
                end
            end
        end
    end
end
datameans = [mean(exc_SWRfiring) mean(inh_SWRfiring) ];
datasems = [(std(exc_SWRfiring)/sqrt(length(exc_SWRfiring)))...
    (std(inh_SWRfiring)/sqrt(length(inh_SWRfiring)))];

bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Percent PFC Ripples Active')
title('PFC Ripple Modulation and PFC Ripple Participation')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p h] = ranksum(exc_SWRfiring,inh_SWRfiring)
%%
datameansFR = [mean(exc_FR) mean(inh_FR)];
datasemsFR = [(std(exc_FR)/sqrt(length(exc_FR)))...
    (std(inh_FR)/sqrt(length(inh_FR)))];

figure
bar([1:2],datameansFR,'k')
hold on
er = errorbar([1:2],datameansFR,datasemsFR);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Baseline Firing Rate (Hz)')
title('CA1 Cell Mean Firing Rate')
xticklabels({'CA1exc','CA1inh','NonsigMod'}); xtickangle(45)

[p1 h1] = ranksum(exc_FR,inh_FR)

%%
datameansFRsws = [mean(exc_FRsws) mean(inh_FRsws)];
datasemsFRsws = [(std(exc_FRsws)/sqrt(length(exc_FRsws)))...
    (std(inh_FRsws)/sqrt(length(inh_FRsws)))];

figure
bar([1:2],datameansFRsws,'k')
hold on
er = errorbar([1:2],datameansFRsws,datasemsFRsws);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Mean Firing Rate SWS (Hz)')
title('CA1 Cell Mean Firing Rate During SWS')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p2 h2] = ranksum(exc_FRsws,inh_FRsws)

figure;
datacombined = [exc_FRsws; inh_FRsws];
g1 = repmat({'CA1exc'},length(exc_FRsws),1);
g2 = repmat({'CA1inh'},length(inh_FRsws),1);
g = [g1;g2];

boxplot(datacombined,g);
title(['SWSFR-p = ' num2str(p2)])


%%
datameansFRrip = [mean(exc_FRrip) mean(inh_FRrip)];
datasemsFRrip = [(std(exc_FRrip)/sqrt(length(exc_FRrip)))...
    (std(inh_FRrip)/sqrt(length(inh_FRrip)))];

figure
bar([1:2],datameansFRrip,'k')
hold on
er = errorbar([1:2],datameansFRrip,datasemsFRrip);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Firing Rate During SWRs (Hz)')
title('CA1 Cell Firing Rate During Sleep PFC Ripples')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p3 h3] = ranksum(exc_FRrip,inh_FRrip)


%%

datameansFRgain = [mean(exc_FRgain) mean(inh_FRgain)];
datasemsFRgain = [(std(exc_FRgain)/sqrt(length(exc_FRgain)))...
    (std(inh_FRgain)/sqrt(length(inh_FRgain)))];

figure
bar([1:2],datameansFRgain,'k')
hold on
er = errorbar([1:2],datameansFRgain,datasemsFRgain);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('SWR Firing Rate Gain (xBaseline Hz)')
title('CA1 Cell Firing Rate Gain - PFC Ripples')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)

[p4 h4] = ranksum(exc_FRgain,inh_FRgain)


figure;
datacombinedGain = [exc_FRgain; inh_FRgain];
g1 = repmat({'CA1exc'},length(exc_FRgain),1);
g2 = repmat({'CA1inh'},length(inh_FRgain),1);
g = [g1;g2];

boxplot(datacombinedGain,g);
title(['FR Gain-p = ' num2str(p4)])

keyboard;