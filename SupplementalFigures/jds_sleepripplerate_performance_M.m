function jds_sleepripplerate_performance_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plots the relationship between ripple rate and raw behavioral performance
%calculated from rewarded/unrewarded trials
%%------------------------------------------------------------------------

day = 1;

allanim_criprate_Corr = [];
allanim_hpriprate_Corr = [];
rvalsHp = [];
rvalsCtx = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SD_Control/%s_direct/', animalprefix);

    load(sprintf('%s%sctxrippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%srippletime_noncoordSWS0%d.mat',dir,animalprefix,day));


    load(sprintf('%s%ssws0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%spos0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%strajinfo0%d.mat',dir,animalprefix,day));
    epochs = [1:2:15];

    cRip_Corr = [];
    hRip_Corr = [];
    pctCorr = [];

    for ep = 1:length(epochs)

        epoch = epochs(ep);

        swsdur = sws{day}{epoch}.total_duration;

        ctxriptimestmp = [ctxripple{day}{epoch}.starttime ctxripple{day}{epoch}.endtime];
        hpriptimestmp = [ripple{day}{epoch}.starttime ripple{day}{epoch}.endtime];

        if swsdur > 30
            c_rate = length(ctxriptimestmp(:,1))/swsdur;
            h_rate = length(hpriptimestmp(:,1))/swsdur;
        else
            c_rate = NaN;
            h_rate = NaN;
        end

        pctCorr = [pctCorr; sum(trajinfo{day}{epoch+1}.rewarded)/...
            length(trajinfo{day}{epoch+1}.rewarded)];
        cRip_Corr = [cRip_Corr; c_rate];
        hRip_Corr = [hRip_Corr; h_rate];
    end

    zscore_crip = (cRip_Corr - nanmean(cRip_Corr))/nanstd(cRip_Corr);
    zscore_hprip = (hRip_Corr - nanmean(hRip_Corr))/nanstd(hRip_Corr);
    [rh ph] = corrcoef(pctCorr,zscore_hprip,'rows','complete');
    [rc pc] = corrcoef(pctCorr,zscore_crip,'rows','complete');
    allanim_criprate_Corr = [allanim_criprate_Corr; [zscore_crip pctCorr]];
    allanim_hpriprate_Corr = [allanim_hpriprate_Corr; [zscore_hprip pctCorr]];
    rvalsHp = [rvalsHp; rh(1,2)];
    rvalsCtx = [rvalsCtx; rc(1,2)];
end
rvals = [rvalsCtx rvalsHp];
figure; hold on
[p h] = ranksum(rvals(:,1),rvals(:,2))
rvalSem = std(rvals)./sqrt(length(rvals(:,1)));
bar(1:2, mean(rvals(:,[1 2])),'k')
errorbar(1:2, mean(rvals(:,[1 2])),rvalSem(1:2),'k','LineStyle','none')
xticks([1 2])
xticklabels({'ctx-cNc','ca1-cNc'})
title(['Rip rate - Performance correlation - p=' num2str(p)])
set(gcf, 'renderer', 'painters')

[rh ph] = corrcoef(allanim_hpriprate_Corr,'rows','complete');
[rc pc] = corrcoef(allanim_criprate_Corr,'rows','complete');

figure; hold on
scatter(allanim_hpriprate_Corr(:,1),allanim_hpriprate_Corr(:,2),'bo')
lsline
scatter(allanim_criprate_Corr(:,1),allanim_criprate_Corr(:,2),'ro')
lsline
fontSize = 9; 
title(['Noncoordrip perf (postsleeprun) corr - pPFC-' num2str(pc(1,2)) ' rPFC-' num2str(rc(1,2))...
    ' pCA1-' num2str(ph(1,2)) ' rCA1-' num2str(rh(1,2))], 'FontSize', fontSize)
xlabel('Rip rate (z)')
ylabel('Performance')
set(gcf, 'renderer', 'painters')

% Significance for difference of corr coeffs
za  = 0.5*(log((1+rc(1,2))/(1-rc(1,2))));
zb  = 0.5*(log((1+rh(1,2))/(1-rh(1,2))));
szab = sqrt(1/(length(find(~isnan(allanim_criprate_Corr(:,1))))-3)...
    + 1/(length(find(~isnan(allanim_hpriprate_Corr(:,1))))-3));
z = (za-zb)/szab;
oneTl = 1-normcdf(z, 0, 1);
twoTl = 2*oneTl;

hpSE = sqrt((1-(rh(1,2)*rh(1,2)))/(length(find(~isnan(allanim_hpriprate_Corr(:,1))))-2));
ctxSE = sqrt((1-(rc(1,2)*rc(1,2)))/(length(find(~isnan(allanim_criprate_Corr(:,1))))-2));
figure;
bar([rc(1,2) rh(1,2)],'k')
hold on
errorbar(1:2, [rc(1,2) rh(1,2)], [ctxSE hpSE],'k','LineStyle','none');
title(['Noncoordrip performance corrcoeff Diff - zpval-' num2str(twoTl)])
set(gcf, 'renderer', 'painters')
xticklabels({'PFC','CA1'})
ylabel('Correlation coefficient')

keyboard
