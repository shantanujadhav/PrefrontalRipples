function jds_sleepripplerate_overtime_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plots ripple rate over time for independent and coordinated events
%%------------------------------------------------------------------------
day = 1;
rvals = [];
for r = 1:3
    allanim_criprate = [];
    allanim_hpriprate = []; %need to better assign variables. confusing
    for a = 1:length(animalprefixlist)
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SD_Control/%s_direct/', animalprefix);
        if r == 1
            load(sprintf('%s%sctxrippletime_coordSWS0%d.mat',dir,animalprefix,day));
            ripple = ctxripple; clear ctxripple
            load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));
        elseif r == 2
            load(sprintf('%s%srippletime_coordSWS0%d.mat',dir,animalprefix,day));
            ctxripple = ripple; clear ripple
            load(sprintf('%s%srippletime_SWS0%d.mat',dir,animalprefix,day));
        elseif r == 3
            load(sprintf('%s%sctxrippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
            load(sprintf('%s%srippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
        end

        load(sprintf('%s%ssws0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%spos0%d.mat',dir,animalprefix,day));

        epochs = [1:2:17];

        cRip_overtime = [];
        hRip_overtime = [];
        ripcnt = 0;
        for ep = 1:length(epochs)

            concat_crip = zeros(1000,1)';
            concat_hrip = zeros(1000,1)';
            epoch = epochs(ep);
            
            swsdur = sws{day}{epoch}.total_duration;

            ctxriptimestmp = [ctxripple{day}{epoch}.starttime ctxripple{day}{epoch}.endtime];
            hpriptimestmp = [ripple{day}{epoch}.starttime ripple{day}{epoch}.endtime];

            if swsdur > 30 % only epochs with sws duration greater than
                c_rate = length(ctxriptimestmp(:,1))/swsdur;
                h_rate = length(hpriptimestmp(:,1))/swsdur;
            else
                c_rate = NaN;
                h_rate = NaN;
            end

            cRip_overtime = [cRip_overtime c_rate];
            hRip_overtime = [hRip_overtime h_rate];
            ripcnt = ripcnt + length(crip_idx);
        end
        % zscore using available epoch data (some entries may be nan)
        zscore_crip = (cRip_overtime - nanmean(cRip_overtime))/nanstd(cRip_overtime);
        zscore_hprip = (hRip_overtime - nanmean(hRip_overtime))/nanstd(hRip_overtime);
        allanim_criprate = [allanim_criprate; zscore_crip];
        allanim_hpriprate = [allanim_hpriprate; zscore_hprip];
    end
    if r == 3 %independent rate over time
        
        tctx = table(allanim_criprate(:,1),allanim_criprate(:,2),...
            allanim_criprate(:,3),allanim_criprate(:,4),allanim_criprate(:,5),...
            allanim_criprate(:,6),allanim_criprate(:,7),allanim_criprate(:,8),...
            allanim_criprate(:,9),'VariableNames',{'t1','t2','t3','t4','t5',...
            't6','t7','t8','t9'});
        rmctx = fitrm(tctx,'t1-t9 ~ 1','WithinDesign',time);
        ranovatblCtx = ranova(rmctx);

        thp = table(allanim_hpriprate(:,1),allanim_hpriprate(:,2),...
            allanim_hpriprate(:,3),allanim_hpriprate(:,4),allanim_hpriprate(:,5),...
            allanim_hpriprate(:,6),allanim_hpriprate(:,7),allanim_hpriprate(:,8),...
            allanim_hpriprate(:,9),'VariableNames',{'t1','t2','t3','t4','t5',...
            't6','t7','t8','t9'});
        rmhp = fitrm(thp,'t1-t9 ~ 1','WithinDesign',time);
        ranovatblHp = ranova(rmhp);
    end

    figure; hold on; title('Sleep PFC Ripple Rate During SWS')
    plot(allanim_criprate','LineWidth',2)
    plot(nanmean(allanim_criprate),'-k','LineWidth',4)
    sem_ctxrip = nanstd(allanim_criprate)./sqrt(length(allanim_criprate(:,1)));
    errorbar([1:9],nanmean(allanim_criprate),sem_ctxrip,'-k','LineWidth',4)
    xlim([0.5 9.5]); ylabel('Ripples/sec','FontSize', 16)
    ylim([min(min(allanim_criprate))-0.1 max(max(allanim_criprate))+0.1]);
    a1 = gca; a1.FontSize = 16;

    figure; hold on; title('Sleep CA1 Ripple Rate During SWS')
    plot(allanim_hpriprate','LineWidth',2)
    plot(nanmean(allanim_hpriprate),'-k','LineWidth',4)
    sem_hprip = nanstd(allanim_hpriprate)./sqrt(length(allanim_hpriprate(:,1)));
    errorbar([1:9],nanmean(allanim_hpriprate),sem_hprip,'-k','LineWidth',4)
    xlim([0.5 9.5]); ylabel('Ripples/sec','FontSize', 16)
    ylim([min(min(allanim_hpriprate))-0.1 max(max(allanim_hpriprate))+0.1]);
    a2 = gca; a2.FontSize = 16;

    figure; hold on;
    p1 = plot(nanmean(allanim_criprate),'-r', 'LineWidth',4)
    errorbar([1:9],nanmean(allanim_criprate),sem_ctxrip,'-r','LineWidth',4)
    p2 = plot(nanmean(allanim_hpriprate),'-k','LineWidth',4)
    errorbar([1:9],nanmean(allanim_hpriprate),sem_hprip,'-k','LineWidth',4)
    legend([p1 p2],{'PFC Ripples','CA1 Ripples'})
    xlim([0.5 9.5]); xlabel('Sleep Session','FontSize', 16); ylabel('Ripples/sec','FontSize', 16)
    a3 = gca; a3.FontSize = 16;

    tmprvals = [];
    % correlation between overall ripple rate and coordinated rate to
    % investigate which area may be driving coordination
    for a = 1:length(allanim_criprate(:,1))
        tmp1 = allanim_criprate(a,:)';
        tmp2 = allanim_hpriprate(a,:)';
        [r p] = corrcoef(tmp1,tmp2,'rows','complete');
        tmprvals = [tmprvals; r(1,2)];
    end
    rvals = [rvals tmprvals];
end
figure; hold on
[p h] = ranksum(rvals(:,1),rvals(:,2))
rvalSem = std(rvals)./sqrt(length(rvals(:,1)));
bar(1:2, mean(rvals(:,[1 2])),'k')
errorbar(1:2, mean(rvals(:,[1 2])),rvalSem(1:2),'k','LineStyle','none')
xticks([1 2])
xticklabels({'ctx-cNc','ca1-cNc'})
title(['Rip rate correlation - p=' num2str(p)])
set(gcf, 'renderer', 'painters')

keyboard
