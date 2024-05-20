function jds_leadlagrippletriggered_assemblystrength_M(animalprefixlist,area,state)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plots CA1 or PFC reactivation strength aligned to leading or lagging
%ripples in the opposing area to show differential responses around leading
%and lagging events.
%%------------------------------------------------------------------------
day = 1;

bins = 400; %for 20ms
peakbins = find(abs(-bins:bins)<=10);

allevents_leadriptrig = [];
allevents_lagriptrig = [];

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%s%s_RTimeStrength%sNewSpk_20_%02d.mat',dir,animalprefix,area,state,day));

    %Load ripples
    load(sprintf('%s%srippletime_leadlag%02d.mat',dir,animalprefix,day));

    %Analyze epochs where there are assemblies detected
    epochs = find(~cellfun(@isempty,RtimeStrength));

    for e = 1:length(epochs)
        ep = epochs(e);
        assemblytmp = RtimeStrength{ep}.reactivationStrength;
        if strcmp(area, 'CA1')
            leadripstarts = ripplecoupling{day}{ep}.hpleadriptimes;
            lagripstarts = ripplecoupling{day}{ep}.hplagriptimes;
        elseif strcmp(area, 'PFC')
            leadripstarts = ripplecoupling{day}{ep}.ctxleadriptimes;
            lagripstarts = ripplecoupling{day}{ep}.ctxlagriptimes;
        end

        %Do PFC ripples
        if (size(lagripstarts,1) > 10) && (size(leadripstarts,1) > 10)
            if ~isempty(assemblytmp)
                for ii = 1:length(assemblytmp)
                    react_idx = [];
                    for t = 1:length(leadripstarts(:,1))
                        idxtmp = lookup(leadripstarts(t,1), assemblytmp{ii}(:,1));
                        react_idx = [react_idx; idxtmp];
                    end
                    atmp = [];
                    for r = 1:length(react_idx)
                        strengthstmp = assemblytmp{ii}(:,2);
                        if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                            tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins)); %get vector of reactivation strenths for specified time period
                            atmp = [atmp; tmp'];
                        end
                    end
                    allevents_leadriptrig = [allevents_leadriptrig; zscore(mean(atmp,1))];
                end
                %Do CA1 ripples
                for ii = 1:length(assemblytmp)
                    react_idx = [];
                    for t = 1:length(lagripstarts(:,1))
                        idxtmp = lookup(lagripstarts(t,1), assemblytmp{ii}(:,1));
                        react_idx = [react_idx; idxtmp];
                    end
                    atmp = [];
                    for r = 1:length(react_idx)
                        strengthstmp = assemblytmp{ii}(:,2);
                        %                 stregthstmp = assemblytmp{ii}(:,4); %for run
                        if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                            tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins));
                            atmp = [atmp; tmp'];
                        end
                    end
                    allevents_lagriptrig = [allevents_lagriptrig; zscore(mean(atmp,1))];
                end
            end
        end
    end
end

allevents_leadriptrig_z_mn = mean(allevents_leadriptrig);
allevents_lagriptrig_z_mn = mean(allevents_lagriptrig);

allevents_leadriptrig_z_sem = (std(allevents_leadriptrig)./sqrt(length(allevents_leadriptrig(:,1))));
allevents_lagriptrig_z_sem = (std(allevents_lagriptrig)./sqrt(length(allevents_lagriptrig(:,1))));

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-bins:bins],allevents_lagriptrig_z_mn,'-k','LineWidth',1)
boundedline([-bins:bins],allevents_lagriptrig_z_mn,allevents_lagriptrig_z_sem,'-k');
pl2 = plot([-bins:bins],allevents_leadriptrig_z_mn,'-r','LineWidth',1)
boundedline([-bins:bins],allevents_leadriptrig_z_mn,allevents_leadriptrig_z_sem,'-r');
title(sprintf('Lead or Lag Ripple Triggered %s Reactivation Strength',area))

ylabel('Reactivation Strength')
xlabel('Time from Ripple Onset (s)')
xlim([(-bins) bins]);
% xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
xticklabels({'-2','-1','0','1','2'});
xx = [0 0];
plot(xx,ax1.YLim,'--k')
ylim = (ax1.YLim);
legend([pl1 pl2],{'Lag ripples','Lead ripples'})
set(gcf, 'renderer', 'painters')

keyboard
