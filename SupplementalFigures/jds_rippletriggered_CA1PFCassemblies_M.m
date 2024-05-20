function jds_rippletriggered_CA1PFCassemblies_M(animalprefixlist,area)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plots CA1-PFC joint reactivation strength aligned to different ripples and
%investigates the relationship to determine if similar coactivity patterns
%exist across different events
%%------------------------------------------------------------------------
day = 1;

bins = 160;
peakbins = find(abs(-bins:bins)<=4);
g1 = gaussian(3, 3);
coordriptrig = [];
noncoordriptrig = [];
coordriptrigall = [];
noncoordriptrigall = [];

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%s%s_RTimeStrengthSleepNewSpk_50_%02d.mat',dir,animalprefix,area,day));
    %Load ripples
    load(sprintf('%s%srippletime_coordSWS%02d.mat',dir,animalprefix,day));
    coordripple = ripple; clear ripple;
    load(sprintf('%s%srippletime_noncoordSWS%02d.mat',dir,animalprefix,day));
    noncoordripple = ripple; clear ripple;

    %Which epochs to analyze
    epochs = find(~cellfun(@isempty,RtimeStrength));

    for e = 1:length(epochs)
        ep = epochs(e);
        assemblytmp = RtimeStrength{ep}.reactivationStrength;
        cross = RtimeStrength{ep}.crossmembers; 
        noncoordripstarts = noncoordripple{day}{ep}.starttime;
        coordripstarts = coordripple{day}{ep}.starttime;
        if (length(coordripstarts) >= 10) && (length(noncoordripstarts) >= 10)
            if ~isempty(assemblytmp)
                for ii = 1:length(assemblytmp)
                    if cross(ii) == 0 %if at least 1 member in each area
                        continue
                    end
                    react_idx = [];
                    for t = 1:length(coordripstarts)
                        idxtmp = lookup(coordripstarts(t), assemblytmp{ii}(:,1));
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
                    react_z = zscore(mean(atmp));
                    react_z = mean(react_z(peakbins));
                    coordriptrig = [coordriptrig; react_z];
                    coordriptrigall = [coordriptrigall; smoothvect(zscore(mean(atmp)),g1)];

                    react_idx = [];
                    for t = 1:length(noncoordripstarts)
                        idxtmp = lookup(noncoordripstarts(t), assemblytmp{ii}(:,1));
                        react_idx = [react_idx; idxtmp];
                    end
                    atmp = [];
                    for r = 1:length(react_idx)
                        strengthstmp = assemblytmp{ii}(:,2);
                        if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                            tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins));
                            atmp = [atmp; tmp'];
                        end
                    end
                    react_z = zscore(mean(atmp));
                    react_z = mean(react_z(peakbins));
                    noncoordriptrig = [noncoordriptrig; react_z];
                    noncoordriptrigall = [noncoordriptrigall; smoothvect(zscore(mean(atmp)),g1)];
                end
            end
        end
    end
end
numAssems = size(noncoordriptrigall,1);
topHalf = floor(numAssems/2);
botHalf = topHalf + 1;

allevents_noncoordriptrig_z_mn = mean(noncoordriptrigall);
allevents_coordriptrig_z_mn = mean(coordriptrigall);

allevents_noncoordriptrig_z_sem = (std(noncoordriptrigall)./sqrt(length(noncoordriptrigall(:,1))));
allevents_coordriptrig_z_sem = (std(coordriptrigall)./sqrt(length(coordriptrigall(:,1))));

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-bins:bins],allevents_coordriptrig_z_mn,'-k','LineWidth',1)
boundedline([-bins:bins],allevents_coordriptrig_z_mn,allevents_coordriptrig_z_sem,'-k');
pl2 = plot([-bins:bins],allevents_noncoordriptrig_z_mn,'-r','LineWidth',1)
boundedline([-bins:bins],allevents_noncoordriptrig_z_mn,allevents_noncoordriptrig_z_sem,'-r');

[S I] = sort(coordriptrig,'descend');

allevents_noncoordriptrigTop_z_mn = mean(noncoordriptrigall(I(1:topHalf),:));
allevents_coordriptrigTop_z_mn = mean(coordriptrigall(I(1:topHalf),:));

allevents_noncoordriptrig_z_sem = (std(noncoordriptrigall(I(1:topHalf),:))./...
    sqrt(length(noncoordriptrigall(I(1:100),:))));
allevents_coordriptrig_z_sem = (std(noncoordriptrigall(I(1:topHalf),:))./...
    sqrt(length(noncoordriptrigall(I(1:100),:))));

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-bins:bins],allevents_noncoordriptrigTop_z_mn,'-k','LineWidth',1)
boundedline([-bins:bins],allevents_noncoordriptrigTop_z_mn,allevents_noncoordriptrig_z_sem,'-k');
pl2 = plot([-bins:bins],allevents_coordriptrigTop_z_mn,'-r','LineWidth',1)
boundedline([-bins:bins],allevents_coordriptrigTop_z_mn,allevents_coordriptrig_z_sem,'-r');

allevents_noncoordriptrigBot_z_mn = mean(noncoordriptrigall(I(botHalf:end),:));
allevents_coordriptrigBot_z_mn = mean(coordriptrigall(I(botHalf:end),:));

allevents_noncoordriptrig_z_sem = (std(noncoordriptrigall(I(botHalf:end),:))./...
    sqrt(length(noncoordriptrigall(I(101:end),:))));
allevents_coordriptrig_z_sem = (std(noncoordriptrigall(I(botHalf:end),:))./...
    sqrt(length(noncoordriptrigall(I(101:end),:))));

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-bins:bins],allevents_noncoordriptrigBot_z_mn,'-k','LineWidth',1)
boundedline([-bins:bins],allevents_noncoordriptrigBot_z_mn,allevents_noncoordriptrig_z_sem,'-k');
pl2 = plot([-bins:bins],allevents_coordriptrigBot_z_mn,'-r','LineWidth',1)
boundedline([-bins:bins],allevents_coordriptrigBot_z_mn,allevents_coordriptrig_z_sem,'-r');

figure
scatter(noncoordriptrig,coordriptrig)
hold on
lsline
[r p] = corrcoef(noncoordriptrig,coordriptrig);

keyboard
