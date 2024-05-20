function jds_rippletriggered_compareStrengthAll_M(animalprefixlist,area)
%%------------------------------------------------------------------------
%Justin D. Shin

%Compare reactivation strengths in CA1 or PFC during different ripple types
%%------------------------------------------------------------------------
day = 1;

bins = 400;
peakbins = find(abs(-bins:bins)<=10);

coordriptrig = [];
noncoordriptrig = [];
otherriptrig = [];

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%s%s_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,area,day));
    %Load ripples
    load(sprintf('%s%sctxrippletime_coordSWS%02d.mat',dir,animalprefix,day));
    coordripple = ctxripple; clear ctxripple;
    load(sprintf('%s%sctxrippletime_noncoordSWS%02d.mat',dir,animalprefix,day));
    noncoordripple = ctxripple; clear ctxripple;
    load(sprintf('%s%srippletime_noncoordSWS%02d.mat',dir,animalprefix,day));
    otherripple = ripple;
    %Which epochs to analyze
    epochs = find(~cellfun(@isempty,RtimeStrength));

    for e = 1:length(epochs)
        ep = epochs(e);
        assemblytmp = RtimeStrength{ep}.reactivationStrength;
        noncoordripstarts = noncoordripple{day}{ep}.starttime;
        coordripstarts = coordripple{day}{ep}.starttime;
        otherripstarts = otherripple{day}{ep}.starttime;
        if (length(coordripstarts) >= 10) && (length(noncoordripstarts) >= 10) && (length(otherripstarts) >= 10)
            if ~isempty(assemblytmp)
                for ii = 1:length(assemblytmp)
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

                    react_idx = [];
                    for t = 1:length(otherripstarts)
                        idxtmp = lookup(otherripstarts(t), assemblytmp{ii}(:,1));
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
                    otherriptrig = [otherriptrig; react_z];
                end
            end
        end
    end
end
[p1 h1] = signrank(noncoordriptrig,coordriptrig);
[p2 h1] = signrank(noncoordriptrig,otherriptrig);
[p3 h1] = signrank(coordriptrig,otherriptrig);

group_data = [noncoordriptrig coordriptrig otherriptrig];

[p,tbl,stats] = kruskalwallis(group_data);
[c,m,h,gnames] = multcompare(stats,"CriticalValueType","bonferroni");
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
datacombinedReactivation = [noncoordriptrig; coordriptrig; otherriptrig];
g1 = repmat({'NonCoord'},length(noncoordriptrig),1);
g2 = repmat({'Coord'},length(coordriptrig),1);
g3 = repmat({'Other'},length(otherriptrig),1);
g = [g1;g2;g3];

figure; 
h = boxplot(datacombinedReactivation,g); set(h(7,:),'Visible','off');
title(sprintf(['%s - Coord-Noncoord-Other Rip Reactivation-p = ' num2str(p1)],area))
set(gcf, 'renderer', 'painters')

[rval pval] = corrcoef(coordriptrig,noncoordriptrig);
figure
scatter(coordriptrig,noncoordriptrig,'.k')
hold on; lsline
set(gcf, 'renderer', 'painters')
ylabel('Coord CA1 ripples')
xlabel('Ind PFC ripples')
title(sprintf(['%s - CoordCA1-IndPFCrip corr - p=' num2str(pval(1,2)) ' r=' num2str(rval(1,2))],area))

keyboard
