function jds_CA1relativeSuppression_M(animalprefixlist,area,state)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plots and compares CA1 assembly suppression during independent PFC ripples
%during eaarly, middle, and late sleep sessions
%%------------------------------------------------------------------------
day = 1;

bins = 400; %for 20ms bins used for reactivation
peakbins = find(abs(-bins:bins)<=10);

allevents_ctxriptrig_early = [];
allevents_ctxriptrig_middle = [];
allevents_ctxriptrig_late = [];

shuf = 0;
%epoch separation into early, middle, late
epsep = {[3 5 7],[9 11],[13 15 17]};

for sep = 1:length(epsep)
    epochs = epsep{sep};
    for a = 1:length(animalprefixlist)

        animalprefix = char(animalprefixlist(a));
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

        %Load reactivation strength file for all assemblies and epochs
        load(sprintf('%s%s%s_RTimeStrength%sNewSpk_20_%02d.mat',dir,animalprefix,area,state,day));
        %Load ripples

        load(sprintf('%s%sctxrippletime_noncoordSWS%02d.mat',dir,animalprefix,day));
        allevents_ctxriptrig_tmp = [];
        for e = 1:length(epochs)
            ep = epochs(e);
            assemblytmp = RtimeStrength{ep}.reactivationStrength;
            ctxripstarts = ctxripple{day}{ep}.starttime;

            %Do PFC ripples
            if ~isempty(assemblytmp)
                for ii = 1:length(assemblytmp)
                    react_idx = [];
                    for t = 1:length(ctxripstarts)
                        idxtmp = lookup(ctxripstarts(t), assemblytmp{ii}(:,1));
                        react_idx = [react_idx; idxtmp];
                    end
                    atmp = [];
                    strengthstmp = assemblytmp{ii}(:,2);
                    if shuf == 1
                        strengthstmp = strengthstmp(randperm(length(strengthstmp)));
                    end
                    for r = 1:length(react_idx)
                        if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                            tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins)); %get vector of reactivation strenths for specified time period
                            atmp = [atmp; tmp'];
                        end
                    end
                    allevents_ctxriptrig_tmp = [allevents_ctxriptrig_tmp; zscore(mean(atmp,1))];
                end
            end
        end
        mnPeakCtxRip = mean(allevents_ctxriptrig_tmp(:,peakbins),2);
        if sep == 1
            allevents_ctxriptrig_early = [allevents_ctxriptrig_early; mnPeakCtxRip];
        elseif sep == 2
            allevents_ctxriptrig_middle = [allevents_ctxriptrig_middle; mnPeakCtxRip];
        elseif sep == 3
            allevents_ctxriptrig_late = [allevents_ctxriptrig_late; mnPeakCtxRip];
        end
    end
end

n1 = [];
n1(1:length(allevents_ctxriptrig_early)) = 1;
n1 = string(n1)';
n2 = [];
n2(1:length(allevents_ctxriptrig_middle)) = 2;
n2 = string(n2)';
n3 = [];
n3(1:length(allevents_ctxriptrig_late)) = 3;
n3 = string(n3)';

group_data = [allevents_ctxriptrig_early; allevents_ctxriptrig_middle; allevents_ctxriptrig_late];
grouping = [n1;n2;n3]

[p,tbl,stats] = kruskalwallis(group_data,grouping);
[c,m,h,gnames] = multcompare(stats,"CriticalValueType","bonferroni");
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
datacombinedReactivation = [allevents_ctxriptrig_early; allevents_ctxriptrig_middle; allevents_ctxriptrig_late];
g1 = repmat({'Early'},length(allevents_ctxriptrig_early),1);
g2 = repmat({'Middle'},length(allevents_ctxriptrig_middle),1);
g3 = repmat({'Late'},length(allevents_ctxriptrig_late),1);
g = [g1;g2;g3];

figure; 
h = boxplot(datacombinedReactivation,g); set(h(7,:),'Visible','off');
title('CA1 suppression by phase')
ylabel('Suppression')
set(gcf, 'renderer', 'painters')
keyboard
