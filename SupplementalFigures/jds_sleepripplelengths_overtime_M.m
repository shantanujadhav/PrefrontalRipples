function jds_sleepripplelengths_overtime_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plot ripple lengths in CA1 and PFC over time
%%------------------------------------------------------------------------

day = 1;
allanim_criplength = [];
allanim_hpriplength = [];
epochs = [1:2:17];
for ep = 1:length(epochs)
    cRip_overtime = [];
    hRip_overtime = [];
    for a = 1:length(animalprefixlist)
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SD_Control/%s_direct/', animalprefix);
        load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%srippletime_SWS0%d.mat',dir,animalprefix,day));

        load(sprintf('%s%spos0%d.mat',dir,animalprefix,day));

        epoch = epochs(ep);
        ctxriptimestmp = [ctxripple{day}{epoch}.starttime ctxripple{day}{epoch}.endtime];
        hpriptimestmp = [ripple{day}{epoch}.starttime ripple{day}{epoch}.endtime];

        if (~isempty(ctxriptimestmp)) && (~isempty(hpriptimestmp))
            c_length = (ctxriptimestmp(:,2)-ctxriptimestmp(:,1));
            h_length = (hpriptimestmp(:,2)-hpriptimestmp(:,1));
            
            %contrain to lengths <400 ms
            c_length = c_length(find(c_length < 0.4));
            h_length = h_length(find(h_length < 0.4));

            cRip_overtime = [cRip_overtime; c_length];
            hRip_overtime = [hRip_overtime; h_length];
        end
    end
    allanim_criplength{ep} = cRip_overtime;
    allanim_hpriplength{ep} = hRip_overtime;
end

n1 = [];
n1(1:length(allanim_criplength{1})) = 1;
n1 = string(n1)';

n2 = [];
n2(1:length(allanim_criplength{2})) = 2;
n2 = string(n2)';

n3 = [];
n3(1:length(allanim_criplength{3})) = 3;
n3 = string(n3)';

n4 = [];
n4(1:length(allanim_criplength{4})) = 4;
n4 = string(n4)';

n5 = [];
n5(1:length(allanim_criplength{5})) = 5;
n5 = string(n5)';

n6 = [];
n6(1:length(allanim_criplength{6})) = 6;
n6 = string(n6)';

n7 = [];
n7(1:length(allanim_criplength{7})) = 7;
n7 = string(n7)';

n8 = [];
n8(1:length(allanim_criplength{8})) = 8;
n8 = string(n8)';

n9 = [];
n9(1:length(allanim_criplength{9})) = 9;
n9 = string(n9)';

group_data = [allanim_criplength{1}; allanim_criplength{2}; allanim_criplength{3};...
    allanim_criplength{4}; allanim_criplength{5}; allanim_criplength{6};...
    allanim_criplength{7}; allanim_criplength{8}; allanim_criplength{9}];
grouping = [n1;n2;n3;n4;n5;n6;n7;n8;n9];

[p,tbl,stats] = kruskalwallis(group_data,grouping);
[c,m,h,gnames] = multcompare(stats,"CriticalValueType","bonferroni");
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

meanCtx = cellfun(@mean,allanim_criplength);
stdCtx = cellfun(@std,allanim_criplength);
numEvsCtx = cellfun(@length,allanim_criplength);
semCtx = stdCtx./sqrt(numEvsCtx);

figure; hold on; title('PFC Ripple Lengths During SWS')
plot(meanCtx,'-k')
errorbar([1:9],meanCtx,semCtx,'-k','LineStyle','none','Marker','o')
xlim([0.5 9.5]); ylabel('PFC ripple lengths (s)','FontSize', 16)
xlabel('Epoch')
set(gcf, 'renderer', 'painters')

figure;
h = boxplot(group_data,grouping); set(h(7,:),'Visible','off');
title('PFC Ripple Lengths During SWS')
ylabel('PFC ripple lengths (s)','FontSize', 16)
ylim([-0.02 0.2])
xlabel('Epoch')
set(gcf, 'renderer', 'painters')
%%

n1 = [];
n1(1:length(allanim_hpriplength{1})) = 1;
n1 = string(n1)';

n2 = [];
n2(1:length(allanim_hpriplength{2})) = 2;
n2 = string(n2)';

n3 = [];
n3(1:length(allanim_hpriplength{3})) = 3;
n3 = string(n3)';

n4 = [];
n4(1:length(allanim_hpriplength{4})) = 4;
n4 = string(n4)';

n5 = [];
n5(1:length(allanim_hpriplength{5})) = 5;
n5 = string(n5)';

n6 = [];
n6(1:length(allanim_hpriplength{6})) = 6;
n6 = string(n6)';

n7 = [];
n7(1:length(allanim_hpriplength{7})) = 7;
n7 = string(n7)';

n8 = [];
n8(1:length(allanim_hpriplength{8})) = 8;
n8 = string(n8)';

n9 = [];
n9(1:length(allanim_hpriplength{9})) = 9;
n9 = string(n9)';

group_data = [allanim_hpriplength{1}; allanim_hpriplength{2}; allanim_hpriplength{3};...
    allanim_hpriplength{4}; allanim_hpriplength{5}; allanim_hpriplength{6};...
    allanim_hpriplength{7}; allanim_hpriplength{8}; allanim_hpriplength{9}];
grouping = [n1;n2;n3;n4;n5;n6;n7;n8;n9];

[p,tbl,stats] = kruskalwallis(group_data,grouping);
[c2,m,h,gnames] = multcompare(stats,"CriticalValueType","bonferroni");
tbl2 = array2table(c2,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

meanHp = cellfun(@mean,allanim_hpriplength);
stdHp = cellfun(@std,allanim_hpriplength);
numEvsHp = cellfun(@length,allanim_hpriplength);
semHp = stdHp./sqrt(numEvsHp);

figure; hold on; title('CA1 Ripple Lengths During SWS')
plot(meanHp,'-k')
errorbar([1:9],meanHp,semHp,'-k','LineStyle','none','Marker','o')
xlim([0.5 9.5]); ylabel('CA1 ripple lengths (s)','FontSize', 16)
xlabel('Epoch')
set(gcf, 'renderer', 'painters')

figure;
h = boxplot(group_data,grouping); set(h(7,:),'Visible','off');
title('CA1 Ripple Lengths During SWS')
ylabel('CA1 ripple lengths (s)','FontSize', 16)
ylim([-0.02 0.35])
xlabel('Epoch')
set(gcf, 'renderer', 'painters')

keyboard
