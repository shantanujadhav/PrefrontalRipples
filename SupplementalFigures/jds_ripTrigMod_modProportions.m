%%------------------------------------------------------------------------
%Justin D. Shin

%Plots the proportions of modulated cells for different ripple types
%%------------------------------------------------------------------------
clear all
close all

scrsz = get(0,'ScreenSize');
savedir = '/Volumes/JUSTIN/SingleDay/ProcessedDataNew/';
area='PFC';
riptype = {'noncoord','Coord','All'}; %ripple types to evaluate

savefig1=0;
set(0,'defaultaxesfontsize',16);

allallripmodhists = {};
allallripmodhistssig = {};
allallripmodhistssigInd=  {};
allallripmodhiststype = {};
allallripmodpval = {};
allallripmodMI = {};

b=gaussian(20,61); %smoothing
totalCells = [];
for ripInd=1:length(riptype)
    rtype=riptype{ripInd};

    load([savedir sprintf('Allanim_%s250ca1ripplemod_by250mscrit_sleep_%s_alldata_largewin_sepeps_gather_X6.mat',rtype,area)])
    allripplemod_idx=[];
    for w=1:length(allripplemod)
        allripplemod_idx=[allripplemod_idx;allripplemod(w).index];
    end

    allripmodhists=[];
    allripmodhistssig=[];
    allripmodhistssigInd=[];
    allripmodonset3=[];
    allripmodhiststype=[];
    allripmodpvals = [];
    allripmodMI = [];
    totalCells(ripInd) = length(allripplemod);
    for i=1:length(allripplemod)
        if allripplemod(i).rasterShufP2<0.05 %get all significant
            allripmodhists=[allripmodhists; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
            curranim = allripplemod(i).index(1);
            allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
            allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
            allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
            allripmodpvals = [allripmodpvals; allripplemod(i).rasterShufP2];
            allripmodMI = [allripmodMI; allripplemod(i).Dm];
        end
    end

    allallripmodhists{ripInd}=allripmodhists;
    allallripmodhistssig{ripInd}=allripmodhistssig;
    allallripmodhistssigInd{ripInd}=allripmodhistssigInd;
    allallripmodhiststype{ripInd}=allripmodhiststype;
    allallripmodpval{ripInd} = allripmodpvals;
    allallripmodMI{ripInd} = allripmodMI;

    cntcellsRip=length(allripplemod);
end
ncpsths=allallripmodhistssig{1};
ncinds=allallripmodhistssigInd{1};
ncpvals = allallripmodpval{1};
ncMI = allallripmodMI{1};
allpsths=allallripmodhistssig{2};
allinds=allallripmodhistssigInd{2};
allpvals = allallripmodpval{2};
allMI = allallripmodMI{2};

figure;
propData = [length(allallripmodhiststype{1})/totalCells(1) length(allallripmodhiststype{2})/totalCells(2) length(allallripmodhiststype{3})/totalCells(3); ...
    length(find(allallripmodhiststype{1} == 0))/length(allallripmodhiststype{1}) length(find(allallripmodhiststype{2} == 0))/length(allallripmodhiststype{2}) length(find(allallripmodhiststype{3} == 0))/length(allallripmodhiststype{3}); ...
    length(find(allallripmodhiststype{1} == 1))/length(allallripmodhiststype{1}) length(find(allallripmodhiststype{2} == 1))/length(allallripmodhiststype{2}) length(find(allallripmodhiststype{3} == 1))/length(allallripmodhiststype{3})];

bar(propData)
xticklabels({'Prop. Modulated', 'Prop. Inh', 'Prop. Exc'})
title(sprintf('%s - NC vs Coord vs All mod',area))
ylabel('Proportion')
set(gcf, 'renderer', 'painters')
pvs = [];
chiStats = [];

%get sig prop diff for two pops of overall mod prop
for m = 1:3
    for n = 1:3
        if (m ~= n) && (n > m)
            n1 = length(allallripmodhiststype{m}); N1 = totalCells(m);
            n2 = length(allallripmodhiststype{n}); N2 = totalCells(n);
            x1 = [repmat('a',N1,1); repmat('b',N2,1)];
            x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
            [tbl,chi2stat,pval] = crosstab(x1,x2)
            pvs = [pvs; pval];
            chiStats = [chiStats; chi2stat];
        end
    end
end

for i = 1:2
    mod = i-1;
    for m = 1:3
        for n = 1:3
            if (m ~= n) && (n > m)
                n1 = length(find(allallripmodhiststype{m} == mod)); N1 = length(allallripmodhiststype{m});
                n2 = length(find(allallripmodhiststype{n} == mod)); N2 = length(allallripmodhiststype{n});
                x1 = [repmat('a',N1,1); repmat('b',N2,1)];
                x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
                [tbl,chi2stat,pval] = crosstab(x1,x2)
                pvs = [pvs; pval];
                chiStats = [chiStats; chi2stat];
            end
        end
    end
end

keyboard