%%------------------------------------------------------------------------
%Justin D. Shin

%Plots the relationship between PFC ripple properties (freq, amp, length)
%split into quartiles and CA1 INH modulation during independent PFC.
%Gathers modulation data from each quartile file sequentially and plots the
%relationship.
%ripples
%%------------------------------------------------------------------------
close all;
clear all;

datadir = '/Volumes/JUSTIN/SingleDay/ProcessedDataNew/';
ripProperty = 'Amp'; %Amp, Freq, or Length
allData = [];
allDataMean = [];
allDataSem = [];
for q = 1:4
    tmpTrig = [];

    load(sprintf('%sAllanim_noncoord250ctxripplemod400%sQ%d_by250mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6.mat',...
        datadir, ripProperty, q));
    sigCells = vertcat(allripplemod.rasterShufP2) < 0.05;

    cellMod = vertcat(allripplemod.Dm);
    sigModCells = cellMod(sigCells);

    inh_idx = find(sigModCells < 0);
    inhMod = sigModCells(inh_idx);

    allData = [allData; [inhMod ones(length(inhMod),1)*q]]
    allDataMean = [allDataMean mean(inhMod)];
    allDataSem = [allDataSem std(inhMod)/length(sqrt(inhMod))];
end

[r p] = corrcoef(allData);

figure; hold on
plot(1:4, allDataMean,'ko')
errorbar(1:4, allDataMean, allDataSem,'-k','LineStyle','none');
title(['CA1 INH mod - Ind PFC rip corr r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
set(gcf, 'renderer', 'painters')
xlim([0.5 4.5])
xticks([1:4])
lsline

keyboard
