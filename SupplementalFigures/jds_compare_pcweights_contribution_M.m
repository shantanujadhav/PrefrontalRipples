%%------------------------------------------------------------------------
%Justin D. Shin

%Plots the relationship absolute weight from PCA/ICA analysis and
%contribution metric from Peyrache 2009.
%%------------------------------------------------------------------------
clear all;
close all;
%%
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8'};
day = 1;
sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
weights_contribution = [];
all = [];
exc = [];
inh = [];
excMember = [];
inhMember = [];
memberCells = [];
memberCellsNeg = [];
for a = 1:length(animalprefixlist)
    animalprefix = char(animalprefixlist{a});
    animdir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);


    load(sprintf('%s%sCA1icareactivationcontributionSpk_20_SWS0%d.mat',animdir,animalprefix,day));
    load(sprintf('%s%sCA1_icareactivationtimes20SWSSpk0%d.mat',animdir,animalprefix,day));
    load(sprintf('%s%sCA1ctxripallmodNewWin_epsIncludeHigh0%d.mat',animdir,animalprefix,day));

    epochs = find(~cellfun(@isempty,cell_contribution));
    ep19idx = find(epochs == 19);
    epochs(ep19idx) = [];

    for e = 1:length(epochs)
        ep = epochs(e);
        modcells = epochModulation.cellidx;
        eps = [epochs(e) epochs(e)];

        epsleep = eps(2);

        ep2 = find(sleeps(:,2) == epsleep);
        modvals = epochModulation.modVals(:,ep2);
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = modvals(find(epochModulation.modMat(:,ep2) == -1));
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = modvals(find(epochModulation.modMat(:,ep2) == 1));

        allmodcells = [inhcells; exccells];
        for ii = 1:length(cell_contribution{ep}.contribution)
            weightcontrib_tmp = [];
            thiscell = cell_contribution{ep}.contribution{ii}.contrib_SWS;
            thiscellidx = cell_contribution{ep}.contribution{ii}.cellidx;
            [q, idxloc] = ismember(thiscellidx, cell_contribution{ep}.cellWeights{1}(:,[1 2]), 'rows');
            [q, idxloc2] = ismember(thiscellidx, allmodcells(:,[1 2]), 'rows');
            for l = 1:length(cell_contribution{ep}.cellWeights)
                meanW = mean(cell_contribution{ep}.cellWeights{l}(:,3));
                stdW = std(cell_contribution{ep}.cellWeights{l}(:,3));
                thresh = meanW + 2*stdW;
                weight_tmp = abs(cell_contribution{ep}.cellWeights{l}(idxloc,3));
                weight_tmp2 = cell_contribution{ep}.cellWeights{l}(idxloc,3);
                contrib_tmp = thiscell{l};
                if weight_tmp2 > thresh
                    memberCells = [memberCells; weight_tmp contrib_tmp];
                else
                    weightcontrib_tmp = [weightcontrib_tmp; weight_tmp contrib_tmp];
                end
            end
            all = [all; weightcontrib_tmp];
            if idxloc2 ~= 0
                if allmodcells(idxloc2,3) > 0
                    exc = [exc; weightcontrib_tmp];
                elseif allmodcells(idxloc2,3) < 0
                    inh = [inh; weightcontrib_tmp];
                end
            else
                weights_contribution = [weights_contribution; weightcontrib_tmp];
            end
        end
    end
end

[r, p] = corrcoef(weights_contribution);
scatter(weights_contribution(:,1),weights_contribution(:,2))
hold on
l=lsline
l.LineWidth = 3; l.Color = 'r';

figure
f = fit(weights_contribution(:,1),weights_contribution(:,2),'poly2');
f2 = fit(exc(:,1),exc(:,2),'poly2');
f3 = fit(inh(:,1),inh(:,2),'poly2');
f4 = fit(all(:,1),all(:,2),'poly2');
figure;
plot(f4,all(:,1),all(:,2),'.k')
hold on
plot(memberCells(:,1),memberCells(:,2),'bo')
f4(1,1).LineWidth = 4; f4(1,1).Color = 'r';

plot(f,weights_contribution(:,1),weights_contribution(:,2),'.k')
hold on
plot(f,exc(:,1),exc(:,2),'ro')
plot(f,inh(:,1),inh(:,2),'bo')
f2(1,1).LineWidth = 4; f2(1,1).Color = 'm';
keyboard