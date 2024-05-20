function jds_CA1mod_stability_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates and compares field stability for CA1 mod cells
%%------------------------------------------------------------------------

stabilityExc = [];
stabilityInh = [];
coherenceExc = [];
coherenceInh = [];
day = 1;
epochs = [2:2:14];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %%
    %-----match neurons across two epochs-----%
    for e = 1:length(epochs)
        eps = epochs(e):(epochs(e) + 2);
        epsleep = eps(2);
        [ctxidx, hpidx] =  matchidx_acrossep_singleday(dir, animalprefix, day, eps, []); %(tet, cell)
        hpnum = length(hpidx(:,1));

        load(sprintf('%s%sCA1ctxripallmodNewWin_epsIncludeHigh0%d.mat',dir,animalprefix,day));

        sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];

        ep2 = find(sleeps(:,2) == epsleep);
        modcells = epochModulation.cellidx;
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = -1;
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = 1;

        allmodcells = [inhcells; exccells];

        if ~isempty(allmodcells)
            load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
            for i = 1:length(allmodcells(:,1))
                [b, cellLoc] = ismember(allmodcells(i,[1 2]),hpidx,'rows','legacy');
                if cellLoc ~= 0
                    cind = allmodcells(i,:);
                    cellmod = cind(3);
                    pre_ep = eps(1);
                    post_ep = eps(3);
                    preLinTmp = [];
                    postLinTmp = [];
                    for t = 1:4
                        tmp1 = linfields{day}{pre_ep}{cind(1)}{cind(2)}{t}(:,5);
                        tmp2 = linfields{day}{post_ep}{cind(1)}{cind(2)}{t}(:,5);
                        pos1 = linfields{day}{pre_ep}{cind(1)}{cind(2)}{t}(:,1);
                        pos2 = linfields{day}{post_ep}{cind(1)}{cind(2)}{t}(:,1);
                        mi = min(min(pos1),min(pos2));
                        mx = max(max(pos1),max(pos2));
                        md = min(pos1(2)-pos1(1),pos2(2)-pos2(1));
                        xi = mi:md:mx;
                        tmp1 = interp1(pos1,tmp1,xi);
                        tmp2 = interp1(pos2,tmp2,xi);

                        preLinTmp = [preLinTmp; tmp1'];
                        postLinTmp = [postLinTmp; tmp2'];
                    end
                    stability = corrcoef(preLinTmp,postLinTmp,'rows','pairwise');
                    coherence = 0.5*log((1+stability(1,2))/(1-stability(1,2)));

                    if cellmod == 1
                        stabilityExc = [stabilityExc; stability(1,2)];
                        coherenceExc = [coherenceExc; coherence];
                    elseif cellmod == -1
                        stabilityInh = [stabilityInh; stability(1,2)];
                        coherenceInh = [coherenceInh; coherence];
                    end
                end
            end
        end
    end
end

datameans = [mean(stabilityExc) mean(stabilityInh)];
datasems = [(std(stabilityExc)/sqrt(length(stabilityExc)))...
    (std(stabilityInh)/sqrt(length(stabilityInh)))];

datameansCoh = [mean(coherenceExc) mean(coherenceInh)];
datasemsCoh = [(std(coherenceExc)/sqrt(length(coherenceExc)))...
    (std(coherenceInh)/sqrt(length(coherenceInh)))];

figure
bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Stability')

xticklabels({'CA1exc','CA1inh'}); xtickangle(45)
[p1 h1] = ranksum(stabilityExc,stabilityInh)

datacombinedStability = [stabilityExc; stabilityInh];
g1 = repmat({'EXC'},length(stabilityExc),1);
g2 = repmat({'INH'},length(stabilityInh),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedStability,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['CA1 Field Stability-p = ' num2str(p1)])

figure
bar([1:2],datameansCoh,'k')
hold on
er = errorbar([1:2],datameansCoh,datasemsCoh);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Spatial Coherence')
xticklabels({'CA1exc','CA1inh'}); xtickangle(45)
[p1 h1] = ranksum(coherenceExc,coherenceInh)

keyboard
