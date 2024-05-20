function jds_CA1mod_getFieldMetricsPrePost_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates number of place fields, field width, in/out spike ratio for CA1
%mod cells
%%------------------------------------------------------------------------

%%
FieldsCntExc = [];
FieldsCntInh = [];
FieldsCntNM = [];

RatioExc = [];
RatioInh = [];
RatioNM = [];

FieldWidthExc = []; %If multiple fields, average width of fields
FieldWidthInh = [];
FieldWidthNM = [];

pre = 1;
post = 0;

day = 1;
epochs = [2:2:16];
numshufs = 5000;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    load(sprintf('%s%sCA1ctxripallmodNewWin_epsIncludeHigh0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields

    for e = 1:length(epochs)
        if pre == 1
            eps = [epochs(e) (epochs(e) + 1)];
        elseif post == 1
            eps = [epochs(e) (epochs(e) - 1)];
        end
        epsleep = eps(2);
        [ctxidx, hpidx] =  matchidx_acrossep_singleday(dir, animalprefix, day, eps, []); %(tet, cell)
        hpnum = length(hpidx(:,1));

        sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];

        ep2 = find(sleeps(:,2) == epsleep);
        modcells = epochModulation.cellidx;
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = -1;
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = 1;
        nmcells = modcells(find(epochModulation.modMat(:,ep2) == 0),:);
        nmcells(:,3) = 0;

        allmodcells = [inhcells; exccells; nmcells];

        %%
        if ~isempty(allmodcells)
            for i = 1:length(allmodcells(:,1))
                [b, cellLoc] = ismember(allmodcells(i,[1 2]),hpidx,'rows','legacy');
                if cellLoc ~= 0
                    ca_run = [];
                    cind = allmodcells(i,:);
                    cellmod = cind(3);
                    run_ep = eps(1);
                    NumSpks = [];
                    LinTmp = [];
                    OccTmp = [];
                    SmOccTmp = [];
                    fields = [];
                    for t = 1:4
                        tmpfield1 = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,5); %field
                        pos1 = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,1); %position bin
                        occ1 = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,2); %not smoothed occupancy
                        smocc1 = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,6); %smoothed occ
                        NumSpksTmp = linfields{day}{run_ep}{cind(1)}{cind(2)}{t}(:,3); %num spks per bin

                        stdLength1 = linspace(pos1(1),pos1(end),100);
                        tmp1 = lookup(stdLength1,pos1);
                        tmpfield1std = tmpfield1(tmp1);

                        fields{t} = tmpfield1std;

                        OccTmp = [OccTmp; occ1];
                        SmOccTmp = [SmOccTmp; smocc1];
                        NumSpks = [NumSpks; NumSpksTmp];
                        LinTmp = [LinTmp; tmpfield1];
                    end

                    shuf = [];
                    for s = 1:numshufs
                        sampledData = distsample(sum(NumSpks),OccTmp);
                        spksPerBin = histcounts(sampledData,length(SmOccTmp));

                        shufLinField = gaussSmooth(spksPerBin',2)./SmOccTmp;
                        lowocc = find(SmOccTmp < 2*.1); %very low occupancy bins should be excluded
                        shufLinField(lowocc) = nan;

                        shuf = [shuf; shufLinField'];
                    end

                    fieldVec = zeros(1,length(LinTmp));
                    for l = 1:length(LinTmp)
                        shufTmp = shuf(:,l);
                        pvalue = 1 - sum(shufTmp < LinTmp(l))/length(shufTmp);
                        if pvalue <= 0.01
                            fieldVec(l) = 1;
                        end
                    end

                    binaryVector = fieldVec;
                    % Label each region with a label - an "ID" number.
                    [labeledVector, numRegions] = bwlabel(binaryVector);
                    % Measure lengths of each region and the indexes
                    measurements = regionprops(labeledVector, LinTmp', 'Area', 'PixelValues');
                    % Find regions where the area (length) are 5 or greater and
                    % put the values into a cell of a cell array
                    for k = 1:numRegions
                        if measurements(k).Area >= 5
                            % Width is 5 bins or greater, so store the values.
                            ca_run{k} = measurements(k).PixelValues;
                        end
                    end

                    fieldIdx = [];
                    for f = 1:length(ca_run)
                        tmp = ca_run{f};
                        if ~isempty(tmp)
                            idx = find(labeledVector == f);
                            fieldIdx = [fieldIdx; idx'];
                        end
                    end

                    if ~isempty(ca_run)
                        SpksInField = nanmean(LinTmp(fieldIdx));
                        [C notInIdx] = setdiff(1:length(LinTmp),fieldIdx);
                        SpksOutField = nanmean(LinTmp(notInIdx));

                        tmpWidths = [];
                        NumFields = 0;
                        for w = 1:length(ca_run)
                            tmp = ca_run{w};
                            if tmp ~= 0
                                NumFields = NumFields + 1;
                                tmpW = length(tmp)*2;
                                tmpWidths = [tmpWidths; tmpW];
                            end
                        end
                        if ~isempty(find(tmpWidths < 10))
                            keyboard
                        end

                        FieldWidth = mean(tmpWidths);

                        spkRatio = SpksInField/SpksOutField;
                        beh1_4 = corrcoef(fields{1},fields{4});
                        beh2_3 = corrcoef(fields{2},fields{3});
                        runBeh = [beh1_4(1,2);beh2_3(1,2)];
                        eq1_2 = corrcoef(fields{1},fields{2});
                        eq3_4 = corrcoef(fields{3},fields{4});
                        runEquiv = [eq1_2(1,2); eq3_4(1,2)];
                    else
                        SpksInField = [];
                        SpksOutField = [];

                        NumFields = [];
                        FieldWidth = [];

                        spkRatio = [];

                        runBeh = [];
                        runEquiv = [];
                    end

                    if cellmod == 1
                        FieldsCntExc = [FieldsCntExc; NumFields];
                        RatioExc = [RatioExc; spkRatio];
                        FieldWidthExc = [FieldWidthExc; FieldWidth];
                    elseif cellmod == -1
                        FieldsCntInh = [FieldsCntInh; NumFields];
                        RatioInh = [RatioInh; spkRatio];
                        FieldWidthInh = [FieldWidthInh; FieldWidth];
                    elseif cellmod == 0
                        FieldsCntNM = [FieldsCntNM; NumFields];
                        RatioNM = [RatioNM; spkRatio];
                        FieldWidthNM = [FieldWidthNM; FieldWidth];
                    end
                end
            end
        end
    end
end

%num fields
mnFields = [nanmean(FieldsCntExc) nanmean(FieldsCntInh)];
semFields = [nanstd(FieldsCntExc)./(sqrt(length(FieldsCntExc)))...
    nanstd(FieldsCntInh)./(sqrt(length(FieldsCntInh)))];
figure; hold on
bar(mnFields,'k');
errorbar(1:length(mnFields), mnFields, semFields, '-k', 'LineStyle','none');
title('Number Fields')

[p1 h1] = ranksum(FieldsCntExc,FieldsCntInh)

figure;
datacombinedFields = [FieldsCntExc; FieldsCntInh];
g1 = repmat({'CA1exc'},length(FieldsCntExc),1);
g2 = repmat({'CA1inh'},length(FieldsCntInh),1);
g = [g1;g2];

h = boxplot(datacombinedFields,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['NumFields-p = ' num2str(p1)])
ylim([1 8])
set(gcf, 'renderer', 'painters')

%In/Out Spike Ratio
mnRatio = [nanmean(RatioExc) nanmean(RatioInh)];
semRatio = [nanstd(RatioExc)./(sqrt(length(RatioExc)))...
    nanstd(RatioInh)./(sqrt(length(RatioInh)))];
figure; hold on
bar(mnRatio,'k');
errorbar(1:length(mnRatio), mnRatio, semRatio, '-k', 'LineStyle','none');
title('In-Out Spike Ratio')

%Field Width
mnWidth = [nanmean(FieldWidthExc) nanmean(FieldWidthInh)];
semWidth = [nanstd(FieldWidthExc)./(sqrt(length(FieldWidthExc)))...
    nanstd(FieldWidthInh)./(sqrt(length(FieldWidthInh)))];
figure; hold on
bar(mnWidth,'k');
errorbar(1:length(mnWidth), mnWidth, semWidth, '-k', 'LineStyle','none');
title('Field Width')

[p2 h2] = ranksum(FieldWidthExc,FieldWidthInh)

figure;
datacombinedFields = [FieldWidthExc; FieldWidthInh];
g1 = repmat({'CA1exc'},length(FieldWidthExc),1);
g2 = repmat({'CA1inh'},length(FieldWidthInh),1);
g = [g1;g2];

h = boxplot(datacombinedFields,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['Field Width-p = ' num2str(p2)])

keyboard

function  out = gaussSmooth(vector, binrange)

paddinglength = round(binrange*2.5);
padding = ones(paddinglength,1);

out = smoothvect([padding*vector(1) ; vector; padding*vector(end)],gaussian(binrange,binrange*5));
out = out(paddinglength+1:end-paddinglength);


