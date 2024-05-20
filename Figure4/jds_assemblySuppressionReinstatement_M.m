function jds_assemblySuppressionReinstatement_M(animalprefixlist,area,state)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates and plots the relationship between reinstatement (assembly 
%strength comparison between ep(n) and ep(n+1)) and CA1 suppression during \
%independent PFC ripple events
%%------------------------------------------------------------------------
day = 1;

bins = 400; %for 20ms
peakbins = find(abs(-bins:bins)<=10);

reinstatementSupp = [];
AssemCnt = 0;

for a = 1:length(animalprefixlist)
    
    animalprefix = char(animalprefixlist(a));
    
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%s%s_RTimeStrength%sReinstatementSpkZero_20_%02d.mat',dir,animalprefix,area,state,day));
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxrippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    %Which epochs to analyze
    epochs = 2:2:14;
    
    for e = 1:length(epochs)

        ep = epochs(e);
        refepoch = ep;

        ripstarts = ctxripple{day}{ep+1}.starttime;

        if length(ripstarts) < 10
            continue
        end
       
        assemTimes = RtimeStrength{ep}.reactivationStrength{ep+1}{1}(:,1);
        numAssem = length(RtimeStrength{ep}.reactivationStrength{ep+1});
        AssemCnt = AssemCnt + numAssem; 
        react_idx = [];
        for t = 1:length(ripstarts)
            idxtmp = lookup(ripstarts(t), assemTimes);
            react_idx = [react_idx; idxtmp];
        end

        if ~isempty(react_idx)
            for r = 1:numAssem
                rTmp = RtimeStrength{ep}.reactivationStrength{ep+1}{r}; %sleep reactivation between two runs
                
                atmp = [];
                strengthstmp = rTmp(:,2);

                for rr = 1:length(react_idx)
                    if ((react_idx(rr) + bins) < length(strengthstmp)) && ((react_idx(rr) - bins) > 1)
                        tmp = strengthstmp((react_idx(rr) - bins):(react_idx(rr) + bins));
                        atmp = [atmp; tmp'];
                    end
                end
                suppression = zscore(mean(atmp,1));
                suppression = mean(suppression(peakbins));

                assem1 = mean(RtimeStrength{ep}.reactivationStrengthRun{ep}{r}(:,2)); %reference epoch
                assem2 = mean(RtimeStrength{ep}.reactivationStrength{ep+2}{r}(:,2)); %next run epoch
                reinstTmp = assem2-assem1;
                reinstatementSupp = [reinstatementSupp; [reinstTmp suppression]];
            end
        end
    end
end

[r1 p1] = corr(reinstatementSupp,'Type', 'Spearman','Rows','pairwise')
[rPearson pPearson] = corr(reinstatementSupp)
figure
scatter(reinstatementSupp(:,1),reinstatementSupp(:,2),'.k')
hold on; lsline
title(['Reinstatement-Suppression-p=' num2str(p1(1,2))])
xlabel(['Reinstatement-Suppression-r=' num2str(r1(1,2))])

keyboard
