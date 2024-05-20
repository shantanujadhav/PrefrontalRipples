%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates the contribution of single units to activation/reactivation
%strength of assemblies. Metric to compare with cell weights output by
%PCA/ICA reactivation analysis. Method to confirm assembly weight output.
%%------------------------------------------------------------------------
close all;
clear all;
%%
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8'};
excludeRips = 1;
savedata = 1;
duringrips = 0;
duringsws = 0;
PFC = 0;
CA1 = 1;
compileWeightsContribRun = [];
compileWeightsContribSleep = [];
for a = 1:length(animalprefixlist)
    
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    day = 1;
    bin = 20; %ms
    bintime = bin./1000./60; %min
    
    sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
    %%
    load(sprintf('%s%s_spikematrix_ev_allepochallcell20_0%d',dir,animalprefix,day));
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%spos0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%srippletime_ALL0%d.mat',dir,animalprefix,day));
    ripplerun = ripple;
    load(sprintf('%s%srippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxrippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    
    swsdurs = [];
    for t = 1:2:17
        tmp = sws{day}{t}.total_duration;
        swsdurs = [swsdurs tmp];
    end
    
    swsdurs = swsdurs./60;
    idx = find(swsdurs>0); 
    
    %compile epochs to analyze
    eps = [];
    for t = 1:length(idx) %need to modify in case of gaps
        if idx(t) > 1
            sleep = sleeps(idx(t),2);
            run = sleep - 1;
            eps = [eps; run sleep];
        end
    end

    cell_contribution = [];
    for e = 1:length(eps(:,1)) %need to modify in case of gaps
        contrib_tmp = [];
        ep = eps(e,:);
        if PFC == 1
            cellidx = observation_matrix{ep(2)}.ctxidx;
            total_neuron = observation_matrix{ep(2)}.ncortical;
            %---------Sleep spike matrix-------%
            POST_spike = observation_matrix{ep(2)}.ctxdata;
            %---------W1 and W2 spike matrix-------%
            W2_spike = observation_matrix{ep(1)}.ctxdata;
            area = 'PFC';
        elseif CA1 == 1
            cellidx = observation_matrix{ep(2)}.hpidx;
            total_neuron = observation_matrix{ep(2)}.nhp;
            %---------Sleep spike matrix-------%
            POST_spike = observation_matrix{ep(2)}.hpdata;
            %---------W1 and W2 spike matrix-------%
            W2_spike = observation_matrix{ep(1)}.hpdata;
            area = 'CA1';
        end
        %%
        %-----POST SWS-----%
        swsep = sws{day}{ep(2)};
        swslist_POST(:,1) = swsep.starttime;
        swslist_POST(:,2) = swsep.endtime;
        times_filteeg = observation_matrix{ep(2)}.timeeeg;
        timevect = times_filteeg(1)*1000:bin:times_filteeg(end)*1000;
        POST_sws=zeros(size(timevect));
        for sws_seg =1:length(swslist_POST(:,2))
            indtemp = find(timevect >= swslist_POST(sws_seg,1).*1000 & timevect < swslist_POST(sws_seg,2).*1000);
            POST_sws(indtemp) = 1;
        end
        %%
        %--------ripple time------%
        rippost = ctxripple{day}{ep(2)}; %Use PFC ripple times for the sleep sessions
        riplist_POST(:,1) = rippost.starttime;
        riplist_POST(:,2) = rippost.endtime;
        times_filteeg = observation_matrix{ep(2)}.timeeeg;
        timevect = times_filteeg(1)*1000:bin:times_filteeg(end)*1000;
        POST_rip=zeros(size(timevect));
        for rip_seg =1:length(riplist_POST(:,2))
            indtemp = find(timevect >= riplist_POST(rip_seg,1).*1000 & timevect < riplist_POST(rip_seg,2).*1000);
            POST_rip(indtemp) = 1;
        end
        
        st_end_W2 = [pos{day}{ep(1)}.data(1,1) pos{day}{ep(1)}.data(end,1)];
        ripw2 = ripplerun{day}{ep(1)};
        riplist_W2(:,1) = ripw2.starttime;
        riplist_W2(:,2) = ripw2.endtime;
        times_filteeg = observation_matrix{ep(1)}.timeeeg;
        idx_st = lookup(st_end_W2(1),times_filteeg);
        idx_end = lookup(st_end_W2(2),times_filteeg);
        timevect = times_filteeg(idx_st)*1000:bin:times_filteeg(idx_end)*1000;
        W2_rip=zeros(size(timevect));
        for rip_seg =1:length(riplist_W2(:,2))
            indtemp = find(timevect >= riplist_W2(rip_seg,1).*1000 & timevect < riplist_W2(rip_seg,2).*1000);
            W2_rip(indtemp) = 1;
        end
        
        nriplist_W2(:,1) = ripw2.nstarttime;
        nriplist_W2(:,2) = ripw2.nendtime;
        times_filteeg = observation_matrix{ep(1)}.timeeeg;
        timevect = times_filteeg(1)*1000:bin:times_filteeg(end)*1000;
        W2_nrip=zeros(size(timevect));
        for rip_seg =1:length(nriplist_W2(:,2))
            indtemp = find(timevect >= nriplist_W2(rip_seg,1).*1000 & timevect < nriplist_W2(rip_seg,2).*1000);
            W2_nrip(indtemp) = 1;
        end
        
        %%
        qPOST_spike = POST_spike;
        qW2_spike = W2_spike;

        % Exclude ripple time for RUN epochs
        if excludeRips == 1
            idx_nrip_w = find(W2_rip==0);
            W2_spike_raw = W2_spike;
            qW2_spike = qW2_spike(:,idx_nrip_w);
            qW2_spike_all = W2_spike_raw;
        else
            W2_spike_raw = W2_spike;
            qW2_spike_all = W2_spike_raw;
        end

        qW2_spike_z = zscore(qW2_spike')';
        
        qW2_spike_z(isnan(qW2_spike_z)) = 0; %Zscored spike matrix for awake SWRs
        %%
        %-----------correlation matrix----------%
        cW2_spike = corr(qW2_spike_z'); %Correlation matrix using zscore spike count data
        cW2_spike(find(isnan(cW2_spike)))=0;
        
        cPOST_spike =corr(qPOST_spike');
        cPOST_spike(find(isnan(cPOST_spike)))=0;
        
        %%
        %---------------PCA-----------------%
        [u0,s0,v0] = svd(cW2_spike);
        sdiag0 = diag(s0);
        
        %-----PCs significant test-----%
        lamdamax = (1+sqrt(length(qW2_spike(:,1))/length(qW2_spike(1,:)))).^2;
        diagind = find(sdiag0 > lamdamax);
        numPCs = length(diagind);
        
        %ICA
        Psign = u0(:,[1:numPCs]); %Restrict ICA to significant eigenvalues
        Zproj = Psign'*qW2_spike_z;
        [S, H, iter, W] = robustica(Zproj,{});
        V = Psign*W; %Columns of V are the weight vectors of the assembly patterns
        
        %Scale the weight vectors to unit length and process such that highest
        %absolute value is positive (since sign is arbitrary in ICA)
        vtmp = [];
        for t = 1:length(V(1,:))
            w_vectmp = V(:,t);
            w_vec = w_vectmp/norm(w_vectmp);
            min_w = min(w_vec);
            max_w = max(w_vec);
            if abs(min_w) > max_w
                w_vec = w_vec*(-1);
            end
            vtmp = [vtmp w_vec];
        end
        
        V = vtmp;
        %Calculate projection matrices by taking each column of V and taking
        %the outer product of itself
        
        %Use this section to find high weight cells
        ica_neuron_weights = [];
        cW2_spike_ic = cell(1,numPCs);
        for i = 1:numPCs
            tmpmat = V(:,diagind(i))*V(:,diagind(i))';
            tmpmat2 = tmpmat - diag(diag(tmpmat)); %set diagonal to 0
            cW2_spike_ic{i} = tmpmat2;
            
            [A B]=sort(V(:,diagind(i)),'descend');
            tmpcellidx = [cellidx V(:,diagind(i))];
            ica_neuron_weights{i} = tmpcellidx;
        end
        
        thetaratio = sdiag0(diagind)./lamdamax;
        %%
        %------zscore--------%
        qPOST_spike_rip = zscore(POST_spike')';
        
        qW2_spike_rip = zscore(W2_spike_raw')';
        
        qW2_spike_rip(isnan(qW2_spike_rip)) = 0; %Zscored spike matrix for awake SWRs
        qPOST_spike_rip(isnan(qPOST_spike_rip)) = 0; %Zscored spike matrix for POST SLEEP SWRs
       
        %%
        %This will calculate the reactivation strength using all the cells
        R_W2 = zeros(1,length(qW2_spike_rip(1,:)));
        R_POST = zeros(1,length(qPOST_spike_rip(1,:)));
        
        for pc = 1:length(diagind)
            PC = cW2_spike_ic{pc};
            for t = 1:length(qW2_spike_rip(1,:))
                R1_W2_temp = qW2_spike_rip(:,t)'*PC*qW2_spike_rip(:,t);
                R_W2(pc,t) =  R1_W2_temp;
            end
            for t = 1:length(qPOST_spike_rip(1,:))
                R1_POST_temp = qPOST_spike_rip(:,t)'*PC*qPOST_spike_rip(:,t);
                R_POST(pc,t) =  R1_POST_temp;
            end
        end

        %To get the contribution of each cell to the reactivation strength, the
        %reactivation strength is recalculated without the cell in question.
        compareTmpRun = [];
        compareTmpSleep = [];
        for k = 1:total_neuron
            for pc = 1:length(diagind)
                PC = cW2_spike_ic{pc};
                PC(:,k) = []; %delete cell from projection matrix
                PC(k,:) = [];
                qPOST_spike_rip_del = qPOST_spike_rip;
                qW2_spike_rip_del = qW2_spike_rip;
                qPOST_spike_rip_del(k,:) = []; %delete cell from zscored spike matrix
                qW2_spike_rip_del(k,:) = [];
                R_POST_contrib = zeros(1,length(qPOST_spike_rip(1,:)));
                R_Run_contrib = zeros(1,length(qW2_spike_rip(1,:)));

                for t = 1:length(qPOST_spike_rip(1,:))
                    R1_POST_temp = qPOST_spike_rip_del(:,t)'*PC*qPOST_spike_rip_del(:,t);
                    R_POST_contrib(t) =  R1_POST_temp;
                end
                for t = 1:length(qW2_spike_rip(1,:))
                    R1_Run_temp = qW2_spike_rip_del(:,t)'*PC*qW2_spike_rip_del(:,t);
                    R_Run_contrib(t) =  R1_Run_temp;
                end

                tmp1 = mean(abs(R_POST(pc,:)));
                tmp2 = mean(abs(R_POST_contrib));
                Ik = 0.5*(1-(tmp2/tmp1));

                tmp3 = mean(abs(R_W2(pc,:)));
                tmp4 = mean(abs(R_Run_contrib));
                Ik_run  = 0.5*(1-(tmp4/tmp3));

                compileWeightsContribRun = [compileWeightsContribRun; [abs(V(k,pc)) Ik_run]];
                compileWeightsContribSleep = [compileWeightsContribSleep; [abs(V(k,pc)) Ik]];

                compareTmpRun = [compareTmpRun; [abs(V(k,pc)) Ik_run]];
                compareTmpSleep = [compareTmpSleep; [abs(V(k,pc)) Ik]];
                
                contrib_tmp{k}.cellidx = cellidx(k,:);
                contrib_tmp{k}.contrib_SWS{pc} = Ik;
                contrib_tmp{k}.contrib_Run{pc} = Ik_run;
                contrib_tmp{k}.minusCellreactivationSWS{pc} = compareTmpSleep;
                contrib_tmp{k}.minusCellactivationRun{pc} = compareTmpRun;
            end
        end
        cell_contribution{ep(2)}.contribution = contrib_tmp;
        cell_contribution{ep(2)}.reactivationSWS = R_POST;
        cell_contribution{ep(2)}.activationRun = R_W2;
        cell_contribution{ep(2)}.cellWeights = ica_neuron_weights;
        cell_contribution{ep(2)}.comparisonRun = compareTmpRun;
        cell_contribution{ep(2)}.comparisonSleep = compareTmpSleep;

        clear swslist_POST swslist riplist_POST riplist_PRE riplist_W1 riplist_W2 ...
            nriplist_W2 nriplist_W1
    end
    if savedata == 1
        save(sprintf('%s%s%sicareactivationcontributionSpk_20_SWS%02d.mat', dir,animalprefix,area,day), 'cell_contribution', '-v7.3');
        disp(sprintf('%s done',animalprefix))
    end
end
