%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates reactivation strength for reinstatement analysis
%%------------------------------------------------------------------------
close all;
clear all;
%%
savedata = 1;

PFC = 1;
CA1 = 0;
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8','ER1'};

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    day = 1;
    bin = 20; % ms
    bintime = bin./1000./60; %min
    
    %%
    load(sprintf('%s%s_spikematrix_ev_trackedallepochallcell20_0%d',dir,animalprefix,day));
    load(sprintf('%s%srippletime_ALL0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%spos0%d.mat',dir,animalprefix,day));

    eps = [2 1; 2 3; 2 4; 4 3; 4 5; 4 6; 6 5; 6 7; 6 8; 8 7; 8 9; 8 10;...
        10 9; 10 11; 10 12; 12 11; 12 13; 12 14; 14 13; 14 15; 14 16];

    icareactivationtimes = [];
    RtimeStrength = [];
    
    for e = 1:length(eps(:,1)) %need to modify in case of gaps
        ep = eps(e,:);
        if PFC == 1
            cellidx = observation_matrix{ep(2)}.ctxidx;
            total_ctx_neuron = observation_matrix{ep(2)}.ncortical;
            %---------Sleep spike matrix-------%
            POST_spike = observation_matrix{ep(2)}.ctxdata;
            %---------W1 and W2 spike matrix-------%
            W_spike = observation_matrix{ep(1)}.ctxdata;
            area = 'PFC';
        elseif CA1 == 1
            cellidx = observation_matrix{ep(2)}.hpidx;
            total_ctx_neuron = observation_matrix{ep(2)}.nhp;
            %---------Sleep spike matrix-------%
            POST_spike = observation_matrix{ep(2)}.hpdata;
            %---------W1 and W2 spike matrix-------%
            W_spike = observation_matrix{ep(1)}.hpdata;
            area = 'CA1';
        end
        %%
        %restrict time to posstart and posend
        times_filteeg = observation_matrix{ep(2)}.timeeeg;
        
        st_end_post = [pos{day}{ep(2)}.data(1,1) pos{day}{ep(2)}.data(end,1)];
        
        %-----POST SWS-----%
        idx_st = lookup(st_end_post(1),times_filteeg);
        idx_end = lookup(st_end_post(2),times_filteeg);
        timevect = times_filteeg(idx_st)*1000:bin:times_filteeg(idx_end)*1000;
        %%
        %--------ripple time------%
        rippost = ripple{day}{ep(2)}; %Use PFC ripple times for the sleep sessions
        riplist_POST(:,1) = rippost.starttime;
        riplist_POST(:,2) = rippost.endtime;
        POST_rip=zeros(size(timevect));
        for rip_seg =1:length(riplist_POST(:,2))
            indtemp = find(timevect >= riplist_POST(rip_seg,1).*1000 & timevect < riplist_POST(rip_seg,2).*1000);
            POST_rip(indtemp) = 1;
        end
        
        st_end_W = [pos{day}{ep(1)}.data(1,1) pos{day}{ep(1)}.data(end,1)];
        ripW = ripple{day}{ep(1)};
        riplist_W(:,1) = ripW.starttime;
        riplist_W(:,2) = ripW.endtime;
        times_filteeg = observation_matrix{ep(1)}.timeeeg;
        idx_st = lookup(st_end_W(1),times_filteeg);
        idx_end = lookup(st_end_W(2),times_filteeg);
        timevect = times_filteeg(idx_st)*1000:bin:times_filteeg(idx_end)*1000;
        W_rip=zeros(size(timevect));
        for rip_seg =1:length(riplist_W(:,2))
            indtemp = find(timevect >= riplist_W(rip_seg,1).*1000 & timevect < riplist_W(rip_seg,2).*1000);
            W_rip(indtemp) = 1;
        end
        
        nriplist_W(:,1) = ripW.nstarttime;
        nriplist_W(:,2) = ripW.nendtime;
        times_filteeg = observation_matrix{ep(1)}.timeeeg;
        timevect = times_filteeg(idx_st)*1000:bin:times_filteeg(idx_end)*1000;
        W_nrip=zeros(size(timevect));
        for rip_seg =1:length(nriplist_W(:,2))
            indtemp = find(timevect >= nriplist_W(rip_seg,1).*1000 & timevect < nriplist_W(rip_seg,2).*1000);
            W_nrip(indtemp) = 1;
        end
        
        qPOST_spike = POST_spike;
        qW_spike = W_spike;
        
        % Exclude ripple time for RUN epochs
        idx_nrip_w = find(W_rip==0);
        W_spike_raw = W_spike;
        qW_spike = qW_spike(:,idx_nrip_w);
        %%
        qW_spike_all = W_spike_raw;
        
        qW_spike_z = zscore(qW_spike')';
        
        qW_spike_z(isnan(qW_spike_z)) = 0; %Zscored spike matrix for awake SWRs
        %%
        %-----------correlation matrix----------%
        cW_spike = corr(qW_spike_z'); %Correlation matrix using zscore spike count data
        cW_spike(find(isnan(cW_spike)))=0;
        
        cPOST_spike =corr(qPOST_spike');
        cPOST_spike(find(isnan(cPOST_spike)))=0;
        %%
        %---------------PCA-----------------%
        [u0,s0,v0] = svd(cW_spike);
        sdiag0 = diag(s0);
       
        %-----PCs significant test-----%
        lamdamax = (1+sqrt(length(qW_spike(:,1))/length(qW_spike(1,:)))).^2;
        
        diagind = find(sdiag0 > lamdamax);
        numPCs = length(diagind);
        
        %ICA
        Psign = u0(:,[1:numPCs]); %Restrict ICA to significant eigenvalues
        Zproj = Psign'*qW_spike_z;
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
        cW_spike_ic = cell(1,numPCs);
        for i = 1:numPCs
            tmpmat = V(:,diagind(i))*V(:,diagind(i))';
            tmpmat2 = tmpmat - diag(diag(tmpmat)); %set diagonal to 0
            cW_spike_ic{i} = tmpmat2;
            
            [A B]=sort(V(:,diagind(i)),'descend');
            tmpcellidx = cellidx(B,:);
            tmpcellidx = [tmpcellidx A];
            ica_neuron_weights{i} = tmpcellidx;
        end
        
        thetaratio = sdiag0(diagind)./lamdamax;
        %%
        %------zscore--------%
        %Finding the reactivation strength over the entire epoch
        qPOST_spike_rip = zscore(POST_spike')';
        qW_spike_rip = zscore(W_spike_raw')';
        qW_spike_rip(isnan(qW_spike_rip)) = 0; %Zscored spike matrix for awake SWRs
        qPOST_spike_rip(isnan(qPOST_spike_rip)) = 0; %Zscored spike matrix for POST SLEEP SWRs
        
        %%
        R_W = zeros(1,length(qW_spike_rip(1,:)));
        R_POST = zeros(1,length(qPOST_spike_rip(1,:)));
        
        for pc = 1:length(diagind)
            PC = cW_spike_ic{pc};
            for t = 1:length(qW_spike_rip(1,:))
                R1_W_temp = qW_spike_rip(:,t)'*PC*qW_spike_rip(:,t);
                R_W(pc,t) =  R1_W_temp;
            end
            for t = 1:length(qPOST_spike_rip(1,:))
                R1_POST_temp = qPOST_spike_rip(:,t)'*PC*qPOST_spike_rip(:,t);
                R_POST(pc,t) =  R1_POST_temp;
            end
        end
        
        %%
        R_W(find(R_W < 0)) = 0;
        R_POST(find(R_POST <0)) = 0;

        reactivation_strength_run = R_W';
        trun = (1:length(reactivation_strength_run)).*bin./1000;%s
        reactivation_strength_post = R_POST';
        tpost = (1:length(reactivation_strength_post)).*bin./1000;%s
        
        clr = ['k','b','r','m','g','r','k','b','r','m','g','r','k','b','r',...
            'm','g','r','k','b','r','m','g','r','k','b','r','m','g','r','k','b','r','m','g','r'];
        ylimit = max(max(reactivation_strength_post));
       
        clear riplist_POST riplist_W nriplist_W
        
        times_filteeg_post = observation_matrix{ep(2)}.timeeeg;
        posttimevect = times_filteeg_post(1)*1000:bin:times_filteeg_post(end)*1000;
        
        times_filteeg_w = observation_matrix{ep(1)}.timeeeg;
        wtimevect = times_filteeg_w(1)*1000:bin:times_filteeg_w(end)*1000;
        st_end_W = [pos{day}{ep(1)}.data(1,1) pos{day}{ep(1)}.data(end,1)];
        
        stIdx = lookup(st_end_W(1)*1000,wtimevect);
        endIdx = lookup(st_end_W(2)*1000,wtimevect);
        
        wtimevect = wtimevect(stIdx:endIdx);
        %truncate to get rid of empty bins
        reactivation_strength_run = reactivation_strength_run([stIdx:endIdx],:);
        
        windowWidth = 10; %for 1 second using 100ms bins
        
        rChange = [];
        rvals = [];
        pvals = [];
        for runas = 1:length(reactivation_strength_run(1,:))
            runReactTmp = reactivation_strength_run(:,runas);
            averagedData = zeros(1, floor(length(runReactTmp)/windowWidth));
            binidx = 1:length(averagedData);
            winst = 1;
            for k = 1:floor(length(runReactTmp)/windowWidth)
                averagedData(k) = mean(runReactTmp(winst:winst+windowWidth-1));
                winst = winst + windowWidth;
            end
            [r p] = corrcoef(binidx,averagedData);
            rvals = [rvals; r(1,2)];
            pvals = [pvals; p(1,2)];
            if p(1,2) >= 0.05
                rChange = [rChange; 0];
            elseif (p(1,2) < 0.05) && (r(1,2) > 0)
                rChange = [rChange; 1];
            elseif (p(1,2) < 0.05) && (r(1,2) < 0)
                rChange = [rChange; -1];
            end
        end
        
        react_time_strength_all = [];
        react_time_strength_run = [];
        for ii = 1:length(reactivation_strength_run(1,:))
            
            wpctmp = reactivation_strength_run(:,ii);
            postpctmp = reactivation_strength_post(:,ii);
            
            react_time_strength_all{ii} = [(posttimevect./1000)' postpctmp];
            react_time_strength_run{ii} = [(wtimevect./1000)' wpctmp];
            
            icareactivationtimes{e}.post_strengths{ii} = [(posttimevect./1000)' postpctmp]; %All event strengths, no filter
            icareactivationtimes{e}.run_strengths{ii} = [(wtimevect./1000)' wpctmp];
            icareactivationtimes{e}.runchange = rChange; %increasing, unchanged, decreasing
            icareactivationtimes{e}.runrval = rvals; %correlation rvals
            icareactivationtimes{e}.runpval = pvals; %correlation pvals
            icareactivationtimes{e}.thetaratio = thetaratio;
            icareactivationtimes{e}.pc_weights = ica_neuron_weights;
            icareactivationtimes{e}.cellidx = cellidx;
            icareactivationtimes{e}.timebinsize = bin;
        end
        icareactivationtimes{e}.epochs = ep;
        
        RtimeStrength{ep(1)}.reactivationStrength{ep(2)} = react_time_strength_all;
        RtimeStrength{ep(1)}.reactivationStrengthRun{ep(1)} = react_time_strength_run;
        RtimeStrength{ep(1)}.epochs = ep;
        RtimeStrength{ep(1)}.weights = ica_neuron_weights;
        RtimeStrength{ep(1)}.cellidx = cellidx;
        RtimeStrength{ep(1)}.runchange = rChange;

    end
    if savedata == 1
        save(sprintf('%s%s%s_icareactivationtimesReinstatementSWSSpkZero%02d.mat', dir,animalprefix,area,day), 'icareactivationtimes');
        save(sprintf('%s%s%s_RTimeStrengthSleepReinstatementSpkZero_%02d_%02d.mat', dir,animalprefix,area,bin,day), 'RtimeStrength');
    end
end
