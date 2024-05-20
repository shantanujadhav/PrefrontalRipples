%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates cross-region reactivation strength in sleep based on run
%activity in CA1 and PFC using PCA/ICA. Uses run session preceding sleep to calculate
%reactivation strength. Cells tracked over 2 consecutive run/sleep epochs
%%------------------------------------------------------------------------
clc;
close all;
clear all;
%%
savedata = 1;
numCross = 0;
animalprefixlist = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8','ER1'};
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    day = 1;
    bin = 50; % ms
    bintime = bin./1000./60; %min

    sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
    %%
    load(sprintf('%s%s_spikematrix_ev_allepochallcell50_0%d',dir,animalprefix,day));
    load(sprintf('%s%ssws0%d.mat',dir,animalprefix,day));
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
    idx = find(swsdurs>=0);

    %compile epochs to analyze
    eps = [];
    for t = 1:length(idx) %need to modify in case of gaps
        if idx(t) > 1 %dont use sleep 1
            sleep = sleeps(idx(t),2);
            run = sleep - 1;
            eps = [eps; run sleep];
        end
    end
    icareactivationtimes = [];
    for e = 1:length(eps(:,1))
        ep = eps(e,:);

        %%
        total_ctx_neuron = observation_matrix{ep(1)}.ncortical;
        total_hp_neuron =  observation_matrix{ep(1)}.nhp;
        cellidx = [observation_matrix{ep(2)}.ctxidx;observation_matrix{ep(2)}.hpidx];
        %---------Sleep spike matrix-------%
        POST_spike = [observation_matrix{ep(2)}.ctxdata;observation_matrix{ep(2)}.hpdata];
        %---------W1 and W2 spike matrix-------%
        W2_spike = [observation_matrix{ep(1)}.ctxdata;observation_matrix{ep(1)}.hpdata];
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
        rippost = ripple{day}{ep(2)}; %Use PFC ripple times for the sleep sessions
        riplist_POST(:,1) = rippost.starttime;
        riplist_POST(:,2) = rippost.endtime;
        times_filteeg = observation_matrix{ep(2)}.timeeeg;
        timevect = times_filteeg(1)*1000:bin:times_filteeg(end)*1000;
        POST_rip=zeros(size(timevect));
        for rip_seg =1:length(riplist_POST(:,2))
            indtemp = find(timevect >= riplist_POST(rip_seg,1).*1000 & timevect < riplist_POST(rip_seg,2).*1000);
            POST_rip(indtemp) = 1;
        end

        ripw2 = ripplerun{day}{ep(1)};
        riplist_W2(:,1) = ripw2.starttime;
        riplist_W2(:,2) = ripw2.endtime;
        times_filteeg = observation_matrix{ep(1)}.timeeeg;
        timevect = times_filteeg(1)*1000:bin:times_filteeg(end)*1000;
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
        idx_nrip_w2 = find(W2_rip==0);
        W2_spike_raw = W2_spike;
        qW2_spike = qW2_spike(:,idx_nrip_w2); %get rid of ripple times

        qW2_spike = zscore(qW2_spike')'; %Zscore for PCA
        qW2_spike(isnan(qW2_spike)) = 0;
     
        %%
        %-----------correlation matrix----------%
        cW2_spike =corr(qW2_spike'); %correlation matrix from run session excluding ripples
        %%
        %---------------PCA-----------------%
        [u0,s0,v0] = svd(cW2_spike);
        sdiag0 = diag(s0);

        %-----PCs significant test-----%
        %Is this how the threshold should be defined for cross regional
        %correlation?
        lamdamax = (1+sqrt((total_ctx_neuron + total_hp_neuron)/length(qW2_spike(1,:)))).^2;%+(length(qW2_spike(:,1))).^(-2/3);
        diagind = find(sdiag0 > lamdamax);

        numPCs = length(diagind);
        if numPCs > 0
            %ICA
            Psign = u0(:,[1:numPCs]); %Restrict ICA to significant eigenvalues
            Zproj = Psign'*qW2_spike; %project original zscored spk matrix into PC space spanning sig PCs
            
            %%
            %FastICA
            %[icasig, A, W] = fastica (Zproj); %Use un-mixing matrix W to get V. V = Psign*W

            %%
            %RobustICA - used robustICA due to insensitivity to
            %initialization. robust output
            [S, H, iter, W] = robustica(Zproj,{});
            %%
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
                    w_vec = w_vec*(-1); %flip signs if max is on (-)weight side
                end
                vtmp = [vtmp w_vec];
            end

            V = vtmp;

            %Use this section to find high weight cells
            ica_neuron_weights = [];
            cW2_spike_ic = cell(1,numPCs);
            crossMembers = [];
            for i = 1:numPCs
                wtsTmp = V(:,diagind(i));
                meanWts = mean(wtsTmp);
                wtsThresh = meanWts + std(wtsTmp)*2;
                highIdx = find(wtsTmp > wtsThresh);
                tmpmat = V(:,diagind(i))*V(:,diagind(i))';
                highCtx = find(highIdx <= total_ctx_neuron);
                highHp = find(highIdx >= total_ctx_neuron+1);
                if (~isempty(highCtx)) && (~isempty(highHp))
                    crossMembers = [crossMembers; 1]; %if member cells in CA1 AND PFC
                else
                    crossMembers = [crossMembers; 0];
                end
                
                tmpmat2 = tmpmat - diag(diag(tmpmat)); %set diag to 0
                tmpmat2(1:total_ctx_neuron,1:total_ctx_neuron) = 0;
                tmpmat2(total_ctx_neuron+1:end,total_ctx_neuron+1:end) = 0; %set intraregion coactivation to 0 so only cross-regional events are calculated
                cW2_spike_ic{i} = tmpmat2;
                ica_neuron_weights{i} = [[observation_matrix{ep(2)}.ctxidx V(1:total_ctx_neuron,i) zeros(total_ctx_neuron,1)]; ...
                [observation_matrix{ep(2)}.hpidx V(total_ctx_neuron+1:end,i) ones(total_hp_neuron,1)]];
            end
            numCross = numCross + sum(crossMembers);

            thetaratio = sdiag0(diagind)./lamdamax;
            %%
            %------zscore--------%
            qPOST_spike_z = zscore(POST_spike')';

            qW2_spike_all_z = zscore(W2_spike_raw')'; %all spike matrix during run for projection

            qW2_spike_all_z(isnan(qW2_spike_all_z)) = 0; %Zscored spike matrix for awake ALL
            qPOST_spike_z(isnan(qPOST_spike_z)) = 0; %Zscored spike matrix for POST SLEEP ALL

            %%
            R_W2 = zeros(1,length(qW2_spike_all_z(1,:)));
            R_POST = zeros(1,length(qPOST_spike_z(1,:)));

            for pc = 1:length(diagind)
                PC = cW2_spike_ic{pc};
                for t = 1:length(qW2_spike_all_z(1,:))
                    R1_W_temp = qW2_spike_all_z(:,t)'*PC*qW2_spike_all_z(:,t);
                    R_W2(pc,t) =  R1_W_temp;
                end
                for t = 1:length(qPOST_spike_z(1,:))
                    R1_POST_temp = qPOST_spike_z(:,t)'*PC*qPOST_spike_z(:,t);
                    R_POST(pc,t) =  R1_POST_temp;
                end
            end

            %%
            %%
            reactivation_strength_run = R_W2';
            reactivation_strength_post = R_POST';
           
            clear riplist_POST riplist_W2 nriplist_W2

            times_filteeg_post = observation_matrix{ep(2)}.timeeeg;
            posttimevect = times_filteeg_post(1)*1000:bin:times_filteeg_post(end)*1000;
            times_filteeg_w = observation_matrix{ep(1)}.timeeeg;
            wtimevect = times_filteeg_w(1)*1000:bin:times_filteeg_w(end)*1000;
            %Find the mean and STD of reactivation strength during SWS, Extract
            %times where the signal is +3SD above mean, and record times.

            react_time_strength_run = [];
            react_time_strength_sleep = [];

            for ii = 1:length(reactivation_strength_run(1,:))

                wpctmp = reactivation_strength_run(:,ii);
                postpctmp = reactivation_strength_post(:,ii);

                react_time_strength_run{ii} = [(wtimevect./1000)' wpctmp];

                react_time_strength_sleep{ii} = [(posttimevect./1000)' postpctmp];

                icareactivationtimes{e}.post_strengths{ii} = [(posttimevect./1000)' postpctmp]; %All event strengths, no filter
                icareactivationtimes{e}.run_strengths{ii} = [(wtimevect./1000)' wpctmp];
                icareactivationtimes{e}.thetaratio = thetaratio;
                icareactivationtimes{e}.pc_weights = ica_neuron_weights;
                icareactivationtimes{e}.crossmembers = crossMembers;
                icareactivationtimes{e}.cellidx = cellidx;
                icareactivationtimes{e}.timebinsize = bin;
                icareactivationtimes{e}.descrip = '0 is PFC; 1 is CA1';
            end
            icareactivationtimes{e}.epochs = ep;

            RtimeStrength{ep(2)}.reactivationStrengthRun = react_time_strength_run;
            RtimeStrength{ep(2)}.reactivationStrength = react_time_strength_sleep;
            RtimeStrength{ep(2)}.weights = ica_neuron_weights;
            RtimeStrength{ep(2)}.crossmembers = crossMembers;
            RtimeStrength{ep(2)}.epochs = ep;
            RtimeStrength{ep(2)}.cellidx = cellidx;
            RtimeStrength{ep(2)}.descrip = '0 is PFC; 1 is CA1';
            clear swslist_POST swslist riplist_POST riplist_W2 nriplist_W2
        end
    end
    if savedata == 1
        save(sprintf('%s%sCA1PFC_icareactivationtimes%02dSWSSpk%02d.mat', dir,animalprefix,bin,day), 'icareactivationtimes');
        save(sprintf('%s%sCA1PFC_RTimeStrengthSleepNewSpk_%02d_%02d.mat', dir,animalprefix,bin,day), 'RtimeStrength');
    end
end

keyboard



