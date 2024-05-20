function jds_replay_spikematrix_ev_singleday_M(animalprefixlist,day,ep)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates spatial information for CA1 mod cells
%%------------------------------------------------------------------------

bin = 20; %set time bin size in ms (20 ms for reactivation analysis; 5 ms for MUA analysis)
%%
%----save the results?-----%
savedata = 1;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    animaldir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    exclude_list = [];
    eegdir = [animaldir,'EEG/'];

    %%
    %---------------------spike matrix----------------------%
    for e = 1:length(ep)

        epoch = ep(e);
        
        %specify epochs to match cells over
        if (mod(epoch,2) == 0 || epoch == 1)
            eps = [epoch epoch+1];
        else
            eps = [epoch epoch-1];
        end
        
        %change cell matching function to include other cell types as well
        %(ex. CA1 interneurons)
        [ctxidx, hpidx] = matchidx_acrossep_singleday(animaldir, animalprefix, day, eps, exclude_list); %(tet, cell)
        ctxnum = size(ctxidx,1);
        hpnum = size(hpidx,1);

        tmpflist1 = sprintf('%s%seeg%02d-%02d-%02d.mat', eegdir,animalprefix, day, epoch, ctxidx(1,1));

        load(tmpflist1);
        times_filteeg = geteegtimes(eeg{day}{epoch}{ctxidx(1,1)}) ;
        times_filteeg = times_filteeg(:)';% reference time is eeg time
        timevect = times_filteeg(1)*1000:bin:times_filteeg(end)*1000;
        nbins = length(timevect);
        spikes = loaddatastruct(animaldir, animalprefix, 'spikes', day); % get spikes
        if ctxnum > 0
            for i = 1:ctxnum
                % get cortical neurons
                index = [day,epoch,ctxidx(i,:)];
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                spkstime = ( double(spiketimes)-times_filteeg(1) ).*1000; % ms

                if bin == 1
                    spkstime = round((double(spiketimes)-times_filteeg(1) ).*1000); % ms
                    for ntime = 1:nbins
                        nspks = length(find(spkstime == ntime*bin));
                        spike_matrix_ctx(i,ntime) = nspks;
                    end
                else
                    for ntime = 1:nbins
                        nspks = length(find(spkstime >= (ntime-1)*bin + 1 & spkstime < ntime*bin));
                        spike_matrix_ctx(i,ntime) = nspks;
                    end
                end
            end
        else
            spike_matrix_ctx = [];
        end

        if hpnum >0
            for i = 1:hpnum
                index = [day,epoch,hpidx(i,:)] ;
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                spkstime = ( double(spiketimes)-times_filteeg(1) ).*1000; % ms

                if bin == 1
                    spkstime = round((double(spiketimes)-times_filteeg(1) ).*1000); % ms
                    for ntime = 1:nbins
                        nspks = length(find(spkstime == ntime*bin));
                        spike_matrix_hp(i,ntime) = nspks;
                    end
                else
                    for ntime = 1:nbins
                        nspks = length(find(spkstime >= (ntime-1)*bin + 1 & spkstime < ntime*bin));
                        spike_matrix_hp(i,ntime) = nspks;
                    end
                end
            end
        else
            spike_matrix_hp = [];
        end

        observation_matrix{epoch}.ctxdata =  spike_matrix_ctx;
        observation_matrix{epoch}.hpdata =  spike_matrix_hp;
        observation_matrix{epoch}.hpidx =  hpidx;
        observation_matrix{epoch}.ctxidx =  ctxidx;
        observation_matrix{epoch}.ncortical =  ctxnum;
        observation_matrix{epoch}.nhp =  hpnum;
        observation_matrix{epoch}.timeeeg =  times_filteeg;
        clear spike_matrix_hp spike_matrix_ctx
    end
    %%
    if savedata
        save(sprintf('%s%s_spikematrix_ev_allepochallcell%d_%02d.mat', animaldir,animalprefix,bin,day), 'observation_matrix', '-v7.3');
    end
end