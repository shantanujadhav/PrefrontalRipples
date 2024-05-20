function jds_CA1_PFCrips_xcorr_M(animalprefixlist, epochs)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates cross-correlation between CA1 and PFC ripples. Can be modified
%to plot corss-corr between other events as well. Also calculates the 95%
%CIs based on shuffled data
%%------------------------------------------------------------------------
day = 1;
bin = 0.1;
tmax = 5;
sw1 = bin*3;
shufnum = 1000;
shuf = 0;

ripsxcorr = [];
ripsxcorr_shuf = [];
for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));% get ripple time
    load(sprintf('%s%srippletime_%sSWS0%d.mat',dir,animalprefix,rtype,day));
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));

    for ep=1:length(epochs)

        epoch = epochs(ep);

        ctx = ctxripple{day}{epoch};
        hp = ripple{day}{epoch};

        if ~isempty(sws{day}{epoch}.starttime)
            swsdur = sws{day}{epoch}.total_duration;

            if (~isempty(ctx)) && (~isempty(hp))

                ctxmidtimes = ctx.starttime;
                hpmidtimes = hp.starttime;

                ctxmidtimes_shuf = ctx.starttime;
                hpmidtimes_shuf = hp.starttime;
                iri = diff(ctxmidtimes_shuf); %use distribution of IRIs to jitter ctxriptimes
                iri2 = diff(hpmidtimes_shuf);
                iri(find(iri > 5)) = [];
                iri2(find(iri2 > 5)) = [];
                if shuf
                    for s = 1:shufnum
                        %%
                        %jitter ctx
                        sizeR = [1 length(ctxmidtimes_shuf)] ;
                        num1 = floor(length(ctxmidtimes_shuf)/2);
                        R = zeros(sizeR);  % set all to zero
                        ix = randperm(numel(R)); % randomize the linear indices
                        ix = ix(1:num1); % select the first
                        R(ix) = 1; % set the corresponding positions to 1
                        R(find(R == 0)) = -1;
                        if ~isempty(iri)
                            jittimes = (datasample(iri,length(ctxmidtimes_shuf)).*R')+(randi(5,1,length(R))'.*R');
                        else
                            jittimes = (randi(5,1,length(R))'.*R');
                        end

                        ctxmidtimes_shuf = sort(ctxmidtimes + jittimes);
                        %%
                        sizeR = [1 length(hpmidtimes_shuf)] ;
                        num1 = floor(length(hpmidtimes_shuf)/2);
                        R = zeros(sizeR);  % set all to zero
                        ix = randperm(numel(R)); % randomize the linear indices
                        ix = ix(1:num1); % select the first
                        R(ix) = 1; % set the corresponding positions to 1
                        R(find(R == 0)) = -1;
                        if ~isempty(iri2)
                            jittimes = (datasample(iri2,length(hpmidtimes_shuf)).*R')+(randi(5,1,length(R))'.*R');
                        else
                            jittimes = (randi(5,1,length(R))'.*R');
                        end
                        hpmidtimes_shuf = sort(hpmidtimes + jittimes);
                        %%
                        xc_shuf = spikexcorr(hpmidtimes, ctxmidtimes_shuf, bin, tmax);

                        if ~isempty(xc_shuf.c1vsc2)
                            Zcrosscov_shuf = zscore(xc_shuf.c1vsc2);

                            nstd=round(sw1/(xc_shuf.time(2) - xc_shuf.time(1))); % will be 3 std
                            g1 = gaussian(nstd, nstd);
                            timebase = xc_shuf.time;
                            bins_run = find(abs(timebase) <= tmax); % +/- Corrln window

                            Zcrosscov_sm_shuf = smoothvect(Zcrosscov_shuf, g1);% smoothed

                            ripsxcorr_shuf = [ripsxcorr_shuf; Zcrosscov_sm_shuf];
                        end
                    end
                end

                xc = spikexcorr(hpmidtimes, ctxmidtimes, bin, tmax);

                if ~isempty(xc.c1vsc2)
                    Zcrosscov = zscore(xc.c1vsc2);
                    prob = xc.c1vsc2./sum(xc.c1vsc2);
                    nstd=round(sw1/(xc.time(2) - xc.time(1))); % will be 3 std
                    g1 = gaussian(nstd, nstd);
                    timebase = xc.time;
                    bins_run = find(abs(timebase) <= tmax); % +/- Corrln window

                    Zcrosscov_sm = smoothvect(Zcrosscov, g1);% smoothed
                    prob_sm = smoothvect(prob, g1);
                    ripsxcorr = [ripsxcorr; Zcrosscov_sm];
                end
            end
        end
    end
end

figure
imagesc(ripsxcorr(randperm(length(ripsxcorr(:,1))),:))
colormap(inferno)

figure(1)
boundedline([-(tmax/bin):(tmax/bin)-1],nanmean(ripsxcorr),nanstd(ripsxcorr)./sqrt(size(ripsxcorr,1)),'-r');
xlabel(['Lag (bins) - ' num2str(bin) 'ms bins'])
ylabel('Z-score')
xlim([-40 40])
set(gcf, 'renderer', 'painters')

upperBnd = prctile(ripsxcorr_shuf,95);
lowerBnd = prctile(ripsxcorr_shuf,5);

%Compare to shuffled times
if shufnum > 0
    figure(1); hold on;
    plot([-(tmax/bin):(tmax/bin)-1],upperBnd,'--r')
    plot([-(tmax/bin):(tmax/bin)-1],lowerBnd,'--r')
end
keyboard
