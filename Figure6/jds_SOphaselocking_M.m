function jds_SOphaselocking_M(animalprefixlist, day, epochs)

%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates the preferred slow oscillation phase for ripples and spindles
%%------------------------------------------------------------------------
ctxPh = [];
hpPh = [];
spinPh = [];
ctxPh_sep = [];
hpPh_sep = [];
spinPh_sep = [];
alpha = linspace(-pi, pi, 50)';

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    load(sprintf('%s%sctxspindletime_SWS0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sslowoscdeltatimesSep_SWS%02d.mat', dir,animalprefix,day))
    load(sprintf('%s%srippletime_leadlag0%d.mat',dir,animalprefix,day));% get ripple time
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
    for ep=1:length(epochs)
        epoch = epochs(ep);
        if ~isempty(ripplecoupling{day}{epoch}.ctxleadriptimes)
            ctx_rippletimes = ripplecoupling{day}{epoch}.ctxleadriptimes(:,1);
            ca1_rippletimes = ripplecoupling{day}{epoch}.hplagriptimes(:,1);
        else
            ctx_rippletimes = [];
            ca1_rippletimes = [];
        end
        spindletimes = ctxspindle{day}{epoch}.starttime;
        SOtimes = [slowosc{day}{epoch}.SOstarttime-0.1 slowosc{day}{epoch}.SOendtime+0.1];
        t = slowosc{day}{epoch}.tvec;
        phdata = slowosc{day}{epoch}.phasedata';
        goodctx = isExcluded(ctx_rippletimes, SOtimes);
        goodhp = isExcluded(ca1_rippletimes, SOtimes);
        goodspin = isExcluded(spindletimes, SOtimes);

        if length(ctx_rippletimes)>1
            sph = phdata(lookup(ctx_rippletimes, t));
            sph = double(sph(logical(goodctx))) / 10000;
            ctxPh = [ctxPh; sph];
            [m, ph_ctx] = modulation(sph);
            ctxPh_sep = [ctxPh_sep; ph_ctx];
        end
        if length(ca1_rippletimes)>1
            sph = phdata(lookup(ca1_rippletimes, t));
            sph = double(sph(logical(goodhp))) / 10000;
            hpPh = [hpPh; sph];
            [m, ph_hp] = modulation(sph);
            hpPh_sep = [hpPh_sep; ph_hp];
        end
        if length(spindletimes)>1
            sph = phdata(lookup(spindletimes, t));
            sph = double(sph(logical(goodspin))) / 10000;
            spinPh = [spinPh; sph];
            [m, ph_spin] = modulation(sph);
            spinPh_sep = [spinPh_sep; ph_spin];
        end
    end
end

[p1,U2_1] = watsons_U2_approx_p(ctxPh_sep, hpPh_sep)
[p2,U2_2] = watsons_U2_approx_p(ctxPh_sep, spinPh_sep)
[p3,U2_3] = watsons_U2_approx_p(hpPh_sep, spinPh_sep)

[thetahat_ctx, kappa_ctx] = circ_vmpar(ctxPh);
[thetahat_hp, kappa_hp] = circ_vmpar(hpPh);
[thetahat_spin, kappa_spin] = circ_vmpar(spinPh);
[pdf_ctx] = circ_vmpdf(alpha,thetahat_ctx,kappa_ctx);
[pdf_hp] = circ_vmpdf(alpha,thetahat_hp,kappa_hp);
[pdf_spin] = circ_vmpdf(alpha,thetahat_spin,kappa_spin);

figure; hold on
pdf_ctx = [pdf_ctx; pdf_ctx(2:end)];
pdf_hp = [pdf_hp; pdf_hp(2:end)];
pdf_spin = [pdf_spin; pdf_spin(2:end)];

plot(normalize(pdf_ctx,'range'))
plot(normalize(pdf_hp,'range'))
plot(normalize(pdf_spin,'range'))
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')

figure; hold on
plot(pdf_ctx)
plot(pdf_hp)
plot(pdf_spin)
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
legend({'ctxrip','ca1rip','spin'})

keyboard