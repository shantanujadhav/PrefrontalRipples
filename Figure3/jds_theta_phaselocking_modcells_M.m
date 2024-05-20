function jds_theta_phaselocking_modcells_M(animalprefixlist)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates and plot theta phase locking during run and REM for CA1 modulated cells
%%------------------------------------------------------------------------
day = 1;

daystring = '01';

%%
epochs = [2:2:16];

phaseExc = [];
phaseInh = [];

phaseExc_run = [];
phaseInh_run = [];
sphExc_run = [];
sphInh_run = [];
phaseExc_rem = [];
phaseInh_rem = [];
phaseExc_runAll = [];
phaseInh_runAll = [];
phaseExc_remAll = [];
phaseInh_remAll = [];
pdfExc_run = [];
pdfInh_run = [];
pdfExc_rem = [];
pdfInh_rem = [];
pdfExc_runAll = [];
pdfInh_runAll = [];
pdfExc_remAll = [];
pdfInh_remAll = [];

pre = 1;
post = 1;

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    %%
    %-----create the event matrix during SWRs-----%
    load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day)); % get spikes
    load(sprintf('%s%sCA1INTctxripallmodNewWin_epsIncludeHigh%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sthetatime%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%srem%02d.mat',dir,animalprefix,day));
    rem = rem;
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    modcells = epochModulation.cellidx;
    
    tets = tetinfo{1}{epochs(1)}; 
    
    reftet = []; 
    for t = 1:length(tets)
        tmp = tets{t};
        if isfield(tmp, 'descrip')
            if isequal(tmp.descrip, 'CA1Ref')
                reftet = [reftet; t];
            end
        end
    end
    
    sleeps = [1 1; 2 3; 3 5; 4 7; 5 9; 6 11; 7 13; 8 15; 9 17];
    
    allcellphase_ep = [];
    sigmod_ep = [];
    for ep = 1:length(epochs)

        if pre == 1
            eps = [epochs(ep) (epochs(ep) + 1)];
        elseif post == 1
            eps = [epochs(ep) (epochs(ep) - 1)];
        end
        
        ep2 = find(sleeps(:,2) == eps(2));
        
        inhcells = modcells(find(epochModulation.modMat(:,ep2) == -1),:);
        inhcells(:,3) = -1;
        exccells = modcells(find(epochModulation.modMat(:,ep2) == 1),:);
        exccells(:,3) = 1;
        
        allmodcells = [inhcells; exccells];
        allcellphase = [];
        sigmod = [];
        
        cellnum = length(allmodcells(:,1));

        celldata = [];
        spikecounts = [];

        eprun = eps(1);
        epsleep = eps(2);
        
        if (eprun <10) && (isequal(animalprefix, 'ZT2'))
            epochstring = ['0',num2str(eprun)];
        else
            epochstring = num2str(eprun);
        end
        
        if (eps(2) <10) && (isequal(animalprefix, 'ZT2'))
            epochstringrem = ['0',num2str(eps(2))];
        else
            epochstringrem = num2str(eps(2));
        end
        
        thetalist = thetatime{day}{eprun};
        thetalist = [thetalist.starttime thetalist.endtime];
        remlist = rem{day}{eps(2)};
        remlist = [remlist.starttime remlist.endtime];
        
        for cellcount = 1:cellnum %get spikes for each cell
            cellmod = allmodcells(cellcount, 3);
            index = [day,eprun,allmodcells(cellcount,[1 2])];
            index2 = [day,epsleep,allmodcells(cellcount,[1 2])];
            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            else
                spiketimes = [];
            end
            if ~isempty(spikes{index2(1)}{index2(2)}{index2(3)}{index2(4)})
                spiketimes2 = spikes{index2(1)}{index2(2)}{index2(3)}{index2(4)}.data(:,1);
            else
                spiketimes2 = [];
            end
            
            goodspikes = isExcluded(spiketimes, thetalist);
            goodspikes2 = isExcluded(spiketimes2, remlist);
            
            tet = allmodcells(cellcount,1);
            
            if (reftet<10)
                reftetstring = ['0',num2str(reftet)];
            else
                reftetstring = num2str(reftet);
            end
            
            curreegfile = [dir,'/EEG/',animalprefix,'thetagnd', daystring,'-',epochstring,'-',reftetstring];
            load(curreegfile);
            thetagnd_run = thetagnd; clear thetagnd;
            
            curreegfileRem = [dir,'/EEG/',animalprefix,'thetagnd', daystring,'-',epochstringrem,'-',reftetstring];
            load(curreegfileRem);
            thetagnd_rem = thetagnd; clear thetagnd;
            
            phasedata = thetagnd_run{day}{eprun}{reftet}.data(:,2);
            
            phasedata2 = thetagnd_rem{day}{epsleep}{reftet}.data(:,2);
            
            t = geteegtimes(thetagnd_run{day}{eprun}{reftet});
            t2 = geteegtimes(thetagnd_rem{day}{epsleep}{reftet});
            
            if length(goodspikes)~=0
                sph = phasedata(lookup(spiketimes, t));
                sph = double(sph(logical(goodspikes))) / 10000;  % If no spikes, this will be empty
            else
                sph = [];
            end
            
            if length(goodspikes2)~=0
                sph2 = phasedata2(lookup(spiketimes2, t2));
                sph2 = double(sph2(logical(goodspikes2))) / 10000;  % If no spikes, this will be empty
            else
                sph2 = [];
            end
            
            if (length(sph)>1) && (length(sph2)>1)
                
                % Rayleigh and Modulation: Originally in lorenlab Functions folder
                stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                [m, ph] = modulation(sph);
                phdeg = ph*(180/pi);
                % Von Mises Distribution - From Circular Stats toolbox
                [thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                thetahat_deg = thetahat*(180/pi);
                
                
                [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                
                stats2 = rayleigh_test(sph2); % stats.p and stats.Z, and stats.n
                [m2, ph2] = modulation(sph2);
                phdeg2 = ph2*(180/pi);
                % Von Mises Distribution - From Circular Stats toolbox
                [thetahat2, kappa2] = circ_vmpar(sph2); % Better to give raw data. Can also give binned data.
                thetahat_deg2 = thetahat2*(180/pi);
                
                [prayl2, zrayl2] = circ_rtest(sph2); 
                
                % Make finer polar plot and overlay Von Mises Distribution Fit.
                % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                % -------------------------------------------------------------
                nbins = 50;
                bins = -pi:(2*pi/nbins):pi;
                count = histc(sph, bins);
                count2 = histc(sph2, bins);
                
                % Make Von Mises Fit
                alpha = linspace(-pi, pi, 50)';
%                 alpha = linspace(-pi, pi, 100)';
                [pdf] = circ_vmpdf(alpha,thetahat,kappa);
                [pdf2] = circ_vmpdf(alpha,thetahat2,kappa2);

                if (prayl < 0.05) && (prayl2 < 0.05) %get proportion shifters here
                    if cellmod == 1
                        run_rem_ph = [ph ph2];
                        phaseExc = [phaseExc; run_rem_ph];
                    elseif cellmod == -1
                        run_rem_ph = [ph ph2];
                        phaseInh = [phaseInh; run_rem_ph];
                    end 
                end
                if prayl < 0.05 
                    if cellmod == 1
                        phaseExc_run = [phaseExc_run; ph];
                        pdfExc_run = [pdfExc_run; pdf'];
                        sphExc_run = [sphExc_run; mean(sph)];
                    elseif cellmod == -1
                        phaseInh_run = [phaseInh_run; ph];
                        pdfInh_run = [pdfInh_run; pdf'];
                        sphInh_run = [sphInh_run; mean(sph)];
                    end 
                end
                if prayl2 < 0.05
                    if cellmod == 1
                        phaseExc_rem = [phaseExc_rem; ph2];
                        pdfExc_rem = [pdfExc_rem; pdf2'];
                    elseif cellmod == -1
                        phaseInh_rem = [phaseInh_rem; ph2];
                        pdfInh_rem = [pdfInh_rem; pdf2'];
                    end 
                end
                
                if cellmod == 1
                    phaseExc_runAll = [phaseExc_runAll; ph];
                    phaseExc_remAll = [phaseExc_remAll; ph2];
                    pdfExc_runAll = [pdfExc_runAll; pdf'];
                    pdfExc_remAll = [pdfExc_remAll; pdf2'];
                elseif cellmod == -1
                    phaseInh_runAll = [phaseInh_runAll; ph];
                    phaseInh_remAll = [phaseInh_remAll; ph2];
                    pdfInh_runAll = [pdfInh_runAll; pdf'];
                    pdfInh_remAll = [pdfInh_remAll; pdf2'];
                end
            end
        end
    end
end
[p1,U2_1] = watsons_U2_approx_p(phaseInh_run, phaseInh_rem)
[p2,U2_2] = watsons_U2_approx_p(phaseExc_run, phaseExc_rem)
[p3,U2_3] = watsons_U2_approx_p(phaseInh_rem, phaseExc_rem)

shiftExc = phaseExc*(180/pi);
shiftInh = phaseInh*(180/pi);

shiftAll = [shiftExc; shiftInh];
changeAll = abs(shiftAll(:,2)-shiftAll(:,1));
nonShiftCells = shiftAll(find(changeAll<90),2);
shiftCells = shiftAll(find(changeAll>90),2);
nonShiftCells = nonShiftCells*(pi/180);
shiftCells = shiftCells*(pi/180);

figure; histogram(phaseExc_run*(180/pi),25);
hold on
histogram(phaseInh_run*(180/pi),25);

figure; histogram(phaseExc_run,25);
hold on
histogram(phaseInh_run,25);

pdfexc = [pdfExc_run pdfExc_run(:,(2:end))];
pdfinh = [pdfInh_run pdfInh_run(:,(2:end))];

pdfexc_rem = [pdfExc_rem pdfExc_rem(:,(2:end))];
pdfinh_rem = [pdfInh_rem pdfInh_rem(:,(2:end))];

figure;
shadedErrorBar([1:99],mean(pdfexc),nanstd(pdfexc)./sqrt(size(pdfexc,1)),'-b',0); hold on
shadedErrorBar([1:99],mean(pdfexc_rem),nanstd(pdfexc_rem)./sqrt(size(pdfexc_rem,1)),'-r',0);
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
title(['CA1mod Theta Phase Locking PDF - EXC - RUN vs REM p=' num2str(p2) 'U=' num2str(U2_2)])
set(gcf, 'renderer', 'painters')
figure;
shadedErrorBar([1:99],mean(pdfinh),nanstd(pdfinh)./sqrt(size(pdfinh,1)),'-b',0); hold on
shadedErrorBar([1:99],mean(pdfinh_rem),nanstd(pdfinh_rem)./sqrt(size(pdfinh_rem,1)),'-r',0);
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
title(['CA1mod Theta Phase Locking PDF - INH - RUN vs REM p=' num2str(p1) 'U=' num2str(U2_1)])
set(gcf, 'renderer', 'painters')
figure;
shadedErrorBar([1:99],mean(pdfexc),nanstd(pdfexc)./sqrt(size(pdfexc,1)),'-b',0); hold on
shadedErrorBar([1:99],mean(pdfinh),nanstd(pdfinh)./sqrt(size(pdfinh,1)),'-r',0);
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
title('CA1mod Theta Phase Locking PDF')

figure;
shadedErrorBar([1:99],mean(pdfexc_rem),nanstd(pdfexc_rem)./sqrt(size(pdfexc_rem,1)),'-b',0); hold on
shadedErrorBar([1:99],mean(pdfinh_rem),nanstd(pdfinh_rem)./sqrt(size(pdfinh_rem,1)),'-r',0);
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
title('CA1mod Theta Phase Locking PDF')

%% COMPARE RUN AND SLEEP DISTRIBUTIONS
figure
subplot(2,1,1)
plot(mean(pdfinh_rem),'LineWidth',2);
hold on
plot(mean(pdfinh),'LineWidth',2);
title('CA1inh Theta Modulation - Run/REM')
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
legend({'REM','RUN'})
subplot(2,1,2)
plot(mean(pdfexc_rem),'LineWidth',2);
hold on
plot(mean(pdfexc),'LineWidth',2);
title('CA1exc Theta Modulation - Run/REM')
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')

%%
figure
hold on
plot(mean(pdfinh),'b','LineWidth',2);
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
plot(mean(pdfexc),'r','LineWidth',2);
title('CA1Int Theta Modulation - Run/REM')
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
legend({'INH','EXC'})

keyboard

