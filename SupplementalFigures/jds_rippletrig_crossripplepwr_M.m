function jds_rippletrig_crossripplepwr_M(animalprefixlist,rtype)
%%------------------------------------------------------------------------
%Justin D. Shin

%Plots ripple triggered ripple power in the opposing area to show increase
%and flat (or decrease) in ripple power aligned to coordinated and
%independent events, respectively.
%%------------------------------------------------------------------------

allenv1 = [];
allenv2 = [];

g1 = gaussian(50,50);
epochs = 1:2:17;
days = 1;

for r = 1:2
    for a = 1:length(animalprefixlist)

        animalprefix = animalprefixlist{a};
        
        for d = 1:length(days)
            day = days(d);
            dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

            daystring = ['0',num2str(day)];
            
            %Align to these ripples
            load(sprintf('%s%s%sctxrippletime_coordSWS0%d.mat',dir,animalprefix,rtype,day));% get ripple time

            load(sprintf('%s%stetinfo.mat',dir,animalprefix));

            tets = tetinfo{1}{1};

            tetrodes = [];
            for t = 1:length(tets)
                tmp = tets{t};
                if isfield(tmp, 'descrip')
                    if strcmp(rtype,'PFC')
                        if r == 1
                            if isequal(tmp.descrip, 'ctxriptet')
                                tetrodes = [tetrodes; t];
                            end
                        elseif r == 2
                            if isequal(tmp.descrip, 'riptet')
                                tetrodes = [tetrodes; t];
                            end
                        end
                    else
                        if r == 1
                            if isequal(tmp.descrip, 'riptet')
                                tetrodes = [tetrodes; t];
                            end
                        elseif r == 2
                            if isequal(tmp.descrip, 'ctxriptet')
                                tetrodes = [tetrodes; t];
                            end
                        end
                    end
                end
            end

            for ep=1:length(epochs)
                epoch = epochs(ep);

                if epoch <10
                    epochstring = ['0',num2str(epoch)];
                    epochstring2 = ['0',num2str(epoch)];
                else
                    epochstring = num2str(epoch);
                    epochstring2 = num2str(epoch);
                end

                if strcmp('ER1', animalprefix)
                    epochstring = num2str(epoch);
                end

                rips = ctxripple{day}{epoch};
                
                if strcmp('ER1', animalprefix)
                    curreegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring2,'-','01']; %use tet1 to extract starttime of epoch
                else
                    curreegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring,'-','01']; %use tet1 to extract starttime of epoch
                end
                load(curreegfile);

                times = geteegtimes(eegref{day}{epoch}{1});

                riplist = [rips.starttime rips.endtime];

                if length(riplist(:,1)) > 20
                    allenv = [];
                    for ii = 1:length(tetrodes)
                        allenvtmp = [];
                        tetrode = tetrodes(ii);
                        if tetrode <10
                            tetstring = ['0',num2str(tetrode)];
                        else
                            tetstring = num2str(tetrode);
                        end

                        currripfile = [dir,'/EEG/',animalprefix,'ripple', daystring,'-',epochstring,'-',tetstring];
                        ripeeg = load(currripfile);
                        env = zscore(double(ripeeg.ripple{day}{epoch}{tetrode}.data(:,3)));
                        for i = 1:length(riplist(:,1))
                            riptimetmp = [riplist(i,1) riplist(i,2)];
                            tmp_idx = find(times >= riptimetmp(1) & times <= riptimetmp(2));
                            tmp_idx2 = lookup(riptimetmp(1), times);
                            envtmp = max(env(tmp_idx));
                            if ((tmp_idx2-1000) > 1) && ((tmp_idx2+1000) < length(times))
                                envtmp2 = smoothvect(env((tmp_idx2-1000):(tmp_idx2+1000)),g1);
                                allenvtmp = [allenvtmp; envtmp2'];
                            end
                        end
                        allenv = [allenv; zscore(mean(allenvtmp,1))];
                    end
                    if r == 1
                        allenv1 = [allenv1; allenv];
                    elseif r == 2
                        allenv2 = [allenv2; allenv];
                    end
                end
                clear riplist
            end
        end
    end
end
figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-1000:1000],mean(allenv1),'-k','LineWidth',1)
boundedline([-1000:1000],mean(allenv1),(std(allenv1)./sqrt(length(allenv1(:,1)))),'-k');
pl2 = plot([-1000:1000],mean(allenv2),'-r','LineWidth',1)
boundedline([-1000:1000],mean(allenv2),(std(allenv2)./sqrt(length(allenv2(:,1)))),'-r');
legend([pl1 pl2],{'Intra','Other'})
xlim([-500 500])

keyboard

