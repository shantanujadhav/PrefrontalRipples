function jds_relativeRippleRate_M(animalprefixlist, eps)
%%------------------------------------------------------------------------
%Justin D. Shin

%Calculates the rate of ripples in the second half of sleep relative to the
%first half. Also calculates the change in the difference between
%independent and coordinated ripple rate for first and second half of sleep
%(within an epoch)
%%------------------------------------------------------------------------
ncChange = [];
cChange = [];

diffFirst = [];
diffSecond = [];

day = 1; 
daystring = '01';
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    % get ripple times
    load(sprintf('%s%srippletime_noncoordSWS0%d.mat',dir,animalprefix,day));
    ncrip = ripple; clear ripple
    load(sprintf('%s%srippletime_coordSWS0%d.mat',dir,animalprefix,day));
    crip = ripple;
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get sws time

    for ep = eps
        if ep <10
            epochstring = ['0',num2str(ep)];
        else
            epochstring = num2str(ep);
        end
        swsDur = sws{day}{ep}.total_duration;
        halfsleep = swsDur/2; %amount to div # of ripples by
        halftstamps = halfsleep*1500; %1500 is sampling rate
        criptimes = crip{day}{ep}.starttime;
        ncriptimes = ncrip{day}{ep}.starttime;

        if halfsleep > 60 %at least 30 seconds for each half ripple rate calculation
            swslist = [sws{day}{ep}.starttime sws{day}{ep}.endtime];

            curreegfile = [dir,'/EEG/',animalprefix,'eeg', daystring,'-',epochstring,'-','01']; %use any tetrode
            load(curreegfile);
            time1 = geteegtimes(eeg{day}{ep}{1}); % construct time array

            [~,swsvec] = wb_list2vec(swslist,time1);
            eqn_c = cumsum(swsvec) == halftstamps;
            midtstamp = time1(find(eqn_c==1));
            ncratefirst = length(find(ncriptimes<midtstamp))/halfsleep;
            ncratesecond = length(find(ncriptimes>midtstamp))/halfsleep;
            ncrelative = ((ncratesecond-ncratefirst)/ncratefirst)*100;

            cratefirst = length(find(criptimes<midtstamp))/halfsleep;
            cratesecond = length(find(criptimes>midtstamp))/halfsleep;
            crelative = ((cratesecond-cratefirst)/cratefirst)*100;

            diffFirst = [diffFirst; cratefirst-ncratefirst];
            diffSecond = [diffSecond; cratesecond-ncratesecond];

            ncChange = [ncChange; ncrelative];
            cChange = [cChange; crelative];
        end
    end
end
figure
hold on
bar([mean(ncChange) mean(cChange)],'k')
errorbar([1:2],[mean(ncChange) mean(cChange)],...
    [std(ncChange)./sqrt(length(ncChange)) ...
    std(cChange)./sqrt(length(cChange))],'k','LineStyle','none')
keyboard
