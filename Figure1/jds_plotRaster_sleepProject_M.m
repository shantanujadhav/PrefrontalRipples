%%------------------------------------------------------------------------
%Justin D. Shin

%Plots a raster during specified epoch - Includes CA1 and PFC cells, CA1
%and PFC EEG with ripples highlighted
%%------------------------------------------------------------------------
clear; close all;
prefix='ZT2';
day=1; 
epochs=[11]; % epochs to plot

% -------------------------------
win = 10; % 10 second window
overlap = 5; % with 5 second overlap of previous

% Get Data Location
% -------------------------------
switch prefix
    case 'ZT2'
        directoryname = '/Volumes/JUSTIN/SingleDay/ZT2_direct/';
        dire = '/Volumes/JUSTIN/SingleDay/ZT2_direct/';
        animdirect = directoryname;
        riptetlist = [10 11 12 14 16 17 18 19 29 24 25 27 32 36];  % No Need if no ripples
        maineegtet = 14;  % CA1 tet % No need if no EEG
        peegtet = 40; % PFCtet % No need if no EEG

    case 'KL8'
        directoryname = '/Volumes/JUSTIN/SingleDay/KL8_direct/';
        dire = '/Volumes/JUSTIN/SingleDay/KL8_direct/';
        animdirect = directoryname;
        riptetlist = [10];  % No Need if no ripples
        maineegtet = 10;  % CA1 tet % No need if no EEG
        peegtet = 5; % PFCtet % No need if no EEG

    case 'ER1'
        directoryname = '/Volumes/JUSTIN/SingleDay/ER1_direct/';
        dire = '/Volumes/JUSTIN/SingleDay/ER1_direct/';
        animdirect = directoryname;
        riptetlist = [13];  % No Need if no ripples
        maineegtet = 13;  % CA1 tet % No need if no EEG
        peegtet = 27; % PFCtet % No need if no EEG
end

currdir = pwd;
if (directoryname(end) == '/')
    directoryname = directoryname(1:end-1);
end
if (dire(end) == '/')
    dire = dire(1:end-1);
end

if (day < 10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

% SET DATA
eegtets = riptetlist;

% Also get a PFC eeg
eegtets = [eegtets, peegtet];
peegidx= find(eegtets==peegtet);
maineegidx = find(eegtets==maineegtet);
saveg = 0;

for e = 1:length(epochs)
    epoch = epochs(e);
    % Get data files
    % Spike data
    %-----------
    spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
    load(spikefile);
    tetinfofile = sprintf('%s/%stetinfo.mat', directoryname, prefix);
    load(tetinfofile);
    cellinfofile = sprintf('%s/%scellinfo.mat', directoryname, prefix);
    load(cellinfofile);

    % Load ripple files
    %-------------------------------
    ripfile = sprintf('%s/%srippletime_noncoordSWS%02d.mat', directoryname, prefix, day);
    load(ripfile);
    nc_ripple = ripple;

    coordripfile = sprintf('%s/%sripplecoordinationSWS%02d.mat', directoryname, prefix,day);
    load(coordripfile);
    c_ripple = ripple;

    pripfile = sprintf('%s/%sctxrippletime_noncoordSWS%02d.mat', directoryname, prefix, day);
    load(pripfile);

    posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
    load(posfile);
    
    % Get cells
    % ---------
    filterString = 'strcmp($tag2, ''CA1Pyr'') && ($numspikes > 100)';
    
    cellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    cellsi = [repmat([day epoch], size(cellindices,1),1 ), cellindices];
    usecellsi = 1:size(cellsi,1);
   
    filterString = 'strcmp($area, ''PFC'') && ($numspikes > 100)';
    pcellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    cellsp = [repmat([day epoch], size(pcellindices,1),1 ), pcellindices]; 
    usecellsp = 1:size(cellsp,1);
   
    % GET Spike data
    %-----------
    
    % CA1 cells
    for i=1:size(cellsi,1)
        eval(['spiketimei{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,1);']);
        eval(['spikeposi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,2:3);']);
        eval(['spikeposidxi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,7);']);
    end
    
    % PFC cells
    for i=1:size(cellsp,1)
        eval(['spiketimep{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,1);']);
        eval(['spikeposp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,2:3);']);
        eval(['spikeposidxp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,7);']);
    end
    
    %% Independent ca1 ripples
    riptimes = [nc_ripple{day}{epoch}.starttime nc_ripple{day}{epoch}.endtime]; 
    rip_starttime = riptimes(:,1);
    triggers = rip_starttime;
    
    %% Coordinated ripples
    c_riptimes = [ripplecoordination{day}{epoch}.starttime ripplecoordination{day}{epoch}.endtime]; 
    c_rip_starttime = c_riptimes(:,1);
    c_triggers = c_rip_starttime;
    
    %% Independent PFC ripples
    p_riptimes = [ctxripple{day}{epoch}.starttime ctxripple{day}{epoch}.endtime]; 
    p_rip_starttime = p_riptimes(:,1);
    p_triggers = p_rip_starttime;

    %% 
    cellcountthresh = 0; %set to 0
    riptets = riptetlist;
    
    currriptet=riptets(1);
    % SPIKE COUNT THRESHOLD - MOVE TO ANOTHER FILE?
    filterString = ' (strcmp($tag2, ''CA1Pyr'') || strcmp($tag2, ''CA1Int''))';
    cellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    indices = [repmat([day epoch], size(cellindices,1),1 ), cellindices];
    spikecounts = [];   celldata = [];

    %go through each cell and calculate the binned spike counts
    for cellcount = 1:size(indices,1)
        index = indices(cellcount,:);
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
        else
            spiketimes = [];
        end
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes);
            spikebins = spikebins(validspikes);
        end
        
        if ~isempty(spiketimes)
            tmpcelldata = [spiketimes spikebins];
            tmpcelldata(:,3) = cellcount;
        else
            tmpcelldata = [0 0 cellcount];
        end
        celldata = [celldata; tmpcelldata];
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        
        spikecounts = [spikecounts; spikecount];
    end
    celldata = sortrows(celldata,1); %sort all spikes by time
    cellcounts = sum((spikecounts > 0));
    
    %Find all events with enough cells
    eventindex = find(cellcounts >= cellcountthresh);
    rip_ncells = cellcounts(eventindex);
    riptimes_keep = riptimes(eventindex,:);
    rip_starttime = 1000*riptimes_keep(:,1);  % in ms
    rip_endtime = 1000*riptimes_keep(:,2);
    
    c_rip_starttime = c_rip_starttime*1000;
    c_rip_endtime = c_riptimes(:,2)*1000;
    
    p_rip_starttime = p_rip_starttime*1000;
    p_rip_endtime = p_riptimes(:,2)*1000;
    
    lenidx = find(rip_endtime-rip_starttime > 50);
    rip_starttime = rip_starttime(lenidx);
    rip_endtime = rip_endtime(lenidx);
   
    lenidx2 = find(c_rip_endtime-c_rip_starttime > 50);
    c_rip_starttime = c_rip_starttime(lenidx2);
    c_rip_endtime = c_rip_endtime(lenidx2);
    
    lenidx3 = find(p_rip_endtime-p_rip_starttime > 50);
    p_rip_starttime = p_rip_starttime(lenidx3);
    p_rip_endtime = p_rip_endtime(lenidx3);
   
    triggers = rip_starttime; triggers_end = rip_endtime;
   
    pt = triggers; pt_sec = pt./1000;
    pt_end = triggers_end; pt_endsec = pt_end./1000;
    
    c_pt_sec = c_rip_starttime./1000;
    c_pt_endsec = c_rip_endtime./1000;
    
    p_pt_sec = p_rip_starttime./1000;
    p_pt_endsec = p_rip_endtime./1000;
    
    nrip = length(pt);
    maxcells = max(rip_ncells);
    
    % EEG and ripple data
    % ------------------------
    eegpt = pt./1000; %in secs
    tets=[riptetlist,peegtet];
    for te=1:length(tets)
        
        currtet=tets(te);
        cnteeg=0; cntrip=0;
        
        % Load EEG and ripple LFP file
        %-------------------------
        EEGfile = sprintf('%s/EEG/%seegref%02d-%02d-%02d.mat', directoryname, prefix, day,epoch,currtet);
        load(EEGfile);
        eeg=eegref;
        e = eeg{day}{epoch}{currtet};
        if te==1
            teeg = geteegtimes(e);
            eind = lookup(pt, teeg);
            e.samprate=round(e.samprate);
            eegstart = eeg{day}{epoch}{currtet}.starttime; % in secs - Epoch start
            eegend = teeg(end); % in secs - Epoch end
        end
        ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,currtet);
        cripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,peegtet);
        load(ripfile);
        ripamp = ripple{day}{epoch}{currtet}.data(:,1);
        
        load(cripfile);
        cripamp = double(ripple{day}{epoch}{peegtet}.data(:,1)); %This is ripple amp on PFC tetrode
        
        eegalltet{te} = e.data;
        ripampalltet{te} = double(ripamp);
    end
    
    % ------------------------------
    % Figure Parametersand Font Sizes
    % ------------------------------
    forppr = 0;
    
    set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
    
    if forppr==1
        set(0,'defaultaxesfontsize',16);
        tfont = 18; % title font
        xfont = 16;
        yfont = 16;
    else
        set(0,'defaultaxesfontsize',24);
        tfont = 28;
        xfont = 20;
        yfont = 20;
    end
    clr = {'b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],...
        'r','b','g','y','b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],...
        'r','b','g','y','b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y'};
    
    clr1='k';
    clr2='r';
    
    % Set directory to save rasters
    figdir = '';
    
    %% MANUAL WINDOW SETTING
    winst = 16608.5; %EXAMPLE PLOT FOR PAPER
    winend = 16618.5; % secs
    %% 
    
%     winst = pos{day}{epoch}.data(1,1);
%     winend = winst + win; % secs
    epochend = pos{day}{epoch}.data(end,1);
    ii = 1;
    rastnum = 1;
    while winend <= epochend 
        
        ripidx = find(pt_sec>=winst & pt_sec<=winend);
        riptimes = pt_sec(ripidx);
        
        c_ripidx = find(c_pt_sec>=winst & c_pt_sec<=winend);
        c_riptimes = c_pt_sec(c_ripidx);
        
        p_ripidx = find(p_pt_sec>=winst & p_pt_sec<=winend);
        p_riptimes = p_pt_sec(p_ripidx);
        
        %only plot if ripples in window (this also confines to NREM)
        if (isempty(riptimes)) || (isempty(c_riptimes)) || (isempty(p_riptimes))
            winst = winst + win - overlap;
            winend = winst + win;
            continue
        else
            rastnum
            rastnum = rastnum + 1;
            
            figure(1); xlim([0 win]); hold on;
            redimscreen;
            
            % Get eeg times axis
            eind1 = lookup(winst, teeg);
            eind2 = lookup(winend, teeg);
            taxisEEG = teeg(eind1:eind2);
            taxisEEG = taxisEEG - winst;
            
            taxis = winst:winend; 
            taxis = taxis - winst; % Start axis from 0
            
            winst_ms = winst*1000;
            winend_ms = winend*1000;
            
            % Find ripples within this window
            rend = pt_endsec(ripidx);
            ripend_win = rend - winst;
            riptimes_win = riptimes - winst;
            
            c_rend = c_pt_endsec(c_ripidx);
            c_ripend_win = c_rend - winst;
            c_riptimes_win = c_riptimes - winst;
            
            p_rend = p_pt_endsec(p_ripidx);
            p_ripend_win = p_rend - winst;
            p_riptimes_win = p_riptimes - winst;

            baseline = 0;
            
            % PFC cells
            % ---------------------------------------------
            cnt = 0;
            activepfc = 0;
            [B, I] = sort(cellfun(@length,spiketimep),'descend');
            for c=usecellsp
                cc = I(c);
                eval(['currspkt = spiketimep{',num2str(cc),'};']);
                currspkt = currspkt;
                currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));
                
                if ~isempty(currspkt)
                    % set spike times relative to winst
                    currspkt = currspkt - winst;
                end
                figure(1); xlim([0 win]); hold on; 
                if ~isempty(currspkt)
                    activepfc = activepfc+1;
                    cnt=cnt+1;
                    if size(currspkt,2)~=1
                        currspkt=currspkt';
                    end
                    
                    plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,...
                        [],'Color','r','LineWidth',1);
                end
            end
            
            baseline = baseline + (activepfc)*2;
            baseline = baseline+1;
            
            % CA1 spikes
            % ---------------
            cnt = 0;
            activeca1cnt = 0;
            [B, I] = sort(cellfun(@length,spiketimei),'descend');
            for c=usecellsi
                cc = I(c);
                eval(['currspkt = spiketimei{',num2str(cc),'};']);
                currspkt = currspkt;
                currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));
                
                if ~isempty(currspkt)
                    currspkt = currspkt - winst;
                end
                
                figure(1); hold on; 
                if ~isempty(currspkt)
                    activeca1cnt = activeca1cnt+1;
                    cnt=cnt+1;
                    if size(currspkt,2)~=1
                        currspkt=currspkt'; 
                    end
                    plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,...
                        [],'Color','k','LineWidth',1);
                end
            end
            
            % Update baseline to give a little vertical gap
            baseline = baseline + (activeca1cnt)*2;
            baseline = baseline+1;
            %
            % Then PFC EEG
            % ------------
            n = peegidx;
            eegtet = eegalltet{n};
            curreeg = eegtet(eind1:eind2);
            % Plot
            % ----
            eegscale = max(curreeg)-min(curreeg);
            downeeg = baseline; upeeg = downeeg+8;
            plotscale = 8;
            curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
            plot(taxisEEG,curreeg,'r-','LineWidth',1);
            
            baseline = baseline + 8;
            
            % EEG in PFC ripple band
            % --------------------------------------
            n = maineegidx;
            cripamp_plot = cripamp;
            curreeg = cripamp_plot(eind1:eind2);
            % Plot
            % ----
            eegscale = max(curreeg)-min(curreeg);
            downeeg = baseline; upeeg = downeeg+8;
            plotscale = 8;
            curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
            plot(taxisEEG,curreeg,'r-','LineWidth',1);
            
            % Update baseline
            baseline = baseline + 9;

            %Divider line
            xpts = 0:1:win;
            ypts = baseline*ones(size(xpts));
            plot(xpts , ypts, 'k--','Linewidth',1);
            
            % Update baseline to give a little vertical gap
            baseline = baseline+0.5;
           
            n = maineegidx;
            eegtet = eegalltet{n};
            curreeg = eegtet(eind1:eind2);
            % Plot EEG
            % ----
            eegscale = max(curreeg)-min(curreeg);
            downeeg = baseline; upeeg = downeeg+8;
            plotscale = 8;
            curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
            plot(taxisEEG,curreeg,'k-','LineWidth',1);
            
            % Update baseline
            baseline = baseline + 8;
            
            % EEG in ripple band
            % --------------------------------------
            n = maineegidx;
            riptet = ripampalltet{n};
            curreeg = riptet(eind1:eind2);
            % Plot EEG
            % ----
            eegscale = max(curreeg)-min(curreeg);
            downeeg = baseline; upeeg = downeeg+8;
            plotscale = 8;
            curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
            plot(taxisEEG,curreeg,'k-','LineWidth',1);
            
            %Update baseline
            baseline = baseline + 1;
            
            set(gca, 'YTick',[]);

            xlabel('Time (secs)','FontSize',18,'Fontweight','normal');
            title([prefix ' - Day',num2str(day) 'Ep' num2str(epoch) ' Window Time: '...
                num2str(roundn(winst,-1)) ' - ' num2str(roundn(winend,-1))...
                '  hpctet - ' num2str(maineegtet) ' pfctet - ' num2str(peegtet)],...
                'FontSize',18,'Fontweight','normal');

            ylim = get(gca,'YLim');
            ypts = ylim(1):ylim(2);
            
            % Mark Ripple Times in Current Window
            for s=1:length(riptimes_win)
                DIOt = riptimes_win(s);
                ripend = ripend_win(s);
                xaxis = DIOt:0.005:ripend; % 300ms
                jbfill(xaxis,ylim(2)*ones(size(xaxis)),ylim(1)*ones(size(xaxis)),...
                    'k','k',0.2,0.2);
            end
            
            for s=1:length(c_riptimes_win)
                DIOt = c_riptimes_win(s);
                c_ripend = c_ripend_win(s);
                xaxis = DIOt:0.005:c_ripend; % 300ms
                jbfill(xaxis,ylim(2)*ones(size(xaxis)),ylim(1)*ones(size(xaxis)),'b','b',0.2,0.2);
            end
            
            for s=1:length(p_riptimes_win)
                DIOt = p_riptimes_win(s);
                p_ripend = p_ripend_win(s);
                xaxis = DIOt:0.005:p_ripend; % 300ms
                jbfill(xaxis,ylim(2)*ones(size(xaxis)),ylim(1)*ones(size(xaxis)),'r','r',0.2,0.2);
            end
            
            baseline = baseline+1;
            set(gca,'XLim',[0 winend-winst]);
            set(gca,'YLim',[0 baseline+10]);
            
            keyboard;
            
            saveg=0;
            if saveg==1
                figfile = [figdir,prefix,'Day',num2str(day),'Ep',num2str(epoch),...
                    'RasterEEGegNo',num2str(ii),'Window',num2str(win),'Overlap',num2str(overlap)];
                print('-djpeg', figfile);
            end
            
            ii = ii+1;
            close all
            
            % Update winst and winend
            winst = winst + win - overlap;
            winend = winst + win;
        end
    end
end

keyboard;

