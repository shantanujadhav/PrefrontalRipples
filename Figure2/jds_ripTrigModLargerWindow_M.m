%%------------------------------------------------------------------------
%Justin D. Shin

%Plots ripple modulation for the selected file. Modulation calculated by 
%DFSjds_getripalignspiking
%%------------------------------------------------------------------------
clear all
close all

numAnims = 8; %number of animals
scrsz = get(0,'ScreenSize');
savedirX = ''; %point to directory where modulation files are

plotIndividualCells=0;
areas1={'CA1','PFC'};

figdir = '';

savefig1=0;
set(0,'defaultaxesfontsize',16);

histOns=[];
histOnsW=[];
allxcorrLags={};
allallripmodhists={};
allallripmodhistsW={};
allallripmodhistssig={};
allallripmodhistsWsig={};
allanimdayvecRip={};
allanimdayvecRipW={};
allallripmodhistssigInd={};
allallripmodhiststype={};

allallsiginds={};

% how to smooth psths
b=gaussian(20,61);
animIdx = {};

for areaInd=1:2
    area1=areas1{areaInd};
    switch area1
        
        case 'CA1'
            
            load([savedirX 'Allanim_noncoord250ctxripplemod400_by250mscrit_sleep_CA1INT_alldata_largewin_sepeps_gather_X6.mat'])
            allripplemod_idx=[];
            for w=1:length(allripplemod)
                allripplemod_idx=[allripplemod_idx;allripplemod(w).index];
            end
            
        case 'PFC'
            load([savedirX 'Allanim_noncoord250ctxripplemod400_by250mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6.mat'])
            allripplemod_idx=[];
            for w=1:length(allripplemod)
                allripplemod_idx=[allripplemod_idx;allripplemod(w).index];
            end
    end
    allinds=unique([allripplemod_idx],'rows');
    allripmodhists=[];
    allripmodhistssig=[];
    allripmodhistssigInd=[];
    allripmodonset3=[];
    allripmodhiststype=[];
    allanimidx = [];
    
    allsiginds=[];
    
    for i=1:length(allripplemod)
        allripmodhists=[allripmodhists; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];

        doflag=0;
        switch area1
            case 'CA1'                
                if allripplemod(i).rasterShufP2<0.05
                   doflag=1;
                end
            case 'PFC'
                if allripplemod(i).rasterShufP2<0.05
                    doflag=1;
                end
        end
        if doflag==1
            curranim = allripplemod(i).index(1);
            allanimidx = [allanimidx; curranim];
            allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
            allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
            % 1 for EXC, 0 for INH
            allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
            allsiginds = [allsiginds; allripplemod(i).index];
        end
    end
    
    allallripmodhists{end+1}=allripmodhists;
    allallripmodhistssig{end+1}=allripmodhistssig;
    allallripmodhistssigInd{end+1}=allripmodhistssigInd;
    allallripmodhiststype{end+1}=allripmodhiststype;
    animIdx{end+1} = allanimidx;
    
    cntcellsRip=length(allripplemod);
end
clear allripplemod
CA1psths=allallripmodhistssig{1};
CA1inds=allallripmodhistssigInd{1};
PFCpsths=allallripmodhistssig{2};
PFCinds=allallripmodhistssigInd{2};

%% RATE OF SWR-MODULATED CELLS

CA1sleepSWRmodrate=size(CA1psths,1)./size(allallripmodhists{1},1)
PFCsleepSWRmodrate=size(PFCpsths,1)./size(allallripmodhists{2},1)

%% EXC/INH
excinhPFC=allallripmodhiststype{2};
excinhCA1=allallripmodhiststype{1};

excPFCpsths=PFCpsths(excinhPFC>0,:);
inhPFCpsths=PFCpsths(excinhPFC==0,:);

excCA1psths=CA1psths(excinhCA1>0,:);
inhCA1psths=CA1psths(excinhCA1==0,:);

PFCindsExc=PFCinds(excinhPFC>0,:);
PFCindsInh=PFCinds(excinhPFC==0,:);

%% Plot mean curves (animal separated)
figure
hold on

shadedErrorBar(xaxis,mean(CA1psths,1),std(CA1psths)./sqrt(size(CA1psths,1)),'-b',1);
for a = 1:numAnims
    plot(xaxis, mean(CA1psths(find(animIdx{1}==a),:),1),'-k')
end
xlim([-500 500])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')
figure
hold on
shadedErrorBar(xaxis,mean(PFCpsths,1),std(PFCpsths)./sqrt(size(PFCpsths,1)),'-r',1);
for a = 1:numAnims
    plot(xaxis, mean(PFCpsths(find(animIdx{2}==a),:),1),'-k')
end
xlim([-500 500])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

%% Plot mean EXC/INH modulation for CA1 and PFC
figure
hold on
shadedErrorBar(xaxis,mean(inhCA1psths,1),std(inhCA1psths)./sqrt(size(inhCA1psths,1)),'-k',1);
shadedErrorBar(xaxis,mean(excCA1psths,1),std(excCA1psths)./sqrt(size(excCA1psths,1)),'-r');
xlim([-1000 1000])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
title('CA1 modulation')

figure
hold on
shadedErrorBar(xaxis,mean(inhPFCpsths,1),std(inhPFCpsths)./sqrt(size(inhPFCpsths,1)),'-r',0);
shadedErrorBar(xaxis,mean(excPFCpsths,1),std(excPFCpsths)./sqrt(size(excPFCpsths,1)),'-b');
xlim([-1000 1000])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
title('PFC modulation')

if savefig1 == 1
    print('-dpdf', figfile); 
    print('-dpng', figfile, '-r300'); 
    saveas(gcf,figfile,'fig'); 
    print('-depsc2', figfile); 
    print('-djpeg', figfile);
end

%% Plot sorted modulated cells
% CA1 all ripmod cells
[tmp timingindCA1a]=min(inhCA1psths(:,801:1300)');
[tmp timingindCA1b]=max(excCA1psths(:,801:1300)');
[A1 B1]=sort(timingindCA1a);
[A2 B2]=sort(timingindCA1b);
figure;
subplot(2,1,1)
imagesc(xaxis,1:length([B1';B2']), [inhCA1psths(B1,:); excCA1psths(B2,:)]);caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title(sprintf('CA1 all CA1rip-mod cells - %dINH %dEXC',length(inhCA1psths(:,1)), length(excCA1psths(:,1))))
xlim([-500 500])
hold on
subplot(2,1,2); hold on
boundedline(xaxis,mean(excCA1psths,1),std(excCA1psths)./sqrt(size(excCA1psths,1)),'-r');
boundedline(xaxis,mean(inhCA1psths,1),std(inhCA1psths)./sqrt(size(inhCA1psths,1)),'-b');
xlim([-500 500])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

figure;
subplot(2,1,1)
[tmp timingindPFCa]=min(inhPFCpsths(:,801:1300)');
[tmp timingindPFCb]=max(excPFCpsths(:,801:1300)');
[A1 B1]=sort(timingindPFCa);
[A2 B2]=sort(timingindPFCb);
subplot(2,1,1)
imagesc(xaxis,1:length([B1';B2']), [inhPFCpsths(B1,:); excPFCpsths(B2,:)]);caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title(sprintf('PFC all CA1rip-mod cells - %dINH %dEXC',length(inhPFCpsths(:,1)), length(excPFCpsths(:,1))))
xlim([-500 500])
hold on
subplot(2,1,2); hold on
boundedline(xaxis,mean(excPFCpsths,1),std(excPFCpsths)./sqrt(size(excPFCpsths,1)),'-r');
boundedline(xaxis,mean(inhPFCpsths,1),std(inhPFCpsths)./sqrt(size(inhPFCpsths,1)),'-b');
xlim([-500 500])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

keyboard
if savefig1 == 1
    print('-dpdf', figfile); 
    print('-dpng', figfile, '-r300'); 
    saveas(gcf,figfile,'fig'); 
    print('-depsc2', figfile); 
    print('-djpeg', figfile);
end

%% Peak timing for PFC vs CA1 during events
bins1=-1000:10:1000;
[tmp timingindCA1exc]=max(excCA1psths(:,801:1300)');
[tmp timingindCA1inh]=min(inhCA1psths(:,801:1300)');

[tmp timingindPFCexc]=max(excPFCpsths(:,801:1300)');
[tmp timingindPFCinh]=min(inhPFCpsths(:,801:1300)');

[pfc1 t]=hist(timingindPFCexc-250,bins1);
[pfc2 t]=hist(timingindPFCinh-250,bins1);

[CA11 t]=hist(timingindCA1exc-250,bins1);
[CA12 t]=hist(timingindCA1inh-250,bins1);
figure
plot(bins1,CA12/sum(CA12),'b','linewidth',2); %CA1 inh timing
plot(bins1,pfc1/sum(pfc1),'k','linewidth',2); %PFC exc timing
ylabel('fraction of units')
xlabel('Peak time relative to PFCripple onset')

p1peak = ranksum(timingindCA1inh,timingindPFCexc) %Main timing figure

%% Plot for all cells
figure
plot(bins1,(CA11+CA12)./sum(CA11+CA12),'k','linewidth',2)
hold on;
plot(bins1,(pfc1+pfc2)./sum(pfc1+pfc2),'r','linewidth',2);
ylabel('fraction of units')
xlabel('Peak time relative to PFCripple onset')
p5peak = ranksum([timingindCA1exc timingindCA1inh],[timingindPFCinh timingindPFCexc])
title(sprintf(['CA1 pyr-int IndPFCrip timing - p=' num2str(p5peak)]))
set(gcf, 'renderer', 'painters')


