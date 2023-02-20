function plot_TNT_and_Latency(corr_repeat_trial_session,all_labels,dirs,...
    itype,across_lookback,figs_to_run)

numshuff = 1000;
lick_cutoff_ms = 500;

thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'TNT_and_Latency/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

% corr_repeat_trial_session
% 1. PVcorr self, 2. PVcorr target, 3. PVcorr nontarget, 4. PVcorrOthers, 
% 5. trialtype#, 6. trial, 7. session, 8. mouse, 9. accuracy, 
% 10. lookback, 11. lick latency, 12. number of neurons

dat2 = corr_repeat_trial_session(corr_repeat_trial_session(:,5)>2,:);
delind = dat2(:,11)<lick_cutoff_ms;
dat2(delind,:) = [];
pbs = unique(dat2(:,5));

%tested if true for T/NT too, it is

ALB = across_lookback(corr_repeat_trial_session(:,5)>2,:);
ALB(delind,:) = [];

% plot across lookback paramters
% one-sided with the hypothesis that when the repesentation is less 
% similar, the mouse will be more accurate and slower to respond.
Ps = NaN(length(all_labels.LB),1);
for iLB = 1:length(all_labels.LB)
    Ps(iLB,1) = ranksum(ALB(dat2(:,9)==1,iLB),ALB(dat2(:,9)==0,iLB),...
        'tail','left');
end

figure; hold on
set(gcf,'Position',[ 2362         284         455         388])
plot(all_labels.LB,Ps,'ok')
xl = [min(all_labels.LB)-1 max(all_labels.LB)+1];
plot(xl,[0.05 0.05],'--r')
set(gca,'xlim',xl)
ylabel('P-value')
xlabel('Number of trials back to average over')
helper_saveandclosefig([thisdir 'LatencyVsPVcorr_Lookback_' ...
    all_labels.typelab{itype} all_labels.addon])

%Correct trials tend to have longer latencies
for iprobe = 1:length(pbs)+1
    
    if iprobe<=length(pbs) && sum(ismember(figs_to_run,0))==0
        continue
    end

    if iprobe > length(pbs)
        dat = dat2;
    else
        dat = dat2(dat2(:,5)==pbs(iprobe),:);
    end

    figure; hold on; 
    set(gcf,'Position',[ 680    90   997   888])
    
    subplot(2,2,1); hold on; 
    plot(dat(dat(:,9)==1,11),dat(dat(:,9)==1,10),'ob','MarkerSize',2); 
    plot(dat(dat(:,9)==0,11),dat(dat(:,9)==0,10),'xr')
    xl = get(gca,'xlim'); yl = get(gca,'ylim');

    [r1,p1] = corr(dat(dat(:,9)==1,11),dat(dat(:,9)==1,10),...
        'rows','complete');
    [r0,p0] = corr(dat(dat(:,9)==0,11),dat(dat(:,9)==0,10),...
        'rows','complete');
    [r,p] = corr(dat(:,11),dat(:,10),'rows','complete');    

    rps = NaN(numshuff,2);
    for ishuff = 1:numshuff
        [rps(ishuff,1),rps(ishuff,2) ] = ...
            corr(dat(randperm(size(dat,1)),11),...
            dat(randperm(size(dat,1)),10),'rows','complete');
    end

    %two-sided p = 9.99e-4
    Pshuff2 = (1+sum(abs(rps(:,1))>abs(r)))./(1+numshuff); 
     %one-sided p = 9.99e-4
    Pshuff = (1+sum(rps(:,1)<r))./(1+numshuff);

    text(xl(1)+range(xl)/4,yl(1)+(range(yl)/40)*10,['r = ' ...
        num2str(round(r0,2,'significant')) ', p = ' ...
        num2str(round(p0,2,'significant'))],'Color','r')  

    text(xl(1)+range(xl)/4,yl(1)+(range(yl)/40)*8,['r = ' ...
        num2str(round(r1,2,'significant')) ', p = ' ...
        num2str(round(p1,2,'significant'))],'Color','b')

    text(xl(1)+range(xl)/4,yl(1)+(range(yl)/40)*6,['n = ' ...
        num2str(sum(~isnan(dat(:,11)) & ~isnan(dat(:,10)))) ', r = '...
        num2str(round(r,2,'significant')) ', p = '...
        num2str(round(p,2,'significant'))])

    text(xl(1)+range(xl)/4,yl(1)+(range(yl)/40)*4,['Shuffp' ...
        num2str(numshuff) ' = ' num2str(round(Pshuff,2,'significant'))])
    
    ylabel('PV correlation T-NT')
    xlabel('Latency (ms)')    
    legend('Correct','Incorrect')
    
    subplot(2,2,3); hold on
    histogram(dat(dat(:,9)==1,11),xl(1):50:xl(2))
    histogram(dat(dat(:,9)==0,11),xl(1):50:xl(2))
    
    xlabel('Latency (ms)')
    ylabel('Count')    
    yl2 = get(gca,'ylim');
    
    plot([median(dat(dat(:,9)==1,11),'omitnan') ...
        median(dat(dat(:,9)==1,11),'omitnan')],yl2,'b-')
    plot([median(dat(dat(:,9)==0,11),'omitnan') ...
        median(dat(dat(:,9)==0,11),'omitnan')],yl2,'r-')
    
    I1 = dat(dat(:,9)==1,11);
    I2 = dat(dat(:,9)==0,11);
    P = ranksum(I1,I2,'tail','right');

    text(1100,yl2(2)-range(yl2)/2,['Nc = ' ...
        num2str(sum(~isnan(I1))) ', Nin = ' ...
        num2str(sum(~isnan(I2))) ', P = ' ...
        num2str(round(P,3,'significant'))])
    legend('Correct','Incorrect','Correct','Incorrect','Location','Best')
    
    
    subplot(2,2,2); hold on; 
    histogram(dat(dat(:,9)==1,10),(yl(1):.01:yl(2)),...
        'Orientation','horizontal'); 
    histogram(dat(dat(:,9)==0,10),(yl(1):.01:yl(2)),...
        'Orientation','horizontal'); 

    ylabel('PV correlation T-NT')
    xlabel('Count')
    xl2 = get(gca,'xlim'); 
    yl2 = get(gca,'ylim');

    plot(xl2,[median(dat(dat(:,9)==1,10),'omitnan') ...
        median(dat(dat(:,9)==1,10),'omitnan')],'b-')
    plot(xl2,[median(dat(dat(:,9)==0,10),'omitnan') ...
        median(dat(dat(:,9)==0,10),'omitnan')],'r-')
    legend('Correct','Incorrect','Correct','Incorrect')

    I1 = dat(dat(:,9)==1,10);
    I2 = dat(dat(:,9)==0,10);
    P = ranksum(I1,I2,'tail','left');
    
    text(xl2(1)+range(xl2)/4,yl2(1)+range(yl2)/4,['Nc = '...
        num2str(sum(~isnan(I1))) ', Nin = ' ...
        num2str(sum(~isnan(I2))) ', P = ' ...
        num2str(round(P,3,'significant'))])

    %correct-incorrect
    r = median(I1)-median(I2);
    rps = NaN(numshuff,1);
    for ishuff = 1:numshuff
        shufcor = dat(randperm(size(dat,1)),9);
        rps(ishuff,1) = median(dat(shufcor==1,10))...
            -median(dat(shufcor==0,10));        
    end
    %two-sided, p =0.0318
    Pshuff2 = (1+sum(abs(rps(:,1))>abs(r)))./(1+numshuff); 
     %one-sided, p = 0.0155
    Pshuff = (1+sum(rps(:,1)<r))./(1+numshuff);
 
    text(xl2(1)+range(xl2)/4,yl2(1)+range(yl2)/3,['Shuffp' ...
        num2str(numshuff) ' = ' num2str(round(Pshuff,3,'significant'))])
 
    if iprobe > length(pbs)
        helper_saveandclosefig([thisdir 'LatencyVsPVcorr_' ...
            all_labels.typelab{itype} all_labels.addon])
    else
        helper_saveandclosefig([thisdir 'LatencyVsPVcorr_' ...
            all_labels.typelab{itype} all_labels.addon '_probe ' ...
            num2str(iprobe)])
    end
end

