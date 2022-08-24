function plot_PV_Sessions(trialdatsession,typelab,mn,dirs,isession,sessionstouse,corr_repeat_trial_session0,thelicktime,all_labels,figs_to_run)

%%% Example population vector correlations and look back calculation 
% for each session
%%% Used in figure 7c,e

[~,~,ss] = unique(trialdatsession(:,1));

runningTNT = NaN(size(ss));
runningTNT(ss==1) = corr_repeat_trial_session0(ss==1,3);
runningTNT(ss==2) = corr_repeat_trial_session0(ss==2,2);
linTNT = fillmissing(runningTNT,'movmean',[all_labels.PVCorr_lookback 0]);

lick_cutoff_ms = 500;
linTNT(thelicktime<=lick_cutoff_ms) = NaN;
thelicktime(thelicktime<=lick_cutoff_ms) = NaN;


thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'PV_Sessions/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end


if sum(ismember(figs_to_run,[0 7]))>0
    figure; hold on; 
    plot(find(ss<=2),runningTNT(ss<=2),'.k'); 
    plot(find(ss>2 &  trialdatsession(:,2)==1),linTNT(ss>2 &...
        trialdatsession(:,2)==1),'ob')
    plot(find(ss>2 &  trialdatsession(:,2)==0),linTNT(ss>2 &...
        trialdatsession(:,2)==0),'xr')
    legend('raw TNT data',...
        'correct probe trial','incorrect probe trial')  
    ylabel('PV correlation T-NT')
    xlabel('Trial Number')
    axis tight
    helper_saveandclosefig([thisdir 'ProbesOverTrials_' ...
        num2str(all_labels.PVCorr_lookback) '_Sessions_' ...
        typelab '_' mn '_' ...
        num2str(sessionstouse(isession,2)) '_summary' all_labels.addon])
end


if sum(ismember(figs_to_run,0))>0
    figure; hold on;
    ind = ss>2 &  trialdatsession(:,2)==1;
    plot(thelicktime(ind),linTNT(ind),'.');        
    [r1,p1] = corr(thelicktime(ind),linTNT(ind),'rows','complete');
    ind = ss>2 &  trialdatsession(:,2)==0;
    plot(thelicktime(ind),linTNT(ind),'.');     
    if sum(ind)>1
        [r0,p0] = corr(thelicktime(ind),linTNT(ind),'rows','complete');
        pp = ranksum(linTNT(ss>2 &  trialdatsession(:,2)==1),linTNT(ss>2 &...
            trialdatsession(:,2)==0),'tail','left');
    else
        r0= NaN;p0 = NaN; pp = NaN;
    end
    xl = get(gca,'xlim'); yl = get(gca,'ylim');
    ind = ss>2;
    [r,p] = corr(thelicktime(ind),linTNT(ind),'rows','complete');      
    
    text(xl(2)-range(xl)/3,yl(1)+(range(yl)/10)*4,...
        ['Ranksum P = ' num2str(round(pp,2,'significant'))])
    text(xl(2)-range(xl)/3,yl(1)+(range(yl)/10)*3,...
        ['r = ' num2str(round(r0,2,'significant')) ', p = ' ...
        num2str(round(p0,2,'significant'))],'Color','r')    
    text(xl(2)-range(xl)/3,yl(1)+(range(yl)/10)*2,['r = ' ...
        num2str(round(r1,2,'significant')) ', p = ' ...
        num2str(round(p1,2,'significant'))],'Color','b')
    text(xl(2)-range(xl)/3,yl(1)+(range(yl)/10)*1,['r = ' ...
        num2str(round(r,2,'significant')) ', p = ' ...
        num2str(round(p,2,'significant'))])
    ylabel('PV correlation T-NT')
    xlabel('Latency (ms)')    
    legend('Correct','Incorrect')
    
    helper_saveandclosefig([thisdir 'LatencyVsPVcorr_' ...
        num2str(all_labels.PVCorr_lookback) '_Sessions_' ...
        typelab '_' mn '_' ...
        num2str(sessionstouse(isession,2)) all_labels.addon])
end