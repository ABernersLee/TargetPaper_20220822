function [corr_repeat_trial_session0,across_lookback0,corr_repeats0] = ...
    get_pv_corr(ns,trialdatsession,dirs,mn,stu,tl,ts,...
    isession,sessionstouse,all_labels,figs_to_run)
 
% Get the correlation across all neurons across all trials
PVall = corrcoef(ns');
PVall2 = PVall;

% changed 4/20/22 after reading more about RSA
PVall2 = zscore(PVall2);

PVall2(eye(size(PVall))==1) = NaN;

%to plot a single session
if all_labels.plotPVsessions && (sum(ismember(figs_to_run,0))>0 || ...
        (sum(ismember(figs_to_run,7))>0 && strcmp('T',mn) && stu==26))
    calculate_plot_PVcorr_plot_crosscorr_sessions...
        (PVall2,trialdatsession,dirs.figdir,mn,stu,tl, all_labels)
end

PVtrial = NaN(size(ns,1),length(ts)+1);
for it = 1:length(ts)+1
     if it>length(ts)                     
         ind1 = trialdatsession(:,1)>2;
     else
        ind1 = trialdatsession(:,1)==ts(it);
     end
     for it2 = 1:length(ts)+1
         if it2>length(ts)                                       
             ind2 = trialdatsession(:,1)>2;
         else
            ind2 = trialdatsession(:,1)==ts(it2);
         end
         look = PVall2(ind1,ind2);
         PVtrial(ind1,it2) = mean(look,2,'omitnan');
     end
 end
 eachtrial_self = NaN(size(PVtrial,1),1);
 for it = 1:length(ts)
     eachtrial_self(trialdatsession(:,1)==ts(it),1) = ...
         PVtrial(trialdatsession(:,1)==ts(it),it);
 end         
 eachtrial_target = PVtrial(:,1);
 eachtrial_nontarget = PVtrial(:,2);         
 eachtrial_others = PVtrial(:,end);
 [~,~,ss] = unique(trialdatsession(:,1));
 

runningTNT = NaN(size(ss));
runningTNT(ss==1) = eachtrial_nontarget(ss==1);
runningTNT(ss==2) = eachtrial_target(ss==2);
linTNT = fillmissing(runningTNT,'movmean',[all_labels.PVCorr_lookback 0]);


 %Correlations broken down by trial type
 corr_repeat_trial_session0 = [eachtrial_self eachtrial_target ...
     eachtrial_nontarget eachtrial_others ...            
            ss (1:size(eachtrial_self,1))' ... %repeat, trial
            ones(size(eachtrial_self))*...
            (sessionstouse(isession,2)-...
            min(sessionstouse(sessionstouse(isession,1) == ...
            sessionstouse(:,1),2))+1) ... %session
            ones(size(eachtrial_self)) * ...
            sessionstouse(isession,1) trialdatsession(:,2) linTNT];
            %mouse, correct, average running corr of T/NT

if length(unique(ss))>3
    nk = nchoosek(unique(ss),2);
    all_vals = NaN(1,size(nk,1));
    for ink = 1:size(nk,1)    
        ind = ss==nk(ink,1);
        all_vals(ink) = mean(PVtrial(ind,nk(ink,2)));
    end
    
     %ABL adding 7/14/22 to look at all combinations of RSA averages 
     % RSA values, trialtype, session, mouse, accuracy
    corr_repeats0 = [all_vals ... %all combinations
                (sessionstouse(isession,2)-...
                min(sessionstouse(sessionstouse(isession,1) == ...
                sessionstouse(:,1),2))+1) ... %session            
                sessionstouse(isession,1)]; %mouse
else
    corr_repeats0 = [];
end
            
            

across_lookback0 = [];
for iLB = 1:length(all_labels.LB)

    PVCorr_lookback = all_labels.LB(iLB);
    runningTNT = NaN(size(ss));
    runningTNT(ss==1) = eachtrial_nontarget(ss==1);
    runningTNT(ss==2) = eachtrial_target(ss==2);
    linTNT = fillmissing(runningTNT,'movmean',[PVCorr_lookback 0]);

    across_lookback0 = cat(2,across_lookback0,linTNT);
end
     