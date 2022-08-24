function plot_Decoding_SigShuffle_FromShuffle(decode_fold_shuffle, ...
    decode_fold_label, decode_fold_cat_shuffle,decode_fold_label_cat,...
    dirs,all_labels,decode_fold_shuffle_prob,...
    decode_fold_cat_shuffle_prob,itype,figs_to_run)


% currently using only 1 of these figures for the paper (figure 6a):
    % Categorical decoding
    % All trials (not just correct or just incorrect) and not using the 
        % plot of the shuffle distribution with the true value of trial 
        % types averaged
%(The booleans will only generate one if figs_to_run does not have zero.)


thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'Decoding_SigShuffle_FromShuffle/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

trialts = 1:3;
ps = NaN(length(trialts)+2,2);
numshuff = size(decode_fold_shuffle,2);

%get p values with Monte-Carlo
for tt = 1:length(trialts)
    ps(tt,1) = (sum(mean(decode_fold_shuffle(...
        decode_fold_label(:,1)==trialts(tt),:),'omitnan') >= ...
        mean(decode_fold_label(decode_fold_label(:,1) == ...
        trialts(tt),2),'omitnan'))+1)./(numshuff+1);
    ps(tt,2) = (sum(mean(decode_fold_cat_shuffle(...
        decode_fold_label(:,1)==trialts(tt),:),'omitnan') >= ...
        mean(decode_fold_label_cat(decode_fold_label(:,1) == ...
        trialts(tt),2),'omitnan'))+1)./(numshuff+1);
end
ps(tt+1,1) = (sum(mean(decode_fold_shuffle,'omitnan') >= ...
    mean(decode_fold_label(:,2)),'omitnan')+1)./(numshuff+1);
ps(tt+1,2) = (sum(mean(decode_fold_cat_shuffle,'omitnan') >= ...
    mean(decode_fold_label_cat(:,2),'omitnan'))+1)./(numshuff+1);
ps(tt+2,1) = (sum(mean(decode_fold_shuffle(decode_fold_label(:,1)>2,:)...
    ,'omitnan') >= mean(decode_fold_label(decode_fold_label(:,1)>2,2),...
    'omitnan'))+1)./(numshuff+1);
ps(tt+2,2) = (sum(mean(decode_fold_cat_shuffle(decode_fold_label(:,1) >...
    2,:),'omitnan') >= mean(decode_fold_label_cat(decode_fold_label(:,1)...
    > 2,2),'omitnan'))+1)./(numshuff+1);

%Improvement over shuffle
Improvement = [decode_fold_label(:,3)-mean(decode_fold_shuffle_prob,2) ... 
    decode_fold_label_cat(:,3)-mean(decode_fold_cat_shuffle_prob,2)];
RawDec = [decode_fold_label(:,2) decode_fold_label_cat(:,2)];

[~,~,u] = unique(decode_fold_label(:,5)+100*decode_fold_label(:,6));
udat = NaN(4,max(u)); udat2 = udat;
for iu = 1:max(u)
    inddat = u==iu  & decode_fold_label(:,7)==1; %only correct now
    udat(:,iu) = [mean(Improvement(inddat,:),'omitnan') ...
        mean(decode_fold_label(inddat,8),'omitnan') ...
        unique(decode_fold_label(inddat,6))];            
    udat2(:,iu) = [mean(RawDec(inddat,:),'omitnan') ...
        mean(decode_fold_label(inddat,8),'omitnan') ...
        unique(decode_fold_label(inddat,6))];
end

for idd = 1:2
    for irew = 1:3
        if (sum(ismember(figs_to_run,[0 6]))>0 && irew~=2 && idd==2) ...
                || (sum(ismember(figs_to_run,0))>0)

            if idd == 1
                d1 = decode_fold_label; d2 = decode_fold_shuffle; 
                ch = length(trialts);                
                catlabel = 'Individual'; xxl = [0 1];
            else
                d1 = decode_fold_label_cat; d2 = decode_fold_cat_shuffle; 
                ch = 2;
                catlabel = 'Categorical'; xxl = [0 1];
            end
    
            if irew==1
                rewlab = '';
            elseif irew==2
                rewlab = '_rRewardedOnly';
                d1(decode_fold_label(:,7)~=1,:) = [];                    
                d2(decode_fold_label(:,7)~=1,:) = [];
            elseif irew == 3
                rewlab = '_rIncorrectOnly';
                d1(decode_fold_label(:,7)~=0,:) = [];                    
                d2(decode_fold_label(:,7)~=0,:) = [];
            end
            
            aa = figure; hold on        
            for itt = 1:3
                if itt<3
                    dat = d1(d1(:,1)==trialts(itt),2);
                    [~,~,dat2] = unique(d1(d1(:,1)==trialts(itt),6));
                    shuffdat = d2(d1(:,1)==trialts(itt),:);
                    lab1 = all_labels.lablabs{itt}; 
                else
                    dat = d1(d1(:,1)>2,2);
                    %mouse for each trial
                    [~,~,dat2] = unique(d1(d1(:,1)>2,6)); 
                    shuffdat = d2(d1(:,1)>2,:);
                    lab1 = all_labels.typelab{itype}; 
                end
                p = (sum(mean(d2(d1(:,1)==trialts(itt),:),'omitnan') >= ...
                    mean(d1(d1(:,1)==trialts(itt),2),'omitnan'))+1)...
                    ./(numshuff+1);
    
                figure(aa)
                Violin(mean(shuffdat,'omitnan'),itt,'ViolinColor',[0 0 0]);
                %take the mean for each mouse (of all trials)
                dat22 = splitapply(@mean,dat,dat2); 
                errorbar(itt,mean(dat22),std(dat22)./...
                    sqrt(sum(~isnan(dat22))),'r')
                plot(itt,mean(dat22),'.r','MarkerSize',20)            
                
                if sum(ismember(figs_to_run,0))>0 || irew==3
                    figure; hold on
                    histogram(mean(shuffdat,'omitnan'),'FaceColor','w')
                    yl = get(gca,'ylim');
                    plot([mean(dat,'omitnan') mean(dat,'omitnan')]...
                        ,yl,'k','LineWidth',2)    
                    text(mean(dat,'omitnan'),yl(2)-range(yl)*.5,...
                        ['M = ' num2str(mean(dat,'omitnan'))])
                    plot([1/ch 1/ch],yl,'k--','LineWidth',2)
                    legend({'Shuffle control','Real data','Chance'},...
                        'Location','Best')
                    ylabel('Number of shuffles')
                    xlabel('Proportion correct')
                    xlim(xxl)
                    title([catlabel ', ' lab1 ' trials, p = ' ...
                        num2str(round(p,2,'significant')) rewlab])
        
                    helper_saveandclosefig([thisdir...
                        all_labels.typelab{itype} ...
                        '_Histogram_' lab1 ...
                        '_' catlabel '_' ...
                        all_labels.addon rewlab])      
                end
            end
            
            m = mean(d1(d1(:,1)<4,2),'omitnan');
            n = sum(~isnan(d1(d1(:,1)<4,2)));
            p = (sum(mean(d2(d1(:,1)<4,:),'omitnan') >= ...
                    m)+1)...
                    ./(numshuff+1);
    
            figure(aa)
            axis tight
            set(gca,'xtick',1:3,'xticklabel',all_labels.lablabs)
            ylabel('Proportion correct')
            title([catlabel ' Decoding , p = ' ...
                        num2str(round(p,3,'significant')) ', m = ' ...
                        num2str(m)  ', n = ' num2str(n) rewlab])
            
            helper_saveandclosefig([thisdir...
                    all_labels.typelab{itype} ...
                    '_Histogram_All_' ...
                    catlabel '_' ...
                    all_labels.addon rewlab])       
        end
    end
  
    if sum(ismember(figs_to_run,0))>0
        % Improvement
        [r1,p1] = corr(udat(idd,:)',udat(3,:)','rows','complete');  
        %out of the correct ones now
        [r2,p2] = corr(Improvement(decode_fold_label(:,7)==1,idd), ...
            decode_fold_label_cat(decode_fold_label(:,7)==1,8), ...
            'rows','complete');
        [~,tbl] = anovan(Improvement(:,idd),decode_fold_label(:,[1 4:8]) ...
            ,'varnames',{'TrialType','Trial','Session','Mouse','Rewarded','RT'} ...
            ,'continuous',[2,3,6],'display','off');  
        p3 = tbl{7,7}; %Reaction time
        [~,~,stats] = anovan(Improvement(:,idd),decode_fold_label(:,[1 4:7]), ...
            'varnames',{'TrialType','Trial','Session','Mouse','Rewarded'}, ...
            'continuous',[2,3],'display','off');    
        datuse = stats.resid;                      
        if size(datuse,1) == 1
            datuse= datuse';
        end
        [r4,p4] = corr(datuse,decode_fold_label(:,8),'rows','complete');
        
        figure; hold on
        for imouse = 1:max(udat(4,:))
            plot(udat(idd,udat(4,:)==imouse),udat(3,udat(4,:)==imouse),'o')
        end
        xlabel(['Decoding ' catlabel ', improvement over shuffle'])
        ylabel('Lick reaction time (ms)')
        yl = get(gca,'ylim'); xl = get(gca,'xlim');
        title([catlabel ' , Residual R = ' num2str(round(r4,2,'significant')) ...
            ', P = ' num2str(round(p4,2,'significant')) ', AnovaP = ' ...
            num2str(round(p3,2,'significant'))])
        text(xl(1)+range(xl)*.3,yl(1)+range(yl)*.9,['Sessions: R = ' ...
            num2str(round(r1,2,'significant')) ', P = ' ...
            num2str(round(p1,2,'significant'))])
        text(xl(1)+range(xl)*.3,yl(1)+range(yl)*.8,['All: R = ' ...
            num2str(round(r2,2,'significant')) ', P = ' ...
            num2str(round(p2,2,'significant'))]) 
    
        helper_saveandclosefig([thisdir all_labels.typelab{itype} ...
        '_DecodingImprovement_ReactionTimes_' catlabel '_' all_labels.addon])  
        
        
        figure; hold on
        for imouse = 1:max(udat(4,:))
            plot(Improvement(decode_fold_label(:,7)==1 & u==imouse,idd), ...
                decode_fold_label(decode_fold_label(:,7)==1  & u==imouse,8) ...
                ,'o')
        end
        xlabel(['Decoding ' catlabel ', improvement over shuffle'])
        ylabel('Lick reaction time (ms)')
        yl = get(gca,'ylim'); xl = get(gca,'xlim');
        title([catlabel ' , Residual R = ' num2str(round(r4,2,'significant')) ...
            ', P = ' num2str(round(p4,2,'significant')) ', AnovaP = ' ...
            num2str(round(p3,2,'significant'))])
        text(xl(1)+range(xl)*.3,yl(1)+range(yl)*.9,['Sessions: R = ' ...
            num2str(round(r1,2,'significant')) ', P = ' ...
            num2str(round(p1,2,'significant'))])
        text(xl(1)+range(xl)*.3,yl(1)+range(yl)*.8,['All: R = ' ...
            num2str(round(r2,2,'significant')) ', P = ' ...
            num2str(round(p2,2,'significant'))])    
         helper_saveandclosefig([thisdir all_labels.typelab{itype} ...
                    '_DecodingImprovement_ReactionTimes_EachTrial_' ...
                    catlabel '_' all_labels.addon])  
        
        
        %Raw decoding accuracy
        [r1,p1] = corr(udat2(idd,:)',udat2(3,:)','rows','complete');   
        %out of the correct ones now
        [r2,p2] = corr(RawDec(decode_fold_label(:,7)==1,idd), ...
            decode_fold_label(decode_fold_label(:,7)==1,8),'rows','complete');
        [~,tbl] = anovan(RawDec(:,idd),decode_fold_label(:,[1 4:8]), ...
            'varnames',{'TrialType','Trial','Session','Mouse','Rewarded','RT'}, ...
            'continuous',[2,3,6],'display','off');  
        p3 = tbl{7,7}; %reaction time
        [~,~,stats] = anovan(RawDec(:,idd),decode_fold_label(:,[1 4:7]) ...
            ,'varnames',{'TrialType','Trial','Session','Mouse','Rewarded'}, ...
            'continuous',[2,3],'display','off');    
        datuse = stats.resid;         
        if size(datuse,1) == 1
            datuse= datuse';
        end
        [r4,p4] = corr(datuse,decode_fold_label(:,8),'rows','complete');
        
        figure; hold on
        for imouse = 1:max(udat2(4,:))
            plot(udat2(idd,udat2(4,:)==imouse),udat2(3,udat2(4,:)==imouse),'o')
        end
        xlabel(['Decoding ' catlabel ', improvement over shuffle'])
        ylabel('Lick reaction time (ms)')
        yl = get(gca,'ylim'); xl = get(gca,'xlim');
        title([catlabel ' , Residual R = ' ...
            num2str(round(r4,2,'significant')) ', P = ' ...
            num2str(round(p4,2,'significant')) ', AnovaP = ' ...
            num2str(round(p3,2,'significant'))])
        text(xl(1)+range(xl)*.3,yl(1)+range(yl)*.9,['Sessions: R = ' ...
            num2str(round(r1,2,'significant')) ', P = ' ...
            num2str(round(p1,2,'significant'))])
        text(xl(1)+range(xl)*.3,yl(1)+range(yl)*.8,['All: R = ' ...
            num2str(round(r2,2,'significant')) ', P = ' ...
            num2str(round(p2,2,'significant'))])        
         helper_saveandclosefig([thisdir all_labels.typelab{itype} ...
                    '_DecodingRaw_ReactionTimes_' catlabel...
                    '_' all_labels.addon])  
        
        
        figure; hold on
        for imouse = 1:max(udat2(4,:))
            plot(RawDec(decode_fold_label(:,7)==1 & u==imouse,idd), ...
                decode_fold_label(decode_fold_label(:,7)==1  & u==imouse,8), ...
                'o')
        end
        xlabel(['Decoding ' catlabel ', improvement over shuffle'])
        ylabel('Lick reaction time (ms)')
        yl = get(gca,'ylim'); xl = get(gca,'xlim');
        title([catlabel ' , Residual R = ' ...
            num2str(round(r4,2,'significant')) ', P = ' ...
            num2str(round(p4,2,'significant')) ', AnovaP = ' ...
            num2str(round(p3,2,'significant'))])
        text(xl(1)+range(xl)*.3,yl(1)+range(yl)*.9,['Sessions: R = ' ...
            num2str(round(r1,2,'significant')) ', P = ' ...
            num2str(round(p1,2,'significant'))])
        text(xl(1)+range(xl)*.3,yl(1)+range(yl)*.8,['All: R = ' ...
            num2str(round(r2,2,'significant')) ', P = ' ...
            num2str(round(p2,2,'significant'))])
         helper_saveandclosefig([thisdir all_labels.typelab{itype} ...
                    '_DecodingRaw_ReactionTimes_EachTrial_' catlabel '_' ...
                    all_labels.addon])  
    end
                   

end
