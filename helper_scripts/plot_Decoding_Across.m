function plot_Decoding_Across(iplottype,decode_fold_label,...
    decode_fold_label_cat,all_labels,dirs,itype, figs_to_run)

%two figures from here are used for figure 6c-d: 
    % Categorical decoding, correct trials only
    % was using the residuals of the anova, but now using the raw 
    % probabilities.
    % The results are the same either way, we thought the raw would be less
    % confusing to people


numshuffhere = 1000;

% Using only iplottope 2 here, where we look at changes across sessions.
% There were no consistent changes across trials in any of these analyses
if iplottype == 1
    plotind = 4; 
    plotlabel = 'Trial'; 
    plotlabeloth = 'Session'; 
    plotind2 = 5; 
    binsz = 10;
elseif iplottype == 2
    plotind = 5; 
    plotlabel = 'Session'; 
    plotlabeloth = 'Trial'; 
    plotind2 = 4; 
    binsz = 1;
end

thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'Decoding_Across_' plotlabel 's/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

for icat = 1:2
    for iraw = 1:2
       for icor = 1:3
          if (sum(ismember(figs_to_run,[0 6]))>0 && icat==2 && icor==2)...
                  || (sum(ismember(figs_to_run,0))>0)
    
              if icat == 1
                  OGdat = decode_fold_label; catlabel = 'Individual';
              elseif icat == 2
                  OGdat = decode_fold_label_cat; catlabel = 'Categorical';
              end
    
             aa = figure; hold on
             bb = figure; hold on
             xl_all = NaN(3,2);
             for itt = 1:3

                 if itt<3
                    indall = decode_fold_label(:,1)==itt;
                 else             
                    indall = decode_fold_label(:,1)>2;
                 end
    
                 if icor ==1
                     corlab = 'AllTrials';
                 elseif icor == 2
                     corlab = 'CorrectOnly';
                     indall = indall & OGdat(:,7)==1;
                 elseif icor == 3             
                     corlab = 'IncorrectOnly';
                     indall = indall & OGdat(:,7)==0;
                 end

                 dat = OGdat(indall,:);
    
                 [~,tbl] = anovan(dat(:,3),dat(:,[1 4:7]),'varnames',...
                    {'Repeat','Trials','Sessions','Mouse','Rewarded'},...
                    'continuous',[2 3],'display','off');
                
                 if iraw == 1
                     datuse = dat(:,3); rawlab = 'Raw';
                 elseif iraw==2    
                     [~,~,stats] = anovan(dat(:,3),...
                        dat(:,[1 plotind2 6:7]),'varnames',...
                        {'Repeat',plotlabeloth,'Mouse','Rewarded'},...
                        'continuous',2,'display','off');    
                     datuse = stats.resid; rawlab = 'Residuals';
                 end

                 if size(datuse,1)==1
                     datuse = datuse';
                 end
    
                 writecell(tbl,[thisdir all_labels.typelab{itype}...
                     '_' catlabel '_' rawlab '_' corlab all_labels.addon...
                     '_ANOVA_' all_labels.ittlabs{itt+1} '.txt']...
                     ,'Delimiter','tab')  
                 
                 writecell(tbl,[thisdir all_labels.typelab{itype}...
                     '_' catlabel '_' rawlab '_' corlab all_labels.addon...
                     '_ANOVA_' all_labels.ittlabs{itt+1} '.xlsx'])    
    
                 
                 figure(aa)
                 subplot(3,1,itt); hold on;

                 if iplottype==1
                     plot(dat(:,plotind),datuse,'o','Color',[.5 .5 .5])
                 elseif iplottype==2
                     violinplot1(datuse(~isnan(datuse)),...
                         dat(~isnan(datuse),plotind),...
                         'EdgeColor',[0 0 0],'BoxColor',[0 0 0],...
                         'ViolinColor',[1 1 1],'ViolinAlpha',0,...
                         'MedianColor',[0 0 0]);
                 end

                 dat1 = ceil((dat(:,plotind)-min(dat(:,plotind))+.001)...
                     /binsz);
    
                 smd = 2;
                 if sum(diff(sort(unique(dat1)))>1)==0
                     m1 = splitapply(@nanmean,datuse,dat1);                        
                     m2 = smoothdata([m1(end:-1:1) m1 m1(end:-1:1)],...
                         'gaussian',smd);
                     m = m2(length(m1)+1:length(m1)*2)';
                     datss = splitapply(@nanstd,datuse,dat1);
                     h = hist(dat1,1:length(m));
                     sem1 = datss./sqrt(h');
    
                     sem2 = smoothdata([sem1(end:-1:1) sem1 ...
                         sem1(end:-1:1)],'gaussian',smd);
                     sem = sem2(length(sem1)+1:length(sem1)*2)';

                     rev = (m+sem)'; fwd = (m-sem)';
                     tbs = (1:length(m))*binsz;
                     p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd ...
                         rev(end:-1:1)],'black');
                    
                     p1.FaceAlpha=.2; 
                     p1.EdgeAlpha=0;    
                     p1.FaceColor = [0 0 0];
                    
                     plot(tbs,m,'k');
                 end

                 [r,p] = corr(dat(:,plotind),datuse,'rows','complete');    
                                                                
                 N = sum(~isnan(datuse) & ~isnan(dat(:,plotind)));
                 tt = title([all_labels.ittlabs{itt+1} ...
                     ', N = ' num2str(N) ', R = ' ...
                     num2str(round(r,2,'significant')) ', p = ' ...
                     num2str(round(p,2,'significant')) ', ANOVAN: ' ...
                     num2str(round(tbl{iplottype+2,7},2,'significant'))]);
                 if (p<.05 || tbl{iplottype+2,7}<0.05) && r>0
                     tt.Color = 'r';
                 elseif (p<.05 || tbl{iplottype+2,7}<0.05) && r<0
                     tt.Color = 'b';
                 end
                 axis tight
                 xlabel(plotlabel)
                 ylabel(all_labels.ittlabs{itt+1})
                 xl = get(gca,'xlim');
                 xlim([xl(1)-1 xl(2)+1])
                
                 r_shuff = NaN(numshuffhere,1);
                 datjnk = dat(:,plotind);

                 for ishuff = 1:numshuffhere                             
                    r_shuff(ishuff) = ...
                        corr(datjnk(randperm(length(datjnk))),...
                        datuse,'rows','complete');    
                 end

                 p = (1+sum(abs(r_shuff)>=abs(r)))./...
                     (1+all_labels.numshuff);   

                 figure(bb)
                 subplot(3,1,itt); hold on
                 histogram(r_shuff,'FaceColor','k')
                 yy = get(gca,'ylim');
                 plot([r r],yy,'r-','LineWidth',3)
                 axis tight
                 xl_all(itt,:) = get(gca,'xlim');
                 xlabel(['Slope of Decoder over ' plotlabel])
                 ylabel('Shuffles')
                 title([all_labels.ittlabs{itt+1} ...
                     ', p = ' num2str(round(p,2,'significant'))])      
                              
             end
    
             figure(aa)
             suptitle({[all_labels.typelab{itype} ' Decoding ' catlabel]...
                 ;['Over ' plotlabel 's, ' rawlab ', ' corlab]})
             set(gcf,'Position',[2206 173 542 695])
             helper_saveandclosefig([thisdir all_labels.typelab{itype}...
                 '_' catlabel '_' rawlab '_' corlab all_labels.addon])
             
            figure(bb)
            for itt = 1:3
                bb.Children(itt).XLim = [min(xl_all(:,1)) ...
                    max(xl_all(:,2))];
            end
    
            suptitle({[catlabel ' Decoding'];[ rawlab ', ' corlab]})      
            set(gcf,'Position',[2206 173 542 695])
            helper_saveandclosefig([thisdir all_labels.typelab{itype}...
                '_' catlabel '_' rawlab '_' corlab all_labels.addon...
                '_Slope'])
          end
      end
   end
end
