function plot_Selectivity_proportions(singledat,selectivity_params,...
    dirs,figs_to_run)

% 2 figs for figure 5, one for supplemental figure 6

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'Selectivity_proportions/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

sp = selectivity_params;

isProbeRepeatSession(:,:,1) = ...
    any(squeeze(sum(sum(~isnan(...
    singledat.odorselective(:,:,:,:,1,15:17)),4),3))>0,3);
isProbeRepeatSession(:,:,2) = ...
    any(squeeze(sum(sum(~isnan(...
    singledat.odorselective(:,:,:,:,1,18:19)),4),3))>0,3);


for itype = 1:2
    
    A1 = squeeze(singledat.odorselective_all(:,:,:,1,2));       
    A11 = sum(A1<sp.p_cutoff,3);
    indd = sum(~isnan(A1),3);
    sess = repmat(1:size(A11,2),[size(A11,1) 1 size(A11,3)]);        
    ani = repmat(1:size(A11,1),[size(A11,2) 1 size(A11,3)])';
    datA = A11./indd;    
    datA(isProbeRepeatSession(:,:,itype)==0) = NaN;    
    [p2,tbl2,~] = anovan(datA(:),[sess(:) ani(:)],'varnames', ...
        {'Session';'Mouse'},'continuous',1,'display','off');
        
  
    %figure 5e
    if (sum(ismember(figs_to_run,[0 5]))>0 && itype==2) || ...
            (sum(ismember(figs_to_run,0))>0)
        figure; hold on
        for imouse = 1:4
            if sum(~isnan(datA(imouse,:)))>0
                plot(sess(imouse,:),datA(imouse,:),'o-');
            end
        end
        [xF,yF] = helper_plotbestfitline(sess(:),datA(:));
        plot(xF,yF,'Color','k','LineWidth',2)
        [r,p] = corr(sess(:),datA(:),'rows','complete');
        N = sum(~isnan(sess(:)) & ~isnan(datA(:))); 
        axis tight
        xl = get(gca,'xlim');
        set(gca,'xlim',[xl(1)-1 xl(2)+2])
        xlabel('Sessions'); 
        ylabel('Proportion of Target/NT mod. neurons')
           suptitle({['Prop Corr Sig; N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{itype}];...
            ['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})         
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' ...
            sp.rep_probe{itype} ...
            '_OverallP_mouseplotline'])   
         writecell(tbl2,[thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' ...
            sp.rep_probe{itype} ...
            '_OverallP_mouseplotline_ANOVAN.txt'],'Delimiter','tab')  
         writecell(tbl2,[thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' ...
            sp.rep_probe{itype} ...
            '_OverallP_mouseplotline_ANOVAN.xlsx'])  
    end

end



%%%%% paired comparison of trial downsampled target vs probe or repeat    
repnum = [5 7 10 12 15 17;8 9 13 14 18 19];
allp = NaN(size(singledat.odorselective,1),...
    size(singledat.odorselective,2),2,2,3); 
%probe or repeat / target, whether probe or repeat, # of odor
for irep = 1:2

    allp2 = NaN(size(singledat.odorselective,1),...
        size(singledat.odorselective,2),3,3); 

    A1 = squeeze(singledat.odorselective...
        (:,:,:,sp.OdorIndSmall,2,repnum(irep,5):repnum(irep,6)));        
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %supplemental figure 6 - paired downsample
    if sum(ismember(figs_to_run,[0 16]))>0
        A1 = squeeze(singledat.odorselective...
            (:,:,:,sp.OdorIndSmall,2,repnum(irep,1):repnum(irep,2)));
        A2 = squeeze(singledat.odorselective...
            (:,:,:,sp.OdorIndSmall,2,repnum(irep,3):repnum(irep,4)));
        figa = figure; hold on
        figb = figure; hold on
        aa = NaN(size(A1,5),1); ylaa = aa;
        for ip = 1:size(A1,5)            
            inds = sum(sum(~isnan(A1(:,:,:,:,ip)),4)~=0,3);
            if sum(sum(sum(inds)))>0
                A11 = squeeze(sum(any(A1(:,:,:,:,ip)<...
                    (sp.p_cutoff/sum(sp.OdorIndSmall)),4),3)); 
                A22 = squeeze(sum(any(A2(:,:,:,:,ip)<...
                    (sp.p_cutoff/sum(sp.OdorIndSmall)),4),3)); 
                allp(:,:,1,irep,ip) = [A11./inds];
                allp(:,:,2,irep,ip) = [A22./inds];

                A11s = A11; A22s = A22;
                A11s(inds==0) = NaN; A22s(inds==0) = NaN;
                allp2(:,:,1,ip) = A11s;
                allp2(:,:,2,ip) = A22s;
                allp2(:,:,3,ip) = inds;

                figure(figa)
                aa(ip) = subplot(1,size(A1,5),ip); hold on
                for imouse = 1:4
                   plot(imouse,allp(imouse,:,1,irep,ip),'ok')                   
                   plot(imouse+6,allp(imouse,:,2,irep,ip),'or')
                end
                errorbar(2.5,mean(mean(allp(:,:,1,irep,ip)...
                    ,2,'omitnan'),1,'omitnan'),...
                    std(mean(allp(:,:,1,irep,ip),2,'omitnan')...
                    ,1,'omitnan')./sqrt(size(allp,1)),'k')                
                errorbar(8.5,mean(mean(allp(:,:,2,irep,ip)...
                    ,2,'omitnan'),1,'omitnan'),...
                    std(mean(allp(:,:,2,irep,ip),2,'omitnan')...
                    ,1,'omitnan')./sqrt(size(allp,1)),'r')
                set(gca,'xtick',[2.5 8.5],'xticklabel'...
                    ,{sp.rep_probe{irep},'Target'})      
                yl = get(gca,'ylim');
                ylaa(ip) = yl(2);
                xlim([0 11])
                ylabel('Proportion of neurons significant')                
                I1 = A11(:)./inds(:); I2 = A22(:)./inds(:);
                p = signrank(I1,I2,'tail','left');
                xlabel(['One-sided signrank: N =' ...
                    num2str(sum(~isnan(I2))) ', P = ' ...
                    num2str(round(p,2,'significant'))])

                figure(figb); subplot(1,size(A1,5),ip); hold on
                plot(1:2,[sum(A11,2) sum(A22,2)],'.-k','MarkerSize',20)  
                set(gca,'xtick',1:2,'xticklabel',...
                    {sp.rep_probe{irep},'Target'}) 
                xlim([.5 2.5])
                p2 = signrank(sum(A11,2),sum(A22,2),'tail','left');
                xlabel(['One-sided signrank: N =' ...
                    num2str(size(A22,1)) ', P = ' ...
                    num2str(round(p2,2,'significant'))])
            end
            title([sp.rep_probe{irep} ' ' num2str(ip) ' vs Target'])
        end

         figure(figa)
        for ip = 1:size(A1,5)            
            set(aa(ip),'ylim',[0 max(ylaa)])
        end
        
        A11 = squeeze(allp(:,:,1,irep,:));
        A22 = squeeze(allp(:,:,2,irep,:));
        
        p = signrank(A11(:),A22(:),'tail','left');
        suptitle(['One-sided signrank: N =' ...
            num2str(sum(~isnan(A11(:)))) ', P = ' ...
            num2str(round(p,2,'significant'))])
        set(gcf,'Position',[1999 321 965 398])
        helper_saveandclosefig([thisdir 'PairedDownsample_' ...
            sp.ana_labs{3} '_' sp.rep_probe{irep}])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %figure 5d
    if (sum(ismember(figs_to_run,[0 5]))>0 && irep==2) || ...
            (sum(ismember(figs_to_run,0))>0)
         inds = sum(any(any(~isnan(A1(:,:,:,:,:))~=0,4),5),3);
         A11 = sum(any(any(A1(:,:,:,:,:) < ...
             (sp.p_cutoff/sum(sp.OdorIndSmall)),4),5),3);         
         dat = A11./inds;
         dat2 = repmat(1:size(A11,2),[size(A11,1) 1]);
       
        ani = repmat([1:size(A11,1)]',[1 size(A11,2) size(A11,3)]);
        [p2,tbl2,~] = anovan(dat(:),[dat2(:) ani(:)],...
            'varnames',{'Session';'Mouse'},'continuous',1,'display','off');
    
        figure; hold on;
        for imouse = 1:4
            if sum(~isnan(dat(imouse,:)))>0
                plot(dat2(imouse,:),dat(imouse,:),'o-');
            end
        end
        [xF,yF] = helper_plotbestfitline(dat2(:),dat(:));    
        plot(xF,yF,'Color','k','LineWidth',2)
        [r,p] = corr(dat(:),dat2(:),'rows','complete');
         ylabel(['Proportion of neurons significantly ' ...
             sp.rep_probe{irep} ' modulated'])
         xlabel('Sessions')
         N = sum(~isnan(dat(:)) & ~isnan(dat2(:))); 
         axis tight
         xl = get(gca,'xlim');
         set(gca,'xlim',[xl(1)-1 xl(2)+1])
         suptitle({['N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{irep}];...
            ['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})   
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_mouse'])  
        writecell(tbl2,[thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_' sp.rep_probe{irep} ...
            '_AnyChance_mouse_ANOVAN.txt'],'Delimiter','tab')    
        writecell(tbl2,[thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_' sp.rep_probe{irep} ...
            '_AnyChance_mouse_ANOVAN.xlsx'])    
    end
end