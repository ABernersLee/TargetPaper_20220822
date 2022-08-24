function plot_Selectivity_Types(singledat,selectivity_dat,...
    selectivity_params,dirs,figs_to_run) 
%two figures in figure 4

od = singledat.toplotOdorLick(:,1:6000);
cc = varycolor(3);

sp = selectivity_params;
sd = selectivity_dat;

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'Types/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

%%%%% Time course of neurons' responses to odors by category
lablabs = {'Odor on positive neurons','Odor on negative neurons', ...
    'Odor on positive neurons - Targert',...
    'Odor on positive neurons - NonTarget',...
    'Odor on negative neurons - Targert',...
    'Odor on negative neurons - NonTarget',...
    'Odor on positive neurons - Correct',...
    'Odor on positive neurons - Error',...
    'Odor on negative neurons - Correct',...
    'Odor on negative neurons - Error'};   

for iana = 1:2

    if sum(ismember(figs_to_run,4))>0 && iana>1
        continue
    end

    a = figure; hold on
    b = figure; hold on
    set(gcf,'Position',[1958        -307        1003        1633])
    c = figure; hold on
    for ilablab = 1:size(lablabs,2)
        if ilablab==1
            n = find(sd.all_Neurons(:,3)); 
        elseif ilablab==2
            n = find(sd.all_Neurons(:,4));
        elseif ilablab==3
            n = find(sd.all_Neurons(:,3) & sd.all_Neurons(:,5)); 
        elseif ilablab==4
            n = find(sd.all_Neurons(:,3) & sd.all_Neurons(:,6));
        elseif ilablab==5
            n = find(sd.all_Neurons(:,4) & sd.all_Neurons(:,5)); 
        elseif ilablab==6
            n = find(sd.all_Neurons(:,4) & sd.all_Neurons(:,6));
        elseif ilablab==7  
            n = find(sd.all_Neurons(:,3) & sd.all_Neurons(:,7)); 
        elseif ilablab==8
            n = find(sd.all_Neurons(:,3) & sd.all_Neurons(:,8));
        elseif ilablab==9  
            n = find(sd.all_Neurons(:,4) & sd.all_Neurons(:,7)); 
        elseif ilablab==10
            n = find(sd.all_Neurons(:,4) & sd.all_Neurons(:,8));
        end
        
        ntotal = size(sd.neuron_info,1);
        if iana==2
            n = n(ismember(n,find(sd.neuron_info(:,7)==1)));
            ntotal = sum(sd.neuron_info(:,7)==1);
        end
        
        if length(n)<2
            continue
        end
        
        dat = od(n,:);

        if ~ismember(ilablab,[2,5:6,9:10])
            [~,mm] = max(dat(:,1000:3000),[],2);
        else
            [~,mm] = min(dat(:,1000:3000),[],2);        
        end
        [~,ord] = sort(mm);

        dat2 = zscore(dat,[],2);
        datBS = dat2-mean(dat2(:,1:1000),2);

        figure(a)
        subplot(5,2,ilablab); hold on
        imagesc(dat(ord,:)./max(dat(ord,:),[],2))
        set(gca,'ytick',1:size(datBS,1),'yticklabel',n(ord))
        set(gca,'xtick',1000:1000:6000,'xticklabel',0:5)
        colorbar
        axis tight ij            
        title([lablabs{ilablab} ', ' ...
            num2str(round(100*length(n)./ntotal,2,'significant')) '%']) 

        figure(b)
        subplot(5,2,ilablab); hold on    
        tbs = 1:size(datBS,2);
        m1 = mean(datBS); 
        sem1 = std(datBS)./sqrt(size(datBS,1));
        m = smoothdata(m1,'movmean',sp.sm); 
        sem = smoothdata(sem1,'movmean',sp.sm);
        rev = m+sem; fwd = m-sem;
        p1 = patch([tbs tbs(length(tbs):-1:1)],...
            [fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; 
        p1.EdgeAlpha=0;    
        p1.FaceColor = [0 0 0];
        plot(mean(datBS,'omitnan'),'k');
        yl = get(gca,'ylim');
        plot([1000 1000],yl,'r')
        plot([3000 3000],yl,'r')
        set(gca,'xtick',1000:1000:6000,'xticklabel',0:5)    
        title([lablabs{ilablab} ', n = ' num2str(length(n))]) %num2str(round(100*length(n)./ntotal,2,'significant')) '%']) 

        if ismember(ilablab,2:4)
            figure(c)        
            subplot(3,1,1); hold on
            plot(mean(zscore(dat,[],2),'omitnan'),...
                '-','color',cc(ilablab-1,:))
            plot(mean(zscore(dat,[],2),'omitnan')+...
                std(zscore(dat,[],2))./sqrt(size(dat,1)),...
                '--','color',cc(ilablab-1,:))
            plot(mean(zscore(dat,[],2),'omitnan')-...
                std(zscore(dat,[],2))./sqrt(size(dat,1)),...
                '--','color',cc(ilablab-1,:))
            set(gca,'xtick',[1000:1000:6000],'xticklabel',0:5)    
            ylabel('Zscored')

            subplot(3,1,2); hold on
            plot(mean(dat-dat(:,1),'omitnan'),'-','color',cc(ilablab-1,:))
            plot(mean(dat-dat(:,1),'omitnan')+...
                std(dat-dat(:,1))./...
                sqrt(size(dat,1)),'--','color',cc(ilablab-1,:))
            plot(mean(dat-dat(:,1),'omitnan')-std(dat-dat(:,1))./...
                sqrt(size(dat,1)),'--','color',cc(ilablab-1,:))
            set(gca,'xtick',1000:1000:6000,'xticklabel',0:5)    
            ylabel('Baseline subtracted')

            subplot(3,1,3); hold on
            plot(mean(dat-mean(dat,2),'omitnan')...
                ,'-','color',cc(ilablab-1,:))
            plot(mean(dat-mean(dat,2),'omitnan')+...
                std(dat-mean(dat,2))./sqrt(size(dat,1))...
                ,'--','color',cc(ilablab-1,:))
            plot(mean(dat-mean(dat,2),'omitnan')-....
                std(dat-mean(dat,2))./sqrt(size(dat,1))...
                ,'--','color',cc(ilablab-1,:))
            set(gca,'xtick',1000:1000:6000,'xticklabel',0:5)    
            ylabel('Mean subtracted')
        end
    end

    
    if sum(ismember(figs_to_run,[0 4]))>0
        figure(a)
        set(gcf,'Position',[1958 -307 1003 1633])
        helper_saveandclosefig([thisdir 'Imagesc_Neurons' ...
            num2str(ntotal) '_pvalue' num2str(sp.p_cutoff) '_' ...
            sp.ana_labs{iana}])
    
        figure(b)
        helper_saveandclosefig([thisdir 'PlotAvg_Neurons' ...
            num2str(ntotal) '_pvalue' num2str(sp.p_cutoff) '_' ...
            sp.ana_labs{iana}])
    end

    figure(c)
    if sum(ismember(figs_to_run,0))>0
        set(gcf,'Position',[1958 -307 1003 1633])
        helper_saveandclosefig([thisdir 'PlotAvgSubset_Neurons' ...
            num2str(ntotal) '_pvalue' num2str(sp.p_cutoff) '_' ...
            sp.ana_labs{iana}])
    else
        close gcf
    end
end

