function plot_Selectivity_OdorOff(singledat,selectivity_dat,...
    selectivity_params,dirs,figs_to_run)
% one panel for Figure s5

od = singledat.toplotOdorLick(:,1:6000);


sp = selectivity_params;
sd = selectivity_dat;

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'OdorOff/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

  %%%%% Odor-off responses
lablabs = {'Odor on positive neurons','Odor on negative neurons', ...
    'Odor on positive neurons and odor off',...
    'Odor on negative neurons and odor off','Odor-off only'};

if sum(ismember(figs_to_run,[0 15]))>0
    for iana = 1        
        alldat = cat(2,sd.all_Neurons(:,3:4), sd.all_Neurons(:,3) & ...
            sd.all_Neurons(:,13), sd.all_Neurons(:,4) & ...
            sd.all_Neurons(:,13),~sd.all_Neurons(:,3) & ...
            ~sd.all_Neurons(:,4) & sd.all_Neurons(:,13));        
        
        
         ntotal = size(sd.neuron_info,1); 
        
        if iana==2    
            ntotal = sum(sd.neuron_info(:,7));
            alldat = alldat & repmat(sd.neuron_info(:,7)==1,...
                [1 size(alldat,2)]);
        end
        
        figure; hold on
        for idat = 1:5
            subplot(5,2,(idat-1)*2+1); hold on
            n = find(alldat(:,idat));             
            title([lablabs{idat} ', ' num2str(length(n))]) 
            if length(n)>1
                dat = od(n,:);

                dat2 = zscore(dat,[],2);
                datBS = dat2-mean(dat2(:,1:1000),2);

                [~,mm] = max(dat(:,3000:-1:1000),[],2);
                [~,ord] = sort(mm);
                imagesc(dat(ord,:)./max(dat(ord,:),[],2))
                set(gca,'ytick',1:size(dat,1),'yticklabel',n(ord))
                set(gca,'xtick',1000:1000:6000,'xticklabel',[0:5])
                colorbar
                axis tight ij
                
                subplot(5,2,(idat-1)*2+2); hold on
                tbs = 1:size(dat,2);
                m1 = mean(datBS); sem1 = std(datBS)./sqrt(size(dat,1));
                m = smoothdata(m1,'movmean',sp.sm); 
                sem = smoothdata(sem1,'movmean',sp.sm);
                rev = m+sem; 
                fwd = m-sem;
                p1 = patch([tbs tbs(length(tbs):-1:1)],...
                    [fwd rev(end:-1:1)],'black');
                p1.FaceAlpha=.2; p1.EdgeAlpha=0;    
                p1.FaceColor = [0 0 0];
                plot(mean(datBS,'omitnan'),'k');
                set(gca,'xtick',1000:1000:6000,'xticklabel',0:5)
            end
        end
        
        set(gcf,'Position',[1958        -307        1003        1633])
        helper_saveandclosefig([thisdir 'Neurons' num2str(ntotal)...
            '_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana}])
    end      
end

close all