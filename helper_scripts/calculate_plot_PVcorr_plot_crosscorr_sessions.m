function calculate_plot_PVcorr_plot_crosscorr_sessions(PVall2, ...
    trialdatsession,figdir,mouse,ses,tl,all_labels)

thisdir = [figdir '/PopulationVectorsCorrelation/' ...
        'PV_Sessions_Imagesc/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

figlab = {'Target','NT',[tl ' 1'],[tl ' 2'],[tl ' 3']};
[item1,order1] = sort(trialdatsession(:,1));
pinds = find(diff(item1)>0);
pinds1 = pinds;
PVord = PVall2(order1,order1);
PVord2 = PVord;


figure; 
imagesc(PVord); 
axis xy; hold on; 
set(gcf,'Position',[1877         231        1104         780])
plot([pinds1 pinds1]',(ones(size(pinds1,1),1)*[1 size(PVall2,1)])'...
    ,'w','LineWidth',3) 
plot((ones(size(pinds1,1),1)*[1 size(PVall2,1)])',[pinds1 pinds1]', ...
    'w','LineWidth',3)

nk = [nchoosek(1:length(pinds)+1,2); (ones(2,1)*(1:length(pinds)+1))'];
pinds = [0;pinds;size(PVord,2)];

for ii = 1:length(nk)        
    ind1 = pinds(nk(ii,1))+1;
    ind2 = pinds(nk(ii,2)+1);    
    tp = mean(mean(PVord(ind1:ind2,ind1:ind2),2,'omitnan'),1,'omitnan');
    tt = text(mean([pinds(nk(ii,1)) pinds(nk(ii,1)+1)])-5, ...
        mean([pinds(nk(ii,2)) pinds(nk(ii,2)+1)])+2, ...
        num2str(round(tp,1,'significant')));
    tt.Color = 'k'; tt.FontSize = 14; tt.FontWeight = 'bold';
    PVord2(ind1:ind2,ind1:ind2) = tp;
end
tt = title(['Mouse ' mouse ', Session ' num2str(ses)]);
tt.FontSize = 18;
set(gca,'xtick',pinds(2:end)-diff(pinds)/2,'ytick',pinds(2:end)-diff(pinds)/2)
set(gca,'xticklabel',figlab(1:length(pinds)-1),'yticklabel',figlab(1:length(pinds)-1))
helper_saveandclosefig([thisdir 'CrossCorr_Sessions_' tl '_' ...
    mouse '_' num2str(ses) '_all' all_labels.addon])

figure; hold on
imagesc(PVord2); 
axis xy; hold on; 
% plot([pinds1 pinds1],[1 size(PVall2,1)],'w','LineWidth',3) 
% plot([1 size(PVall2,1)],[pinds1 pinds1],'w','LineWidth',3)
plot([pinds1 pinds1]',(ones(size(pinds1,1),1)*[1 size(PVall2,1)])', ...
    'w','LineWidth',3) 
plot((ones(size(pinds1,1),1)*[1 size(PVall2,1)])',[pinds1 pinds1]', ...
    'w','LineWidth',3)

set(gcf,'Position',[1877         231        1104         780])
for ii = 1:length(nk)        
    ind1 = pinds(nk(ii,1))+1;
    ind2 = pinds(nk(ii,2)+1);    
    tp = mean(mean(PVord(ind1:ind2,ind1:ind2),2,'omitnan'),1,'omitnan');
    tt = text(mean([pinds(nk(ii,1)) pinds(nk(ii,1)+1)])-5, ...
        mean([pinds(nk(ii,2)) pinds(nk(ii,2)+1)])+2, ...
        num2str(round(tp,1,'significant')));
    tt.Color = 'k'; tt.FontSize = 14; tt.FontWeight = 'bold';
    PVord2(ind1:ind2,ind1:ind2) = tp;
end
axis tight
tt = title(['Mouse ' mouse ', Session ' num2str(ses)]);
tt.FontSize = 18;
set(gca,'xtick',pinds(2:end)-diff(pinds)/2,'ytick',pinds(2:end)-diff(pinds)/2)
set(gca,'xticklabel',figlab(1:length(pinds)-1),'yticklabel',figlab(1:length(pinds)-1))
helper_saveandclosefig([thisdir 'CrossCorr_Sessions_' tl '_' ...
    mouse '_' num2str(ses) '_summary' all_labels.addon])