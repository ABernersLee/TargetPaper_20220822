function plot_RSA_Averages(corr_repeats,all_labels,dirs,itype)

% plot RSA averages for supplementary figure 7


thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'RSA_Averages/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

if size(corr_repeats,2)==12
    %5-Probe
    nk = nchoosek(1:5,2);
    labs = {'T';'NT';'P1';'P2';'P3'};
elseif size(corr_repeats,2)==8
    %4-Nontarget repeats
    nk = nchoosek(1:4,2);
    labs = {'T';'NT';'NT1';'NT2'};
end

xlab = cell(size(nk,1),1);
figure; hold on
for ii = 1:size(nk,1)
    errorbar(ii-.1,mean(corr_repeats(:,ii)),std(corr_repeats(:,ii))./...
        sqrt(sum(~isnan(corr_repeats(:,ii)))),'k');
    plot(ii+.1,corr_repeats(:,ii),'ok')
    xlab{ii} = [labs{nk(ii,1)} '-' labs{nk(ii,2)}];
end
set(gca,'xtick',1:size(nk,1),'xticklabel',xlab)
plot([0 size(nk,1)+1],[0 0],'r--')
xlim([0 size(nk,1)+1])
ylabel('PV correlation')
helper_saveandclosefig([thisdir all_labels.typelab{itype} '_Average'...
    all_labels.addon])