function plot_RSA_Averages(corr_repeats,all_labels,dirs,itype,prop_selective)

% plot RSA averages for supplementary figure 7


thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'RSA_Averages/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

if size(corr_repeats,2)==13
    %5-Probe
    nk = nchoosek(1:5,2);
    labs = {'T';'NT';'P1';'P2';'P3'};
elseif size(corr_repeats,2)==9
    %4-Nontarget repeats
    nk = nchoosek(1:4,2);
    labs = {'T';'NT';'NT1';'NT2'};
end

xlab = cell(size(nk,1),1);
figure; hold on
for ii = 1:size(nk,1)
%     errorbar(ii-.1,mean(corr_repeats(:,ii)),std(corr_repeats(:,ii))./...
%         sqrt(sum(~isnan(corr_repeats(:,ii)))),'k');
%changed for non-zscored testing
    errorbar(ii-.1,mean(corr_repeats(:,ii),'omitnan'),...
        std(corr_repeats(:,ii),'omitnan')./...
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

if size(corr_repeats,2)==9
    ps_labs = {'TNT','NTrepeat'};
    ps = nan(size(corr_repeats,1),2);
    for ic = 1:size(corr_repeats,1)
        ps(ic,1) = prop_selective.TNT(corr_repeats(ic,8),corr_repeats(ic,7)); %*corr_repeats(ic,9);
        ps(ic,2) = prop_selective.NTrepeat(corr_repeats(ic,8),corr_repeats(ic,7)); %*corr_repeats(ic,9);
    end
    
    for iprop = 1:2
        figure; hold on;
        for ii = 1:size(nk,1)
            subplot(size(nk,1)/2,2,ii)
            hold on
            plot(ps(:,iprop),corr_repeats(:,ii),'ok')
            xlabel(['Proportion ' ps_labs{iprop}])
            ylabel([xlab{ii} ' coeff.'])
            if iprop==1
                [r,p] = corr(ps(:,iprop),corr_repeats(:,ii),"rows","complete");
                title(helper_text_rp(r,p,2))
            else
                [r1,p1] = corr(ps(:,iprop),corr_repeats(:,ii),"rows","complete");
                [r2,p2] = partialcorr(ps(:,iprop),corr_repeats(:,ii),ps(:,1));
                title([helper_text_rp(r1,p1,2) '; ' helper_text_rp(r2,p2,2)])
            end

        end
        suptitle({['Relationship between the proportion of ' ps_labs{iprop} ' selective neurons'],...
            'and the PV correlation coefficients in each session'})
        set(gcf,'Position',[   1198          95         690         869])
        helper_saveandclosefig([thisdir 'RelationshipProp_' ps_labs{iprop} '_'...
            all_labels.typelab{itype} all_labels.addon])
    end
end



figure; hold on

dat1 = corr_repeats(:,1:size(nk,1));
dat2 = repmat(corr_repeats(:,end),[size(nk,1) 1]);

subplot(1,2,1); hold on
plot(dat2(:),dat1(:),'ok')
xlabel('Number of neurons in session')
ylabel('All correlation comparrisons')

dat1 = std(corr_repeats(:,1:size(nk,1)),[],2);
dat2 = corr_repeats(:,end);

subplot(1,2,2); hold on
plot(dat2(:),dat1(:),'ok')
xlabel('Number of neurons in session')
ylabel('Standard deviation of correlation comparisons')
[r,p] = corr(dat2,dat1,'rows','complete');
legend(helper_text_rp(r,p,2))

% set(gcf,'Position',[1389         608        1564         249])

helper_saveandclosefig([thisdir all_labels.typelab{itype} '_Average'...
    all_labels.addon '_num_nuerons'])
