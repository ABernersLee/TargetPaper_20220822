function helper_plot_lick_lat(lickall,dirs)


%from start of decision period
lick_cutoff_ms = 500;

% Target and nontarget
lickall2 = lickall(:,1); 
%probe
lickall2(lickall(:,1)>100) = 3; 
%NT repeat
lickall2(lickall(:,1)>2 & lickall(:,1)<100) = 4; 

if ~isfolder([dirs.figdir '\Behavior\NeuralDataSessions\Licks\'])
    mkdir([dirs.figdir '\Behavior\NeuralDataSessions\Licks\'])
end

%%% lickall = trialtype, latency of first lick, mouse, accuracy, session, 
%%% trial #, type of session   

lickacc = cat(2, lickall(:,3), lickall(:,5), lickall(:,4), lickall(:,2),...
    lickall2(:,1), lickall(:,6:7));
lickacc(lickacc(:,4)<=lick_cutoff_ms,:) = [];
tlabs = {'Target';'NT';'Other'};


%% by trial type

clump = 8; %venki asked to clump together the later latencies

binsize = 100; stnd = [500 2200];
ind = stnd(1):binsize:stnd(2);
Ms = NaN(length(ind),4); SD = Ms; hs = NaN(length(ind)-1,4);
for itt = 1:3
    lickacc2 = lickacc(lickacc(:,5)==itt,:);
    [h,~,c] = histcounts(lickacc2(:,4),ind);
    c(c>(clump-1)) = clump;
    hs(:,itt) = h;
    M = splitapply(@mean,lickacc2(c>0,3),c(c>0));
    Ms(unique(c(c>0)),itt) = M;
    S = splitapply(@std,lickacc2(c>0,3),c(c>0));
    D = splitapply(@sum,ones(sum(c>0),1),c(c>0));
    SD(unique(c(c>0)),itt) = S./sqrt(D);
end


vc = [0 .6 0;0 0 .6; .6 0 0;.6 0 .6];

figure; hold on;
set(gcf,'Position',[87         135        1621         786])

subplot(2,3,[2:3 5:6]); hold on
yyaxis left
for itt = 1:3
    bar(hs(:,itt)./sum(hs(:,itt)),'FaceColor',vc(itt,:),...
        'EdgeColor','w','FaceAlpha',.3)
end
yl = get(gca,'ylim');
set(gca,'ylim',[yl(1) yl(2)+(range(yl)*.5)])
ylabel('Probability')
yyaxis right
for itt = 1:3
    tbs = 1:length(ind);
    rev = (Ms(:,itt)+SD(:,itt))'; fwd = (Ms(:,itt)-SD(:,itt))';
    tbs(isnan(rev)) = []; fwd(isnan(rev)) = []; rev(isnan(rev)) = [];  
    tbs = [tbs(1:end-1)  (tbs(end)+length(ind)-1)/2]; 
    p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
    p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = vc(itt,:);
    plot(tbs,Ms(~isnan(Ms(:,itt)),itt),'.-',...
        'Color',vc(itt,:),'MarkerSize',20)
end
ylabel('Proportion correct')
set(gca,'xtick',1:5:length(ind),'xticklabel',ind(1:5:length(ind)))
xlabel('Lick latency (ms)')
legend('Target','NT','Probe','Location','Best')


for it = 1:2
    subplot(2,3,1+(it-1)*3); hold on
    ind = lickacc(:,5)==it;
    
    histogram(lickacc(lickacc(:,3)==1 & ind,4),'FaceColor',vc(it,:),...
        'EdgeColor','w','FaceAlpha',.3)
    histogram(lickacc(lickacc(:,3)==0 & ind,4),'FaceColor','k',...
        'EdgeColor','w','FaceAlpha',.3)
    
    xlabel('Latency (ms)')
    ylabel('Count')    
    yl2 = get(gca,'ylim');

    plot([median(lickacc(lickacc(:,3)==1 & ind,4),'omitnan') ...
        median(lickacc(lickacc(:,3)==1 & ind,4),'omitnan')],yl2,...
        '-','color',vc(it,:))
    plot([median(lickacc(lickacc(:,3)==0 & ind,4),'omitnan') ...
        median(lickacc(lickacc(:,3)==0 & ind,4),'omitnan')],yl2,...
        '-','color','k')

    legend('Correct','Incorrect','Location','Best')
    title(tlabs{it})
    
    p = anovan(lickacc(ind,4),lickacc(ind,[3 1:2 6]),'Display','off',...
        'varnames',{'Accuracy','Mouse','Day','Trial#'},'continuous',3:4);    
    text(1500,yl2(2)-range(yl2)/2,['p = ' ...
        num2str(round(p(1),2,'significant'))])
end

helper_saveandclosefig([dirs.figdir ...
    '\Behavior\NeuralDataSessions\Licks\SpeedAccuracyTradeoff'])
 