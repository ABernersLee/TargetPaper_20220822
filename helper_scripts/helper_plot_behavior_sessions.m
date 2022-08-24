function helper_plot_behavior_sessions(sessionsall,mice,mousenames,dirs,...
    figs_to_run)

typelab = {'probe','NTrepeat'};

if ~isfolder([dirs.figdir '\Behavior\NeuralDataSessions\Accuracy\'])
    mkdir([dirs.figdir '\Behavior\NeuralDataSessions\Accuracy\'])
end

% used to use in Figure 1 - not used anymore in paper
if sum(ismember(figs_to_run,0))>0
    figure; hold on
    for imouse = 1:length(mice)
        sz = sum(sessionsall(:,1)==mice(imouse));
        plot(ones(sz,1)*imouse+rand(sz,1)*.5-.25,...
            sessionsall(sessionsall(:,1)==mice(imouse),3),...
            'o','Color',[0 .9 0])        
        errorbar(imouse,...
            mean(sessionsall(sessionsall(:,1)==mice(imouse),3),'omitnan'), ...
            std(sessionsall(sessionsall(:,1)==mice(imouse),3),'omitnan')...
            ./sqrt(sz),'-','LineWidth',2,'Color',[0 .6 0])
        plot(ones(sz,1)*imouse+rand(sz,1)*.5-.25,...
            1-sessionsall(sessionsall(:,1)==mice(imouse),4)...
            ,'o','Color',[0 0 .9])        
        errorbar(imouse,...
            mean(1-sessionsall(sessionsall(:,1)==mice(imouse),4),'omitnan'),...
            std(1-sessionsall(sessionsall(:,1)==mice(imouse),4),'omitnan')...
            ./sqrt(sz),'-','LineWidth',2,'Color',[0 0 .6])
    end
    set(gca,'xtick',1:4,'xticklabels',mousenames(1:4))
    xlabel('Mouse')
    ylabel('Choice index')
    ylim([-.1 1.1])
    set(gca,'ytick',0:.2:1)
    title('Sessions im using for analysis')
    aa = gca;
    legend([aa.Children(4) aa.Children(2)], ...
        {'Target';'Non-target'},'Location','best')
    set(gcf,'units','centimeters'); 
    set(gcf,'Position',[0 0 12 9])
    helper_saveandclosefig([dirs.figdir ...
        '\Behavior\NeuralDataSessions\Accuracy\AcrossAllSessionsUsing'])
end

% figure 7
if sum(ismember(figs_to_run,[0 7]))>0
    for itype = 1:2
        if (sum(ismember(figs_to_run,7))>0) && (itype == 2)
            continue
        end
        numtt = 3-(itype-1);
        figure; hold on
        mice = unique(sessionsall(:,1));
        if itype==1
            mice(3) = [];
        end
        for imouse = 1:length(mice)
            sz = sum(sessionsall(:,1)==mice(imouse) & ...
                sessionsall(:,8)==itype);

            plot(ones(sz,1)*imouse+.2*median(1:numtt)+rand(sz,1)*.5-.25...
                ,sessionsall(sessionsall(:,1)==mice(imouse) ...
                & sessionsall(:,8)==itype,3),'o','Color',[0 .9 0])   

            errorbar(imouse+.2*median(1:numtt),...
                mean(sessionsall(sessionsall(:,1)==mice(imouse)...
                & sessionsall(:,8)==itype,3),'omitnan'),...
                std(sessionsall(sessionsall(:,1)==mice(imouse) ...
                & sessionsall(:,8)==itype,3),'omitnan')...
                ./sqrt(sz),'-','LineWidth',2,'Color',[0 .6 0])

            for itt = 1:numtt

                plot(ones(sz,1)*imouse+.2*itt+rand(sz,1)*.1,...
                    1-sessionsall(sessionsall(:,1)==mice(imouse)...
                    & sessionsall(:,8)==itype,itt+4),'o','Color',[0 0 .9])        

                errorbar(imouse+.2*itt,...
                    mean(1-sessionsall(sessionsall(:,1)==mice(imouse) & ...
                    sessionsall(:,8)==itype,itt+4),'omitnan'),...
                    std(1-sessionsall(sessionsall(:,1)==mice(imouse) & ...
                    sessionsall(:,8)==itype,itt+4),'omitnan')...
                    ./sqrt(sz),'-','LineWidth',2,'Color',[0 0 .6])
            end
        end
        tcks = (1:length(mice))'+((1:numtt)*.2);
        tcks = sort(tcks(:))';
        tcksl = (ones(length(mice),1)*(1:numtt))';
        set(gca,'xtick',tcks,'xticklabels',num2str(tcksl(:)))
        xlabel('Mouse')
        ylabel('Choice index')
        ylim([-.1 1.1])
        set(gca,'ytick',0:.2:1)
        title([typelab{itype} ', Sessions im using for analysis'])
        
        aa = gca;
        legend([aa.Children(end) aa.Children(2)], ...
            {'Target';typelab{itype}},'Location','best')
        
        set(gcf,'units','centimeters'); set(gcf,'Position',[0 0 12 9])
        helper_saveandclosefig([dirs.figdir ...
            ['\Behavior\NeuralDataSessions\Accuracy\...' ...
            'AcrossAllSessionsUsing_ProbeRepeat_'] typelab{itype}])
    end
end

%%%% Figure 5 

if sum(ismember(figs_to_run,[0 5]))>0

    itype = 2;
    A = [sessionsall(sessionsall(:,8)==itype,3);...
        sessionsall(sessionsall(:,8)==itype,4);...
        mean(sessionsall(sessionsall(:,8)==itype,5:7),2,'omitnan')];
    C1 = ones(sum(sessionsall(:,8)==itype),1);
    C = cat(1,C1,C1*2,C1*3);

    B = repmat(sessionsall(sessionsall(:,8)==itype,1:2),[3 1]);
    [ps,tbl,~] = anovan(A,cat(2,B,C),'continuous',2,...
        'varnames',{'Mice';'Sessions';'Trialtype'},'display','off');
    [r,p] = corr(A,B(:,2),'rows','complete');


    figure; hold on
    pcol1 = [.7 1 .7; .7 .7 1;.7 .7 .7];
    pcol = [0 .9 0; 0 0 .6;0 0 0];
    pl = NaN(3,1);
    for itrialtype = 1:3
        mdat = NaN(size(sessionsall,1),length(mice));
        for imouse = 1:length(mice)
            ind = sessionsall(:,1)==mice(imouse) & sessionsall(:,8)==itype;
            if itrialtype>2 
                dat1 = mean(sessionsall(ind,5:7),2,'omitnan');
            else
                dat1 = sessionsall(ind,2+itrialtype);
            end
            if itrialtype==1
                dat = dat1;
            else
                dat = 1-dat1;
            end
            %dat instead of sessionsall(ind,2) for excluding nonrecording 
            % days
            mdat(sessionsall(ind,2),imouse) = dat; 
            plot(sessionsall(ind,2),dat,'o-','Color',pcol1(itrialtype,:))
        end   
        pl(itrialtype) = errorbar(mean(mdat,2,'omitnan'), ...
            std(mdat,[],2,'omitnan')./sqrt(sum(~isnan(mdat),2,'omitnan')), ...
            '-','Color',pcol(itrialtype,:)); 
    end
    N = length(unique(sessionsall(sessionsall(:,8)==itype,2)+...
        (sessionsall(sessionsall(:,8)==itype,1)*100)));

    axis tight
    xl = get(gca,'xlim');
    xlim([xl(1)-1 xl(2)+1])
    ylim([-.1 1.1])
    set(gca,'xtick',5:5:25,'ytick',0:.2:1)
    xlabel('Session')
    ylabel('Choice index')
    legend(pl,'target','nontarget','nt repeats','Location','best')
    set(gcf,'units','centimeters'); set(gcf,'Position',[0 0 12 9])
    trp = helper_text_rp(r,p,2);
    title({['N = ' num2str(N) ', '  typelab{itype}];...
        [trp ' ANOVAp = ' num2str(ps(2))]})
    helper_saveandclosefig([dirs.figdir ...
        '\Behavior\NeuralDataSessions\Accuracy\OverSessions_' ...
        typelab{itype}])
    writecell(tbl,[dirs.figdir ...
        '\Behavior\NeuralDataSessions\Accuracy\OverSessions_' ...
        typelab{itype} '.xlsx'])
end