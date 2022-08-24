function TargetPaper_SupplementaryFigure1(dirs)

%%% Accurary for other odor cohorts, with different target mixtures used,
%%% for three other groups of mice
groups(1).dir = 'F:\LizData\behavior\mice_group_first_ID_5and6\';
groups(1).mousenum = [5,6];
groups(1).label = 'A';
groups(2).dir = 'F:\LizData\behavior\mice_groupA\';
groups(2).mousenum = [1,3,4];
groups(2).label = 'B';
groups(3).dir = 'F:\LizData\behavior\mice_groupA\';
groups(3).mousenum = [5,6,7];
groups(3).label = 'C';
%Group 2, mouse 1, days 42-67 has a PID_signal that is the same each day

col = [0 1 0; 0 0 1];
figure; hold on;
set(gcf,'Position',[ 2079         277         428         727])
for igroup = 1:size(groups,2)    
    dat2 = [];
    for imouse = 1:length(groups(igroup).mousenum)
        cd(groups(igroup).dir)
        dd = dir(['*training_mouse_' ...
            num2str(groups(igroup).mousenum(imouse)) '_day_*']);
        
        for ifile = 1:size(dd,1)
            load([dd(ifile).name],'day',...
                'Is_odor_target_or_non_across_all_trials',...
                'times_of_correct_licks_to_port1_across_trials',...
                'times_of_correct_licks_to_port2_across_trials'...
                ,'phase_of_training',...
                'total_correct_licks_NONtarget_across_all_trials',...
                'total_correct_licks_TARGET_across_all_trials',...
                'total_Incorrect_licks_NONtarget_across_all_trials',...
                'total_Incorrect_licks_TARGET_across_all_trials')
            if exist('phase_of_training','var')
                if phase_of_training==1
                    continue
                end
            end
            isTNT = Is_odor_target_or_non_across_all_trials;
            NT_C = total_correct_licks_NONtarget_across_all_trials;
            T_C = total_correct_licks_TARGET_across_all_trials; 
            NT_I = total_Incorrect_licks_NONtarget_across_all_trials;
            T_I = total_Incorrect_licks_TARGET_across_all_trials;

            %%% dont count trailing no lick trials
            nolicks = NT_C==0 & T_C==0 & NT_I==0 & T_I==0;
            cutoff = find(nolicks==0,1,'last');
            isTNT(cutoff+1:end) = NaN;
            T = length(times_of_correct_licks_to_port1_across_trials)...
                ./sum(isTNT==1);
            NT = length(times_of_correct_licks_to_port2_across_trials)...
                ./sum(isTNT==0);

            dat2 = cat(1,dat2,[imouse day phase_of_training T 1-NT]);
            clear('day',...
                'Is_odor_target_or_non_across_all_trials',...
                'times_of_correct_licks_to_port1_across_trials',...
                'times_of_correct_licks_to_port2_across_trials'...
                ,'phase_of_training',...
                'total_correct_licks_NONtarget_across_all_trials',...
                'total_correct_licks_TARGET_across_all_trials',...
                'total_Incorrect_licks_NONtarget_across_all_trials',...
                'total_Incorrect_licks_TARGET_across_all_trials')
        end
    end

    %sort for days not including the freedrop
    dat = NaN(100,length(groups(igroup).mousenum),2);
    for imouse = 1:size(dat,2)
        d = dat2(dat2(:,1)==imouse,:);
        [~,d2] = sort(d(:,2));
        d3 = d(d2,:);
        dat(1:length(d2),imouse,1:2) = d3(:,4:5);
    end

    subplot(size(groups,2),1,igroup);
    hold on;
    clear pp
    for ii = 1:2
        m = squeeze(mean(dat(:,:,ii),2,'omitnan'))';
        tbs = 1:40; 
        m = m(tbs);
        sem = squeeze(std(dat(:,:,ii),[],2,'omitnan')./...
            sqrt(sum(~isnan(dat(:,:,ii)),2)))';
        sem = sem(tbs);
        rev = m+sem; fwd = m-sem;
        p1 = patch([tbs tbs(length(tbs):-1:1)],...
            [fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = col(ii,:);
        pp(ii) = plot(m,'color',col(ii,:));
    end

    ylabel('Choice index')
    xlabel('Session')
    legend(pp,'Target','Nontarget','Location','Best')
    title(['Group ' groups(igroup).label ', N = ' ...
        num2str(length(groups(igroup).mousenum))])
end
helper_saveandclosefig([dirs.figdir '\Behavior\TrainingIncluded\...' ...
    'SupplementaryFigure1_OtherOdorCohorts'])