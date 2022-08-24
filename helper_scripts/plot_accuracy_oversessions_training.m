function plot_accuracy_oversessions_training(dirs,neuron_info,figs_to_run)

%%%% all sessions (not just ones that there are neural recordings in)

plotmice = false;

if ~isfolder([dirs.figdir '\Behavior\TrainingIncluded\IndividualMice\'])
    mkdir([dirs.figdir '\Behavior\TrainingIncluded\IndividualMice\'])
end

lick_cutoff_ms = 500;
mmlab = {'Mean';'Median'};
stagelab = {'Learning';'Over-learning';'Probe'};
stagelab2 = {'Stage 1';'Stage 2';'Stage 3';'Over-learning';'Probe'};
mousenames = {'Q','R','S','T','V'};

training_dirs = cell(5,1);
training_dirs{1} = 'F:\LizData\Behavior_with_tetrodes_Feb2018\Spike_sorted_by_julien\behavior_with_tetrodes_save_all_licks_MouseQ';
training_dirs{2} = 'F:\LizData\Behavior_with_tetrodes_Feb2018\Spike_sorted_by_julien\behavior_with_tetrodes_save_all_licks_MouseR';
training_dirs{3} = 'F:\LizData\Behavior_with_tetrodes_Feb2018\behavior_with_tetrodes_save_all_licks_MouseS\';
training_dirs{4} = 'F:\LizData\Behavior_with_tetrodes_Feb2018\behavior_with_tetrodes_save_all_licks_MouseT\';
training_dirs{5} = 'F:\LizData\Behavior_with_tetrodes_Feb2018\behavior_with_tetrodes_save_all_licks_MouseV\';

vc = [0 .6 0;0 0 .9; .8 0 .8];
exlnan_lab = {'ExclMissNaNs';'IncludeMissNans_exlTrailing'};
exlStg1_lab = {'ExclStage1';'IncludeStage1'};

% exlnan: including/excluding stage one which is just random
% istg: excluding trials where the mouse doesnt lick, or excluding only
%       those trials that are at the end (trailing trials)
for exlnan = 1:2 
    for istg = 1:2

        if exlnan==1 && istg==2
            continue
        end
    
        alldat = NaN(50,size(mousenames,2),5); 
        Tdat = NaN(50,size(mousenames,2),3);     
        alldat_AllProbes = alldat;
        alltrial_AllProbes = NaN(1000,30,size(mousenames,2),5);
        alllick = NaN(50,size(mousenames,2),5,2); alllickc =  alllick; 
        alldatc = alldat; datind = NaN(50,size(mousenames,2),5,6);
        datind_lick = NaN(50,size(mousenames,2),5,2,3);        
        
        for imouse = 1:size(mousenames,2)
    
            
        
            daystotal = [];
    
            trainingdays = [];        
            if exlnan==2
        
                for iday = 1:3
                    if isfile([dirs.dathome '\April\Behavior_training_only\training_mouse_'...
                        num2str(imouse) '_day_' num2str(iday) '.mat'])
                        load([dirs.dathome '\April\Behavior_training_only\training_mouse_'...
                            num2str(imouse) '_day_' num2str(iday) '.mat'],...
                            'Is_odor_target_or_non_across_all_trials',...
                            'times_of_correct_licks_to_port1_across_trials',...
                            'times_of_correct_licks_to_port2_across_trials'...
                            ,'phase_of_training',...
                            'total_correct_licks_NONtarget_across_all_trials',...
                            'total_correct_licks_TARGET_across_all_trials',...
                            'total_Incorrect_licks_NONtarget_across_all_trials',...
                            'total_Incorrect_licks_TARGET_across_all_trials')
    
                        disp(['mouse ' num2str(imouse) ...
                            ' phase_of_training ' ...
                            num2str(phase_of_training)])

                        if phase_of_training==1 && istg==1
                            continue
                        end

                        daystotal = cat(1,daystotal,iday);
                        isTNT = Is_odor_target_or_non_across_all_trials;


                        NT_C = total_correct_licks_NONtarget_across_all_trials;
                        T_C = total_correct_licks_TARGET_across_all_trials; 
                        NT_I = total_Incorrect_licks_NONtarget_across_all_trials;
                        T_I = total_Incorrect_licks_TARGET_across_all_trials;

                        %%% dont count trailing no lick trials
                        nolicks = NT_C==0 & T_C==0 & NT_I==0 & T_I==0;
                        cutoff = find(nolicks==0,1,'last');
                        isTNT(cutoff+1:end) = NaN;
                        T = length(...
                            times_of_correct_licks_to_port1_across_trials)...
                            ./sum(isTNT==1);
                        NT = length...
                            (times_of_correct_licks_to_port2_across_trials)...
                            ./sum(isTNT==0);

                        alldat(iday,imouse,1) = T;
                        alldat(iday,imouse,2) = 1-NT;
                        trainingdays = cat(1,trainingdays,...
                            [iday phase_of_training]);        
                    end
                end
                disp(['mouse ' num2str(imouse)])
    
                cd(training_dirs{imouse})
                d2 = dir("*Day*");
                d2(cell2mat(extractfield(d2,'isdir'))==0) = [];
                for itday = 1:size(d2,1)
                    dayn = str2double(d2(itday).name(4:end));
                    filename_mouse_day = [d2(itday).name '/training_mouse_' ...
                        num2str(imouse) '_day_' num2str(dayn) '.mat'];
                    if (imouse == 4 && (dayn == 5 || dayn == 7)) ||  (imouse == 5 && dayn == 2) 
                        filename_mouse_day = [d2(itday).name '/training_mouse_' ...
                            num2str(1) '_day_' num2str(dayn) '.mat'];
                    elseif  (imouse == 4) && (dayn == 9)
                        filename_mouse_day = [d2(itday).name '/training_mouse_' ...
                            num2str(5) '_day_' num2str(dayn) '.mat'];
                    elseif  (imouse == 5) && (dayn == 3)
                        filename_mouse_day = [d2(itday).name '/training_mouse_' ...
                            num2str(imouse) '_day_' num2str(6) '.mat'];
                    elseif  (imouse == 5) && (dayn == 8)
                        filename_mouse_day = [d2(itday).name '/training_mouse_' ...
                            num2str(imouse) '_day_' num2str(10) '.mat'];
                    end
                        
                    load(filename_mouse_day,...
                                'Is_odor_target_or_non_across_all_trials',...
                                'times_of_correct_licks_to_port1_across_trials',...
                                'times_of_correct_licks_to_port2_across_trials'...
                                ,'phase_of_training')
                    isTNT = Is_odor_target_or_non_across_all_trials;
                    T = length(times_of_correct_licks_to_port1_across_trials)./sum(isTNT==1);
                    NT = length(times_of_correct_licks_to_port2_across_trials)./sum(isTNT==0);
                    Tdat(dayn,imouse,1) = T;
                    Tdat(dayn,imouse,2) = 1-NT;
                    Tdat(dayn,imouse,3) = phase_of_training;
                    trainingdays = cat(1,trainingdays,[dayn phase_of_training]);
                end
            end
        
            %training data
            load([dirs.dathome 'data_structure_for_each_mouse\mouse_' ...
            mousenames{imouse} '_behavior.mat'],...
            'behavior_struct','probe_behavior_struct')
            dat = behavior_struct; dat2 = probe_behavior_struct;
            for itrialtype = 1:size(dat,2)
                days = unique(extractfield(dat(itrialtype).performance,...
                    'day_of_training'));   
                if exlnan==1
                    choiceinx = extractfield(dat(itrialtype).performance,...
                        'avg_session_choice_index');
                elseif exlnan==2

                     %%%% Changed to get rid of trailing no-lick trials
                    choiceinx = arrayfun(@(x) ...
                        mean(dat(itrialtype).performance...
                        (x).within_session_choice_index_every_trial_no_moving_avg...
                        (1:find(~isnan(dat(itrialtype).performance...
                        (x).within_session_choice_index_every_trial_no_moving_avg)...
                        ,1,'last'))==1)...
                        ,1:length(dat(itrialtype).performance));
                end
                lickL = NaN(length(days),2);
                for iday = 1:length(days)
                    theselicks = ...
                        extractfield(dat(itrialtype).lick_latencies(iday),...
                        'lick_latencies_correct_trials');
                    theselicks(theselicks<lick_cutoff_ms) = NaN;
                    theselicks = theselicks-lick_cutoff_ms;
                    lickL(iday,:) = [mean(theselicks,'omitnan') ...
                        median(theselicks,'omitnan')];
                end
                
                % these do not include trials where the mouse didn't lick at all
                % (excludes the NaNs) - just like for figure 4
                alldat(days,imouse,itrialtype) = choiceinx;
                alllick(days,imouse,itrialtype,:) = lickL;
                if itrialtype==1
                    daystotal = cat(1,daystotal,days');
                end
            end
            days2total = [];
                
            %probe data
            for itrialtype = 1:size(dat2,2)
                days2 = unique(extractfield(dat2(itrialtype).performance,...
                    'day_of_training'));
                if exlnan==1
                    choiceinx = extractfield(dat2(itrialtype).performance,...
                        'avg_session_choice_index');
                elseif exlnan==2

                     %%%% Changed to get rid of trailing no-lick trials
                     choiceinx = arrayfun(@(x) ...
                        mean(dat2(itrialtype).performance...
                        (x).within_session_choice_index_every_trial_no_moving_avg...
                        (1:find(~isnan(dat2(itrialtype).performance...
                        (x).within_session_choice_index_every_trial_no_moving_avg)...
                        ,1,'last'))==1)...
                        ,1:length(dat2(itrialtype).performance));
                end
                lickL = NaN(length(days2),2);
                for iday = 1:length(days2)
                    theselicks = ...
                        extractfield(dat2(itrialtype).lick_latencies(iday),...
                        'lick_latencies_correct_trials');
                    theselicks(theselicks<lick_cutoff_ms) = NaN;
                    theselicks = theselicks-lick_cutoff_ms;
                    lickL(iday,:) = [mean(theselicks,'omitnan') ...
                        median(theselicks,'omitnan')];
                    
                    choiceinxtrial = ...
                        extractfield(dat2(itrialtype).performance(iday),...
                        'within_session_choice_index_every_trial_no_moving_avg');
                    alltrial_AllProbes(1:sum(~isnan(choiceinxtrial)),...
                        iday,imouse,itrialtype) = ...
                        choiceinxtrial(~isnan(choiceinxtrial));     
                end
        
                % these do not include trials where the mouse didn't lick at all
                % (excludes the NaNs) - just like for figure 4
                alldat_AllProbes(days2,imouse,itrialtype) = choiceinx;
                alldat(days2,imouse,itrialtype) = choiceinx;
                alllick(days2,imouse,itrialtype,:) = lickL;
                days2total = cat(1,days2total,days2);
            end
        
            % Non-target repeat overlearning % phase_of_training 4s
            NTR = unique(neuron_info(neuron_info(:,3+2)>2 & neuron_info(:,2)==imouse,3),'rows');
        
            % Probe days
            P = unique(neuron_info(neuron_info(:,3+1)>2 & neuron_info(:,2)==imouse,3),'rows');    
        
            if exlnan==2
                NTR2 = sort(trainingdays(trainingdays(:,2)==4,1));
                PH1 = sort(trainingdays(trainingdays(:,2)==1,1));
                PH2 = sort(trainingdays(trainingdays(:,2)==2,1));
                PH3 = sort(trainingdays(trainingdays(:,2)==3,1));
                if isempty(NTR)
                    NTR = NTR2(~ismember(NTR2,[PH1; PH2; PH3]));
                end
            end
    
            % Training days (first ones) 
                % phase_of_training 2s and 3s
            d = [daystotal' days2total(1,:)];    
            OTH = d(~ismember(d,[NTR; P])); 
            OTH2 = OTH(OTH<min([P; NTR])); 
            
            if istg==1 && exlnan==2
                OTH2(ismember(OTH2,PH1)) = [];
                alldat(PH1,imouse,:) = NaN;
                alllick(PH1,imouse,:) = NaN;
            end

            if exlnan==2
                if sum(ismember(trainingdays(:,1),NTR) & trainingdays(:,2)<4)>0
                    error('Using a NTR day that has not gotten to criterion')
                end
                if sum(ismember(trainingdays(:,1),P) & trainingdays(:,2)==4)>0
                    disp('Probe may be mistaken for NTR')
                end
            end
            
            %combine the probes together (could change this)
            alldatc(:,imouse,:) = alldat(:,imouse,:);
            alldatc(:,imouse,3) = mean(alldatc(:,imouse,3:end),3);
            alllickc(:,imouse,:,:) = alllick(:,imouse,:,:);
            alllickc(:,imouse,3,:) = mean(alllickc(:,imouse,3:end,:),3,'omitnan');
        
            if plotmice
                %%% plot individual mice
                figure; hold on
                for itrialtype = 1:3
                    plot(OTH,alldatc(OTH,imouse,itrialtype),'.','color',vc(itrialtype,:))
                end
                for itrialtype = 1:3 %max([size(dat,2) size(dat2,2)])
                    plot(NTR,alldatc(NTR,imouse,itrialtype),'o-','color',vc(itrialtype,:))
                    plot(P,alldatc(P,imouse,itrialtype),'+-','color',vc(itrialtype,:))
                end
                 ylim([-.1 1.1])
                xl = get(gca,'xlim');
                xlim([0 xl(2)+1])
                set(gca,'xtick',1:xl(2),'ytick',0:.2:1)
                xlabel('Session')
                ylabel('Choice index')
                legend('target','nontarget','probe','Location','best')
                title(['Mouse ' mousenames{imouse}])
                helper_saveandclosefig([dirs.figdir ...
                    '\Behavior\TrainingIncluded\IndividualMice\AcrossDaysAllTypes_Mouse' ...
                    mousenames{imouse}])
            
                for imm = 1:2
                    figure; hold on
                    for itrialtype = 1:3
                        plot(OTH,alllickc(OTH,imouse,itrialtype,imm),'.','color', ...
                            vc(itrialtype,:))
                    end
                    for itrialtype = 1:3
                        plot(NTR,alllickc(NTR,imouse,itrialtype,imm),'o-','color', ...
                            vc(itrialtype,:))
                        plot(P,alllickc(P,imouse,itrialtype,imm),'+-','color', ...
                            vc(itrialtype,:))
                    end
                    xl = get(gca,'xlim');
                    xlim([0 xl(2)+1])
                    set(gca,'xtick',1:xl(2))
                    xlabel('Session')
                    ylabel([mmlab{imm} ' Lick Latency'])
                    legend('target','nontarget','probe','Location','best')
                    title(['Mouse ' mousenames{imouse}])
                    helper_saveandclosefig([dirs.figdir ...
                        '\Behavior\TrainingIncluded\IndividualMice\AcrossDaysAllTypes_Mouse'...
                        mousenames{imouse} '_Lick' mmlab{imm}])
                end
            end
        
        
        
            %fill in group data
            datind(1:length(OTH2),imouse,1:3,1) = alldatc(OTH2,imouse,1:3);
            datind(1:length(NTR),imouse,1:3,2) = alldatc(NTR,imouse,1:3);
            datind(1:length(P),imouse,1:3,3) = alldatc(P,imouse,1:3);
    
            if exlnan==2
                datind(1:length(PH1),imouse,1:3,4) = alldatc(PH1,imouse,1:3);
                datind(1:length(PH2),imouse,1:3,5) = alldatc(PH2,imouse,1:3);
                datind(1:length(PH3),imouse,1:3,6) = alldatc(PH3,imouse,1:3);
            end
            
            datind_lick(1:length(OTH2),imouse,1:3,:,1) = alllickc(OTH2,imouse,1:3,:);
            datind_lick(1:length(NTR),imouse,1:3,:,2) = alllickc(NTR,imouse,1:3,:);
            datind_lick(1:length(P),imouse,1:3,:,3) = alllickc(P,imouse,1:3,:);                    
        end
        
        
        datind = datind(:,:,1:3,:);
    
        datind_lick = datind_lick(:,1:4,1:3,:,:);
        
        % figure 1
        if sum(ismember(figs_to_run,[0 1]))>0
            figure; hold on
            for istage = 1:3
                subplot(1,3,istage); hold on
                for itrialtype = 1:3
                    plot(1:size(datind,1),...
                        mean(datind(:,:,itrialtype,istage),2,'omitnan'),'-',...
                        'color',vc(itrialtype,:))
                end
            
                for itrialtype = 1:3
                    erbar = sum(~isnan(datind(:,:,itrialtype,istage)),2)>1;
                    erdot = sum(~isnan(datind(:,:,itrialtype,istage)),2)==1;
                    errorbar(find(erbar),mean(datind(erbar,:,itrialtype,istage),2,'omitnan')...
                        ,std(datind(erbar,:,itrialtype,istage),[],2,'omitnan')...
                        ./sqrt(sum(~isnan(datind(erbar,:,itrialtype,istage)),2)),...
                        'color',vc(itrialtype,:),'LineStyle','none')
                    plot(find(erdot),datind(erdot,:,itrialtype,istage),'.',...
                        'color',vc(itrialtype,:),'MarkerSize',15)
                end
                
                
                axis tight
                ylim([-.1 1.1])
                xl = get(gca,'xlim');
                xlim([0 xl(2)+1])
                set(gca,'xtick',1:xl(2),'ytick',0:.2:1)
                xlabel('Session')
                ylabel('Choice index')
                if istage>2
                    legend('target','nontarget','probe','Location','best')
                else
                    legend('target','nontarget','Location','best')
                end
                title(stagelab{istage}) 
            end
            set(gcf,'units','centimeters'); set(gcf,'Position',[0 0 36 9])
            helper_saveandclosefig([dirs.figdir...
                '\Behavior\TrainingIncluded\AcrossDaysAllTypes_AllMice_addedEarly_'...
                exlnan_lab{exlnan} '_' exlStg1_lab{istg}])
    
            A = cat(3,datind(:,:,1,1:3),1-datind(:,:,2,1:3)); 
            B = ones(size(A)); for i = 1:size(A,1); B(i,:,:,:) = i; end
            C = ones(size(A)); for i = 1:size(A,2); C(:,i,:,:) = i; end
            D = ones(size(A)); for i = 1:size(A,3); D(:,:,i,:) = i; end
            E = ones(size(A)); for i = 1:size(A,4); E(:,:,:,i) = i; end
            [~,tbl] = anovan(A(:),cat(2,B(:),C(:),D(:),E(:)),'varnames',...
                {'Session';'Mice';'T/NT';'Stage'},'display','off',...
                'continuous',1,'model',[eye(4,4);[1 0 0 1]]);
            writecell(tbl,[dirs.figdir...
                '\Behavior\TrainingIncluded\AcrossDaysAllTypes_AllMice_addedEarly_'...
                exlnan_lab{exlnan}  '_' exlStg1_lab{istg} '.xlsx'])
    
            for istage =  1:3
                A1 = A(E==istage); 
                B1 = B(E==istage); 
                C1 = C(E==istage); 
                D1 = D(E==istage); 
                [~,tbl] = anovan(A1(:),cat(2,B1(:),C1(:),D1(:)),...
                    'varnames',{'Session';'Mice';'T/NT'},'display','off',...
                    'continuous',1);
                writecell(tbl,[dirs.figdir...
                    '\Behavior\TrainingIncluded\AcrossDaysAllTypes_AllMice_addedEarly_'...
                    exlnan_lab{exlnan} '_Stage' num2str(istage) '_' ...
                    stagelab{istage}  '_' exlStg1_lab{istg} '.xlsx'])
            end
            
        end
    
         % figure 1 ALT - not in paper
        if sum(ismember(figs_to_run,[0 1]))>0 && exlnan==2
            figure; hold on
            numst = [4 5 6 2 3];
            for istage2 = 1:5
                istage = numst(istage2);
                subplot(1,5,istage2); hold on
                for itrialtype = 1:3
                    plot(1:size(datind,1),...
                        mean(datind(:,:,itrialtype,istage),2,'omitnan'),'-',...
                        'color',vc(itrialtype,:))
                end
            
                for itrialtype = 1:3
                    erbar = sum(~isnan(datind(:,:,itrialtype,istage)),2)>1;
                    erdot = sum(~isnan(datind(:,:,itrialtype,istage)),2)==1;
                    errorbar(find(erbar),mean(datind(erbar,:,itrialtype,istage),2,'omitnan')...
                        ,std(datind(erbar,:,itrialtype,istage),[],2,'omitnan')...
                        ./sqrt(sum(~isnan(datind(erbar,:,itrialtype,istage)),2)),...
                        'color',vc(itrialtype,:),'LineStyle','none')
                    plot(find(erdot),datind(erdot,:,itrialtype,istage),'.',...
                        'color',vc(itrialtype,:),'MarkerSize',15)
                end
                
                
                axis tight
                ylim([-.1 1.1])
                xl = get(gca,'xlim');
                xlim([0 xl(2)+1])
                set(gca,'xtick',1:xl(2),'ytick',0:.2:1)
                xlabel('Session')
                ylabel('Choice index')
                if istage>2
                    legend('target','nontarget','probe','Location','best')
                else
                    legend('target','nontarget','Location','best')
                end
                title(stagelab2{istage2}) 
            end
            set(gcf,'units','centimeters'); set(gcf,'Position',[0 0 36 9])
            helper_saveandclosefig([dirs.figdir...
                '\Behavior\TrainingIncluded\AcrossDaysAllTypes_AllMice_TwoPhases_' ...
                exlnan_lab{exlnan} '_' exlStg1_lab{istg}])
        end
                
        % licks - not in paper
        if sum(ismember(figs_to_run,0))>0
            for imm = 1:2
                figure; hold on
                for istage = 1:3
                    subplot(1,3,istage); hold on
                    for itrialtype = 1:3
                        errorbar(mean(datind_lick(:,:,itrialtype,imm,istage),2,'omitnan') ...
                            ,std(datind_lick(:,:,itrialtype,imm,istage),[],2,'omitnan') ...
                            ./sqrt(sum(~isnan(datind_lick(:,:,itrialtype,imm,istage)),2)) ...
                            ,'color',vc(itrialtype,:))
                    end
                    axis tight
                    ylim([100 1100])
                    xl = get(gca,'xlim');
                    xlim([0 xl(2)+1])
                    set(gca,'xtick',1:xl(2))
                    xlabel('Session')
                    ylabel([mmlab{imm} ' Licks (ms since rew. avail.'])
                    xx = repmat((1:size(datind_lick,1))',[1 size(datind_lick,2)]);
                    xx2 = repmat(1:size(datind_lick,2),[size(datind_lick,1) 1]);
                    yy = mean(datind_lick(:,:,1:2,imm,istage),3,'omitnan');
                    [r,p] = corr(xx(:),yy(:),'rows','complete');
                    p2 = anovan(yy(:),[xx(:) xx2(:)],'varnames',{'Session';'Mouse'},...
                        'continuous',1,'display','off');
                    if istage>2
                        legend('target','nontarget','probe','Location','best')
                    else
                        legend('target','nontarget','Location','best')
                    end
                    title([stagelab{istage} ' ' helper_text_rp(r,p,2) ...
                        ' Ap = ' num2str(round(p2(1),2,'significant'))]) 
                end
                set(gcf,'Position',[ 680         414        1092         564])
                helper_saveandclosefig([dirs.figdir... 
                    '\Behavior\TrainingIncluded\AcrossDaysAllTypes_AllMice_Licks_' ...
                    mmlab{imm} '_' exlnan_lab{exlnan} '_' exlStg1_lab{istg}])
            end
        end
        
        %figure 7
        % across sessions probes
        if sum(ismember(figs_to_run,[0 7]))>0
            figure; hold on
            for itrialtype = 1:size(alldat_AllProbes,3)
                  toplotind = sum(~isnan(...
                      alldat_AllProbes(:,1:4,itrialtype)),2)>1;
                  errorbar(mean(alldat_AllProbes(toplotind,1:4,itrialtype)...
                      ,2,'omitnan')...
                      ,std(alldat_AllProbes(toplotind,1:4,itrialtype),[],...
                      2,'omitnan')...
                      ./sqrt(sum(~isnan(alldat(toplotind,1:4,itrialtype)),...
                      2)),'-','LineWidth',2)
            end
            ylim([-.1 1.1])
            xl = get(gca,'xlim');
            xlim([0 xl(2)+1])
            set(gca,'xtick',1:xl(2),'ytick',0:.2:1)
            xlabel('Day of trianing')
            ylabel('Choice index')
            legend('target','nontarget','probe 1','prboe 2','probe 3',...
                'Location','best')
            title('Only the 4 mouse used in data analysis, all sessions')
            
            set(gcf,'units','centimeters'); set(gcf,'Position',[0 0 12 9])
            helper_saveandclosefig([dirs.figdir ...
                '\Behavior\TrainingIncluded\AcrossDays_probe_FourNeuralMice_'...
                exlnan_lab{exlnan}])
    
            A = alldat_AllProbes(:,:,3:5);
            B = ones(size(A)); for i = 1:size(A,1); B(i,:,:) = i; end
            C = ones(size(A)); for i = 1:size(A,2); C(:,i,:) = i; end
            D = ones(size(A)); for i = 1:size(A,3); D(:,:,i) = i; end
             [~,tbl] = anovan(A(:),cat(2,B(:),C(:),D(:)),'varnames',...
                {'Session';'Mice';'Probe #'},'display','off',...
                'continuous',1);
             writecell(tbl,[dirs.figdir ...
                '\Behavior\TrainingIncluded\AcrossDays_probe_FourNeuralMice_' ...
                exlnan_lab{exlnan}  '_' exlStg1_lab{istg} '.xlsx'])
        end
           
        
        % across trials probes - not in paper
        if sum(ismember(figs_to_run,0))>0
            figure; hold on
            for itrialtype = 1:size(alltrial_AllProbes,3)
                dat = squeeze(mean(alltrial_AllProbes...
                    (1:10,:,1:4,itrialtype),2 ...
                    ,'omitnan'));
                smdat = smoothdata(dat,1,'movmean',3);
                errorbar(mean(smdat,2,'omitnan'),std(smdat,[],2,'omitnan')./ ...
                    sqrt(sum(~isnan(dat),2)),'-','LineWidth',2)
            end
             ylim([-.1 1.1])
                xl = get(gca,'xlim');
                xlim([0 xl(2)+1])
                set(gca,'xtick',1:xl(2),'ytick',0:.2:1)
                xlabel('Trial')
                ylabel('Choice index')
                legend('target','nontarget','probe 1','prboe 2','probe 3','Location','best')
            title('Only the 4 mouse used in data analysis, all sessions')
            
            set(gcf,'units','centimeters'); set(gcf,'Position',[0 0 12 9])
            helper_saveandclosefig([dirs.figdir ...
                '\Behavior\TrainingIncluded\AcrossTrials_probe_FourNeuralMice_' exlnan_lab{exlnan}])
        end
    end
end



%%% Notes:
% Trial using these: 'percent_NONtarget_trials_correct',...
%                       'percent_target_trials_correct')
% Dont!
% These are in 10-trial block increments, have nothing to do
% with the real full session data!
% (see testing_learning_data.m for getting and testing
% the real percents out of these variables)