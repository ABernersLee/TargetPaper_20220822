function plot_accuracy_oversessions(dirs,neuron_info,neuronlistnew,...
    figs_to_run,lick_cutoff_ms)

%%% behavior from sessions/mice with neural recordings


%%% define dirs and labels for plotting and data extraction
if ~isfolder([dirs.figdir '\Behavior\'])
    mkdir([dirs.figdir '\Behavior\'])
end
mousenames = {'Q','R','S','T','V'};

Ls = [111 222 333; 11 22 33];
lab_strs = ...
    {'target','nontarget','probes_comp_1','probes_comp_2','probes_comp_3';...
    'target', 'nontarget nonrepeat', 'nontarget repeat one', ...
    'nontarget repeat two', 'nontarget repeat three';};


%%% get lick and accuracy data out
% there is only one session where there are probes and NT repeats
sessionsall = []; 
lickall = [];
for itype = 1:2 

    %sessions that have either nontarget repeats or probes
    sessionstouse = unique(neuron_info(neuron_info(:,3+...
        itype)>2,2:3),'rows');
    
    trialdatsession = NaN(size(sessionstouse,1),5);   
    
    for isession = 1:size(sessionstouse,1)
        
        %get first neuron of that session's raster data
        neuronlist = find(neuron_info(:,2)==sessionstouse(isession,1) & ...
        neuron_info(:,3)==sessionstouse(isession,2) & ...
        neuron_info(:,3+itype)>2);       
        load(neuronlistnew{neuronlist(1)},'raster_labels')

        %%% Accuracy
        cor = strcmp(raster_labels.correct_vs_error,'correct');
        licked = ~strcmp(raster_labels.lick_choice,'no lick');
        cor_adj = cor & licked;
        % ABL changed 5/10 to only include trials where they made a
        % decision

            %if probe session
        if itype==1
            labels1 = raster_labels.stimulus_category_with_probes;
            %if non-target repeat
        elseif itype==2 
            labels1 = raster_labels.stimulus_category_with_nontarg_repeats;
        end

        trialtypes = NaN(length(labels1),1);
        for itrialtype = 1:5
            ind = strcmp(labels1,lab_strs{itype,itrialtype});

            trialdatsession(isession,itrialtype)=...
                sum(cor_adj(ind)) ./ sum(licked(ind));      

            if itrialtype>2
                trialtypes(ind) = Ls(itype,itrialtype-2);
            else
                trialtypes(ind) = itrialtype;
            end
        end

                
        %%% Lick latency

        %get lick data out for both ports
        L1 = raster_labels.all_licks_port_one; %right
        LL1 = cell2mat(cellfun(@min,L1,'UniformOutput',false));    
        fL1 = ~cellfun(@isempty, L1);
        L2 = raster_labels.all_licks_port_two;  %left
        LL2 = cell2mat(cellfun(@min,L2,'UniformOutput',false));    
        fL2 = ~cellfun(@isempty, L2);

        %fill data in where there is some
        rightleft = NaN(length(labels1),2);
        rightleft(fL1,1) = LL1;
        rightleft(fL2,2) = LL2;

        %take the latency for the correct one (should be first lick anyway)
        mouse_lickright = strcmp(raster_labels.lick_choice,'right');
        mouse_lickleft = strcmp(raster_labels.lick_choice,'left');

        thelicktime = NaN(length(labels1),1);
        thelicktime(mouse_lickright,1) = rightleft(mouse_lickright,1);        
        thelicktime(mouse_lickleft,1) = rightleft(mouse_lickleft,2);

      
        %get things to compare:
        % trialtype, latency of first lick, mouse, accuracy, session, 
        % trial #, type of session   
        lickall = cat(1,lickall,[trialtypes thelicktime ...
            sessionstouse(isession,1)*ones(size(thelicktime)) ...
            cor' sessionstouse(isession,2)*ones(size(thelicktime)) ...
            (1:length(thelicktime))' itype*ones(size(thelicktime))]);  
        
    end
         
    %add session type (itype) to the end and mouse and session #
    %(sessionstouse(:,1:2) to the start
    sessionsall = cat(1,sessionsall,[sessionstouse(:,1:2) ...
        cat(2,trialdatsession,itype*ones(size(trialdatsession,1),1))]);
end


%only mice with neural data on at least one day
mice = 1:4;

%%% figure plots


% licks, speed-accuracy tradeoff for sup. fig 8
if sum(ismember(figs_to_run,[0 18]))>0
    helper_plot_lick_lat(lickall,dirs,lick_cutoff_ms)
end

% licks more ways, figures 5, s8 and s9
if sum(ismember(figs_to_run,[0 5 18 19]))>0
    helper_plot_lick(lickall,mice,mousenames,dirs,figs_to_run,lick_cutoff_ms)
end

% accuracy
if sum(ismember(figs_to_run,[0 5 7]))>0
    helper_plot_behavior_sessions(sessionsall,mice,mousenames,dirs,...
        figs_to_run)
end
