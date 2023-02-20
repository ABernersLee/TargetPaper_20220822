function calculate_plot_PVcorr(dirs,neuron_info,neuronlistnew, ...
    excludeLickNeurons,numNcutoff,postdattimes,ExclLickTime,is_issig_lick...
    ,PVCorr_lookback,runshuffle,plotPVsessions,figs_to_run,...
    tozscore_neuron,tozscore_pv,prop_selective)

%%%%%
%%%%%
%%%%% This script does population-level analyses on the piriform data for
%%%%% each session. It does the population vector correlation (PVcorr) for
%%%%% each session across trials broken down by type. It also does a simple
%%%%% linear discriminant analysis to decode the trial type. It uses
%%%%% shuffles to test the significance of these effects as well as effects
%%%%% of how these aspects change across sessions.
%%%%%
%%%%%

all_labels.plotPVsessions = plotPVsessions;
all_labels.PVCorr_lookback = PVCorr_lookback;
all_labels.runshuffle = runshuffle;

%only figure 6 needs the shuffle
if sum(ismember(figs_to_run,[0 6]))==0
    all_labels.runshuffle = false;
end

%{set up parameters}
all_labels.numshuff = 1000;
all_labels.LB = [3:20 25:5:50];
%10-fold decoding, downsmapling to the fewest 
% of the repeats/probes
all_labels.numfolds = 10;

%%%% if downsampling
all_labels.numfolds2 = 20; 
all_labels.numexlund = 4;
all_labels.todownsample = false;
%%% iterations of randomly doing
% 10-fold decoding with different downsampled populations. if todownsample
% is false then this becomes 1.
all_labels.minhmin = 20; % 15;

%%%% If excluding lick neurons
if excludeLickNeurons>0    
    islick = is_issig_lick(:,1);
    islicksig = is_issig_lick(:,2);
end

%{set up figure labels}
all_labels.mousenames = {'Q';'R';'S';'T'};
all_labels.typelab = {'Probes','NontargetRepeats'};
Ls = [111 222 333; 11 22 33];  % probe and nontarget markers

%%%% Name based on if excluding lick neurons
if excludeLickNeurons==1
    all_labels.addon = [ '_Post_' num2str(postdattimes(1)) '_' ...
        num2str(postdattimes(2)) '_ExlLickNeruons_SessionNeuronCutoff'...
        num2str(numNcutoff)];
elseif excludeLickNeurons==2
    all_labels.addon = [ '_Post_' num2str(postdattimes(1)) '_' ...
        num2str(postdattimes(2)) ...
        '_ExlLickSigNeruons_SessionNeuronCutoff' num2str(numNcutoff)];
else
    all_labels.addon = [ '_Post_' num2str(postdattimes(1)) '_' ...
        num2str(postdattimes(2)) '_SessionNeuronCutoff' ...
        num2str(numNcutoff)];
end

%%%% Name based on if excluding lick time
if ExclLickTime==1
    all_labels.addon = [all_labels.addon  '_ExlLick'];
end

%%%% Name based on normalization
all_labels.addon = [all_labels.addon '_PVCorrLookback_' ...
    num2str(PVCorr_lookback)];

%%%% if downsample
if all_labels.todownsample
    all_labels.addon = [all_labels.addon '_DownsampleTrials' ...
        num2str(minhmin)];
end

%added conditionals 1/5/23 after reviewer comment
all_labels.tozscore_neuron = tozscore_neuron;
all_labels.tozscore_pv = tozscore_pv;

if all_labels.tozscore_pv==1
     all_labels.addon = [all_labels.addon '_PVzscored'];
else
     all_labels.addon = [all_labels.addon '_PVraw'];
end

%added conditionals 1/5/23 after reviewer comment
if all_labels.tozscore_neuron==1
     all_labels.addon = [all_labels.addon '_NSzscored'];
else
     all_labels.addon = [all_labels.addon '_NSraw'];
end

try
    
% there is only one session where there are both probes and NT repeats
% (mouse 4, session 22)
% figure 6 is NTrepeats/over-training data
% figure 7 is Probe data and NTrepeat data
for itype =  2:-1:1 

    if sum(ismember(figs_to_run,[0 7 17]))==0 && itype==1
        continue
    end

    if sum(ismember(figs_to_run,[0 6 7 17]))==0 && itype==2
        continue
    end
    
    % find the sessions that have this type of trial (probe or repeat)
    sessionstouse = ...
        unique(neuron_info(neuron_info(:,3+itype)>2,2:3),'rows');
    
    %define some empty variables
    corr_repeat_trial_session = [];
    corr_repeats = [];
    decode_fold_label = []; decode_fold_label_cat = [];
    decoded_fold_group_all = [];
    decode_fold_shuffle = [];
    decode_fold_cat_shuffle = [];
    decode_fold_shuffle_prob = [];
    decode_fold_cat_shuffle_prob = [];
    across_lookback= [];
    
    for isession = 1:size(sessionstouse,1)
        
        if isession==4 && itype==2
            disp('SKIPPING')
            continue
        end
           
        
        %empty variables to fill across neurons
        neurondatsession = [];
        fr_neuron = [];
        trialdatsession = [];
        
        % get the neurons we need for this session
        neuronlist = find(neuron_info(:,2)==sessionstouse(isession,1)...
            & neuron_info(:,3)==sessionstouse(isession,2) &...
            neuron_info(:,3+itype)>2);
        
        %exclude lick neurons if applicable
        if excludeLickNeurons==1
            neuronlist = neuronlist(~ismember(neuronlist,find(islick)));
        elseif excludeLickNeurons==2            
            neuronlist = neuronlist(~ismember(neuronlist,find(islicksig)));
        end
        
        % If fewer neurons than the cutoff, skip this session
        if length(neuronlist)<numNcutoff             
            disp(['too few, ns = ' num2str(length(neuronlist))...
                ', session ' num2str(isession) ' ' ...
                all_labels.typelab{itype}])
            continue
        end
        


        %%% get data from across all neurons
        for ineuron = 1:length(neuronlist)
       
            %load data
            load(neuronlistnew{neuronlist(ineuron)},...
                'raster_labels','raster_data')
            
            %Get and organize the lick times                        
            thelicktime = get_licks(raster_labels);                
        
            
            % get trial types and correct/incorrect
            clear trialdat
            [trialdat(:,1),trialdat(:,2)] = get_trial_info(raster_labels);
            
            %{get neural data
            %500ms before 
            predat = sum(raster_data(:,4501:5000),2)./.5;     
            %determined by variable (e.g. 500ms after, 100-600)
            postdat = ...
                sum(raster_data(:,postdattimes(1)+5001:...
                postdattimes(2)+5000),2)./(length(postdattimes(1)+5001:...
                postdattimes(2)+5000)/1000); 

            %If excluding lick times, exclude any trials that have a lick
            %before the end of the post-odor period 
            % (e.g. 600ms after odor onset)            
            if ExclLickTime
                exl = thelicktime<postdattimes(2);
                if ineuron==1
                    disp(['Exclude ' num2str(sum(exl))])
                end

                predat(exl,:) = [];
                postdat(exl,:) = [];
                trialdat(exl,:) = [];
                thelicktime(exl) = [];
            end
            
            %making sure that the trialdat looks the same for each neuron
            %loaded within that session
            if ineuron==1
                %this is session's trial dat, set by first neuron
                trialdatsession = trialdat; 
            else
                %check that is the same as all other neurons in that
                %session
                if size(trialdat,1)~=size(trialdatsession,1)                    
                    error('not matching 1')
                elseif sum(sum(isnan(trialdat))) ~=...
                        sum(sum(isnan(trialdatsession)))
                    error('not matching 2')                        
                elseif sum(sum(trialdat~=trialdatsession)) > ...
                        sum(sum(isnan(trialdat)))
                    error('not matching 3')
                end                
            end
          
            %If a neuron has at least 5 spikes in the pre and post period
            %then they are included.
            
            %Any neuron with fewer than five spikes in the post-odor window
            % is excluded later, so this may be redundant
            if sum([predat;postdat])>5
                %then save out the neuron with the session
               neurondatsession =  cat(3, neurondatsession,...
                   [predat postdat]);

               fr = sum(raster_data(:,5201:5600),2)./.4;
               frT = mean(fr(trialdatsession(:,1)==1));
               frNT = mean(fr(trialdatsession(:,1)==2));
               fr_neuron = cat(1,fr_neuron,[frT frNT]);                 
            end
            
        end
        

        %%% applying some preprocessing to the neural data before 
        %%% subsequent decoding or other population analses

        %Using post-odor bin as the analysis window 
            %Since we are using pearsons we should be ok with neurons that
            % tend to fire more or not in general
        ns = squeeze((neurondatsession(:,2,:)));    

        %Some repeat sessions only have 1 or 2 repeats, so exclude the
        %labels of the trial types of 33 and/or 22
        ts = [1 2 Ls(itype,:)]; % The trial types for this type of session
        if itype==2 && ...
                sum(ismember(unique(trialdatsession(~isnan(...
                trialdatsession(:,1)),1)),ts))<=4
            ts = ts(1:end-1);
        end
                        
        %take only the sessions that have the trial types of interest (this
        %is only excluding things in one session that has both trial types)
        ns(~ismember(trialdatsession(:,1),ts),:) = [];               
        thelicktime(~ismember(trialdatsession(:,1),ts)) = [];
        trialdatsession(~ismember(trialdatsession(:,1),ts),:) = [];
        
        
        %Any neuron with fewer than five spikes in the post-odor window is
        %excluded
        ns(:,sum(ns)<5) = [];         
        %The data is z-scored for all analyses  

        ns_save = ns;

        %added conditional 1/5/23 after reviewer comment
        if all_labels.tozscore_neuron==1
%             ns = zscore(ns,1); %scale it for decoding %this was using
%             population standard deviation (assums this is all the values,
%             not a sample)
            ns = zscore(ns,[],1); %scale it for decoding using sample population deviation
        end

        %If there are fewer than the cutoff amount of neurons this session
        %then dont use this session
        if size(ns,2)<numNcutoff 
            disp(['too few, ns = ' num2str(size(ns,2))])
            continue
        end
        

        
        %%% PV CORRELATIONS        

        % run correlations and plot if applicable
        [corr_repeat_trial_session0,across_lookback0,corr_repeats0] = ...
            get_pv_corr(ns,trialdatsession,dirs,...
            all_labels.mousenames{sessionstouse(isession,1)},...
            sessionstouse(isession,2),all_labels.typelab{itype},ts,...
            isession,sessionstouse,all_labels,figs_to_run);

        across_lookback = cat(1,across_lookback,across_lookback0);

        %add correlations broken down by trial type
        corr_repeat_trial_session = cat(1,...
            corr_repeat_trial_session,[corr_repeat_trial_session0 ...
            thelicktime ones(size(thelicktime))*size(ns,2)]);
        %add number of neurons to the end
        if ~isempty(corr_repeats0)
            corr_repeats = cat(1,corr_repeats,[corr_repeats0 size(ns,2)]);
        end
       %plot
        if all_labels.plotPVsessions && itype==1 && ...
                (sum(ismember(figs_to_run,0))>0 ...
                || (sum(ismember(figs_to_run,7))>0 && ...
                strcmp('T',all_labels.mousenames{sessionstouse...
                (isession,1)}) ...
                && sessionstouse(isession,2)==26))

            %only T26 for figure 7

            plot_PV_Sessions(trialdatsession,all_labels.typelab{itype}, ...
                all_labels.mousenames{sessionstouse(isession,1)}, ...
                dirs,isession,sessionstouse,corr_repeat_trial_session0, ...
                thelicktime,all_labels,figs_to_run)
        end



        %%% DECODING
        tic

        [decoded_fold_group, decoded_fold, decoded_fold_cat, ...
            decoded_fold_S, decoded_fold_cat_S, decoded_fold_Sprob, ...
            decoded_fold_cat_Sprob] = ...
            get_decoding_seperate_trialtypes...
            (ns_save,trialdatsession,isession,all_labels,itype);
        %gets zscored within the function
        
        t = toc;
        disp(['Done in ' num2str(t/60) ' min, session ' ...
            num2str(isession) ' ' all_labels.typelab{itype}])
        
        %take the average over the permuations
        decoded_fold_group = mean(mean(decoded_fold_group,3,'omitnan')...
            ,4,'omitnan');
        decoded_fold =  mean(mean(decoded_fold,3,'omitnan'),4,'omitnan'); 
        decoded_fold_cat = mean(mean(decoded_fold_cat,3,'omitnan')...
            ,4,'omitnan');           
        decoded_fold_S =  mean(mean(decoded_fold_S,3,'omitnan')...
            ,4,'omitnan');
        decoded_fold_cat_S = mean(mean(decoded_fold_cat_S,3,'omitnan')...
            ,4,'omitnan');      
        decoded_fold_Sprob =  mean(mean(decoded_fold_Sprob,3,'omitnan')...
            ,4,'omitnan');
        decoded_fold_cat_Sprob = ...
            mean(mean(decoded_fold_cat_Sprob,3,'omitnan'),4,'omitnan');
        
        if size(decoded_fold,1)~=length(thelicktime)
            error('Something is wrong with the lick data')
        end
         
         
%%%%%%add to the variables each session

%repeat type, decoding correct, posterior probability that it is NT [that 
% it is the right one, I make it P(correct) later], trial number, session, 
% mouse, correct, licktime
       decode_fold_label = cat(1,decode_fold_label,...
           [decoded_fold (1:size(decoded_fold,1))' ... 
                    ones(size(decoded_fold,1),1)*...
                    sessionstouse(isession,2) ...
                    ones(size(decoded_fold,1),1)*...
                    sessionstouse(isession,1)...
                    trialdatsession(:,2) thelicktime]); 
       decode_fold_label_cat = cat(1,decode_fold_label_cat,...
           [decoded_fold_cat...
           (1:size(decoded_fold,1))' ... 
                    ones(size(decoded_fold,1),1)*...
                    sessionstouse(isession,2) ... 
                    ones(size(decoded_fold,1),1)*...
                    sessionstouse(isession,1)...
                    trialdatsession(:,2) thelicktime]);

       decoded_fold_group_all = cat(1,decoded_fold_group_all,...
           decoded_fold_group);
       decode_fold_shuffle = cat(1,decode_fold_shuffle,decoded_fold_S);
       decode_fold_cat_shuffle = cat(1,decode_fold_cat_shuffle,...
           decoded_fold_cat_S);
       
       decode_fold_shuffle_prob = cat(1,decode_fold_shuffle_prob,...
           decoded_fold_Sprob);
       decode_fold_cat_shuffle_prob = ...
           cat(1,decode_fold_cat_shuffle_prob,decoded_fold_cat_Sprob);
                     
    end
    
    %%% Making this the probability of the decoder choosing the correct 
    %%% category (was probability of NT choice) so switched it around for 
    %%% the target.
    decode_fold_label_cat(decode_fold_label_cat(:,1)==1,3) = ...
        1-decode_fold_label_cat(decode_fold_label_cat(:,1)==1,3);    
    decode_fold_cat_shuffle_prob(decode_fold_label_cat(:,1)==1,:) = ...
        1-decode_fold_cat_shuffle_prob(decode_fold_label_cat(:,1)==1,:);
   
    %%%%%%%% Labels for plots
    all_labels.ittlabs = {'self';'target';'nontarget';...
        all_labels.typelab{itype}};         
    all_labels.lablabs = {'Target';'NonTarget';...
        [all_labels.typelab{itype} ' 1'];...
        [all_labels.typelab{itype} ' 2'];...
        [all_labels.typelab{itype} ' 3']};        
    all_labels.labss = {'Target';'NonTarget';[all_labels.typelab{itype}]};

        
    %%%%%%%% DECODING FIGURES

    %figure 6 figures
    if all_labels.runshuffle && ...
            ((sum(ismember(figs_to_run,[0 6]))>0 && itype==2) ||...
            (sum(ismember(figs_to_run,0))>0))

        % Plot the decoding difference from shuffle, and significance 
        % in Decoding_SigShuffle_FromShuffle folder
        plot_Decoding_SigShuffle_FromShuffle(decode_fold_shuffle, ...
        decode_fold_label, ...
        decode_fold_cat_shuffle,decode_fold_label_cat,...
        dirs,all_labels,decode_fold_shuffle_prob,...
        decode_fold_cat_shuffle_prob,...
        itype, figs_to_run)
    

        % Anovas for decoding, in Decoding_ANOVAS folder
        plot_Decoding_ANOVAS(decode_fold_label, decode_fold_label_cat,...
        dirs,all_labels,itype, figs_to_run)

        % plotting the decoding probabilities a lot of different ways 
        %Across Sessions
        plot_Decoding_Across(2,decode_fold_label,decode_fold_label_cat...
            ,all_labels,dirs,itype, figs_to_run)
    end
       
                
    
    
    %%%%%%%% POPULATION VECTOR ANALYSIS
    
    % How do latency, T/NT decoding, and accuracy relate?
    % Figure 7
    if sum(ismember(figs_to_run,[0 7]))>0 && itype==1
        plot_TNT_and_Latency...
            (corr_repeat_trial_session,all_labels,dirs,...
            itype,across_lookback,figs_to_run)        
    end 

    % Average RSA values - supplemental figure 7
    if sum(ismember(figs_to_run,[0 17]))>0
        plot_RSA_Averages(corr_repeats,all_labels,dirs,itype,...
            prop_selective) 
        % prop_selective added 1/15/23 for reviewer comments
    end
   
        
end

catch ME
   ME.getReport
   disp('*')
end
