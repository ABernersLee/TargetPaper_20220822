function [neuron_info, neuronlistnew] = get_neuron_info(mousenames,dirs)

neuron_info = []; licks = []; count = 0;
neuronlistnew = cell(1056,1);
repeat_trials = NaN(4,100,2);
probe_trials = NaN(4,100,3); numtrials = [];
for idir = 1:2
    if idir==1; loaddir = dirs.training;
    elseif idir==2; loaddir = dirs.probes; 
    end

    cd(loaddir)
    neuronlist = dir('*.mat');
    n1 = NaN(size(neuronlist,1),5);
    lick1 = NaN(size(neuronlist,1),3);
    for ineuron = 1:size(neuronlist,1)
        load(neuronlist(ineuron).name,'raster_labels') 
        mousenum = find(strcmp(raster_labels.mouse_name,mousenames));
        
        c1 = raster_labels.stimulus_category_with_probes;
        if isfield(raster_labels,'stimulus_category_with_nontarg_repeats')
            c2 = raster_labels.stimulus_category_with_nontarg_repeats;
            probe_repeat_session = ...
                [length(unique(c1(~cellfun(@isempty, c1)))) ...
                length(unique(c2(~cellfun(@isempty, c2))))];
        else           
            probe_repeat_session = ...
                [length(unique(c1(~cellfun(@isempty, c1)))) 0];
        end

        n1(ineuron,:) = [idir mousenum ...
            raster_labels.day_of_behavior probe_repeat_session];

        if probe_repeat_session(:,2)>2
            replabs = unique(c2(~cellfun(@isempty, c2)));
            indlab = find(strncmp('nontarget repeat',replabs,16));
            for ilab = 1:length(indlab)
                repeat_trials(mousenum, raster_labels.day_of_behavior,...
                    ilab) = length(find(...
                    strncmp(replabs(indlab(ilab)), ...
                    c2(~cellfun(@isempty, c2)), ...
                    length(replabs{indlab(ilab)}))));
            end            
        end
            
        if probe_repeat_session(:,1)>2
            probelabs = unique(c1(~cellfun(@isempty, c1)));
            indlab = find(strncmp('probes',probelabs,6));
            for ilab = 1:length(indlab)
                probe_trials(mousenum, raster_labels.day_of_behavior,...
                    ilab) = length(find(...
                    strncmp(probelabs(indlab(ilab)), ...
                    c1(~cellfun(@isempty, c1)),...
                    length(probelabs{indlab(ilab)}))));
            end
        end
        
        L1 = raster_labels.all_licks_port_one; 
        LL1 = cell2mat(cellfun(@min, ...
            L1(~cellfun(@isempty, L1)),'UniformOutput',false));    
        
        L2 = raster_labels.all_licks_port_two;  
        LL2 = cell2mat(cellfun(@min, ...
            L2(~cellfun(@isempty, L2)),'UniformOutput',false));    
       
        lick1(ineuron,:) = [mean([LL1 LL2]) ...
            median([LL1 LL2]) min([LL1 LL2])];
        if rem(ineuron,10)==0
            disp(['idir ' num2str(idir) ', ineuron ' num2str(ineuron)])
        end
        count = count+1;
        neuronlistnew{count} = [loaddir neuronlist(ineuron).name];
        numtrials = cat(1,numtrials,size(c2,1));
    end
    neuron_info = cat(1,neuron_info,n1);
    licks = cat(1,licks,lick1);    
    
end

% this is a repeat neuron with neuron 822 in mouse R session 32. They are
% too eriely similar, maybe clustered two ways
neuron_info(827,:) = [];
neuronlistnew = [neuronlistnew(1:826);neuronlistnew(828:end)];