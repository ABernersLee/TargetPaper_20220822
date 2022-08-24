function [trialtypes, cor,isProbe,isRepeat] = get_trial_info(raster_labels)

isProbe = false; isRepeat = false;

%find out what kind of session this is (probe vs repeat) 
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

% define trialtypes variable
trialtypes = ...
    NaN(length(raster_labels.stimulus_category_with_probes),1);

%If it is a repeat session name the trial types 1,2,11,22,33
if probe_repeat_session(2)>2
    isRepeat = true;
    labels1 = raster_labels.stimulus_category_with_nontarg_repeats;
    t = strcmp(labels1,'target');
    nt = strcmp(labels1,'nontarget nonrepeat');        
    nt1 = strcmp(labels1,'nontarget repeat one');        
    nt2 = strcmp(labels1,'nontarget repeat two');        
    nt3 = strcmp(labels1,'nontarget repeat three');
    trialtypes(t) = 1;  trialtypes(nt) = 2; 
    trialtypes(nt1) = 11; trialtypes(nt2) = 22; 
    trialtypes(nt3) = 33;                
end

%If it is a probe session name the trial types 1,2,111,222,333
if probe_repeat_session(1)>2
    labels1 = raster_labels.stimulus_category_with_probes;
    isProbe = true;
    t = strcmp(labels1,'target');
    nt = strcmp(labels1,'nontarget');        
    p1 = strcmp(labels1,'probes_comp_1');        
    p2 = strcmp(labels1,'probes_comp_2');        
    p3 = strcmp(labels1,'probes_comp_3');

    if probe_repeat_session(2)<=2
        trialtypes(t) = 1;  trialtypes(nt) = 2; 
    end
    trialtypes(p1) = 111; trialtypes(p2) = 222;  trialtypes(p3) = 333;
end           

%if neuron is in a session without repeats or probes (can still include
%it in target/nontarget analyses
if probe_repeat_session(1)<=2 && probe_repeat_session(2)<=2
   labels1 = raster_labels.stimulus_category_with_probes;
   t = strcmp(labels1,'target');
   nt = strcmp(labels1,'nontarget');
   trialtypes(t) = 1;  trialtypes(nt) = 2; 
end

%rewarded/correct trials
cor = strcmp(raster_labels.correct_vs_error,'correct')';

%trialdat is trialtype and whether it was rewarded/correct
% trialdat = [trialtypes cor'];