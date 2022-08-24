function [selectivity_dat,selectivity_params] = get_sig_neurons_all(odorwindowlabel,neuron_info,odorselective,odorselective_all,odorsel,odorsel_all)

% for each timebin the significance threshold is 
% p_cutoff/# of timebins tested
sp.p_cutoff = .05;

%define time periods and significance threshold   
sp.OdorIndLarge = odorwindowlabel(:,1)>=0 & odorwindowlabel(:,2)<=2000;        
sp.OdorIndSmall = odorwindowlabel(:,1)>=0 & odorwindowlabel(:,2)<=1000;
sp.LickInd = odorwindowlabel(:,1)>=-400 & odorwindowlabel(:,2)<=600;    
sp.OdorOffInd = odorwindowlabel(:,1)>=2000 & odorwindowlabel(:,2)<=3000;

%%
%%%%% AVERAGE ACROSS SESSIONS 
% all_MeanSes

%checks if has neuron here with values
sp.dem_num = sum(sum(sum(~isnan(odorselective),6),5),4)==0;  

%lick
%gets the lickspecific trials (rawavgLicktrials) modulation 
% [in meanLick2] and 
meanLick = odorselective(:,:,:,sp.LickInd,:,2);              
lickbins2 = meanLick(:,:,:,:,1); 
lickbins2(meanLick(:,:,:,:,2)>=sp.p_cutoff) = NaN;
[~,minbins3] = max(abs(lickbins2),[],4,'omitnan');
meanLick2 = NaN(size(minbins3));
for ix = 1:size(minbins3,1)
    for iy = 1:size(minbins3,2)
        meanLick2(ix,iy,:) = arrayfun(@(x) ...
            meanLick(ix,iy,x,minbins3(ix,iy,x),3), 1:size(minbins3,3));
    end
end            


%odor-on modulation     
% 1 in dim=6 here is odor
meanOdor = odorselective(:,:,:,sp.OdorIndLarge,:,1); 
%NaN(4,38,30,numbins,4,19); 
%rats, days, neurons, bin, modulation/p-value/rawavgLicktrials/rawavg,
% odor/lick/target/correct/probe/repeat     
% dim=5 here is [mod1 p1 odorbin2 odorbin]; 
    % mod1 = (postdat-predat)./(postdat+predat);
odorbins2 = meanOdor(:,:,:,:,1); 
odorbins2(meanOdor(:,:,:,:,2)>=sp.p_cutoff) = NaN;
%max sig odor modulation across bins for each neuron
[xx,minbins4X1] = max(odorbins2,[],4,'omitnan');
%min sig odor modulation acorss bins for each neuron
[~,minbins4X2] = min(odorbins2,[],4,'omitnan');
minbins4 = minbins4X1;
%take the max ind when it isn't NaN or <=0, then take the ind with the min
%so neurons with a sig + and sig - bin are classified as sig +
minbins4(xx<=0 | isnan(xx)) = minbins4X2(xx<=0 | isnan(xx));

%gets the lickspecific trials (rawavgLicktrials) [in meanOdor2] and 
% whether the modulation is up or down for each neuron [meanLickUpDown] 
% at the ind of its peak modulation 
meanOdor2 = NaN(size(minbins4));        
meanLickUpDown = false(size(minbins4));
for ix = 1:size(minbins4,1)
    for iy = 1:size(minbins4,2)
        meanOdor2(ix,iy,:) = arrayfun(@(x) ...
            meanOdor(ix,iy,x,minbins4(ix,iy,x),3), 1:size(minbins4,3));                
        meanLickUpDown(ix,iy,:) = arrayfun(@(x) ...
            meanOdor(ix,iy,x,minbins4(ix,iy,x),1)>0, 1:size(minbins3,3));
    end
end


% first for upward, then downwardly odor-on modulated neurons, finding
% neurons that have higher average raw firing when lick-aligned than
% odor-aligned, in trials where there is over 800ms (sepcut in
% calculate_plot_singleneuron_selectivity.m) between the odor-on and the
% lick
islickMean = false(size(minbins4));    
islickMean(meanLickUpDown) = (meanLick2(meanLickUpDown) > ...
    (meanOdor2(meanLickUpDown)*1.0))'; 
%played with this cutoff, this looked good
islickMean(~meanLickUpDown & ~sp.dem_num) = ...
    (meanLick2(~meanLickUpDown & ~sp.dem_num) < ...
    (meanOdor2(~meanLickUpDown & ~sp.dem_num)*1.0))'; 
%played with this cutoff (ratio of 1), this looked good

meanOdor = odorselective(:,:,:,sp.OdorIndSmall,:,1); 
meanOdorp = any(meanOdor(:,:,:,:,2)<(sp.p_cutoff/sum(sp.OdorIndSmall)),4);
isSigOdorMean = meanOdorp & ~islickMean;

isSigOdorMean = double(isSigOdorMean);
isSigOdorMean(sp.dem_num==1) = NaN;

odorbins2 = meanOdor(:,:,:,:,1);
odorbins2(meanOdor(:,:,:,:,2)>=(sp.p_cutoff/sum(sp.OdorIndSmall))) = NaN;
isUpMean = ~islickMean & isSigOdorMean==1 & any(odorbins2>0,4);
isDownMean = ~islickMean & isSigOdorMean==1 & any(odorbins2<0,4) & ...
    ~isUpMean;    

islicksigMean = islickMean & any(meanLick(:,:,:,:,2) < ...
    (sp.p_cutoff/sum(sp.LickInd)),4);
islickUPmean = islicksigMean & meanLickUpDown;
islickDOWNmean = islicksigMean & ~meanLickUpDown & ~sp.dem_num;

%Target
Tdat = odorselective(:,:,:,sp.OdorIndSmall,:,3); 
%NaN(4,38,30,numbins,4,14); 
%rats, days, neurons, bin, p-value/modulation/rawavg,odor/
% lick/target/correct/probe/repeat
is_targetMean= any(Tdat(:,:,:,:,2) < ...
    (sp.p_cutoff/sum(sp.OdorIndSmall)) & Tdat(:,:,:,:,1)>0,4);
is_nontargetMean= any(Tdat(:,:,:,:,2) < ...
    (sp.p_cutoff/sum(sp.OdorIndSmall)) & Tdat(:,:,:,:,1)<0,4);
[lookx,lookyz] = find(is_targetMean & is_nontargetMean);
[looky,lookz] = ind2sub(size(Tdat,2:3),lookyz);
[~,mm] = min(Tdat(lookx,looky,lookz,:,2),[],4);
maxT = arrayfun(@(x) ...
    Tdat(lookx(x),looky(x),lookz(x),mm(x),1),1:length(lookx));
is_targetMean(lookx(maxT<0),looky(maxT<0),lookz(maxT<0)) = false;
is_nontargetMean(lookx(maxT>0),looky(maxT>0),lookz(maxT>0)) = false;
is_targetTrialCorrectedMean= any(odorselective(:,:,:,sp.OdorIndSmall,2,10)...
    < (sp.p_cutoff/sum(sp.OdorIndSmall)),4);

%Error (same except first line)
Tdat = odorselective(:,:,:,sp.OdorIndSmall,:,4); 
%NaN(4,38,30,numbins,3,9); 
%rats, days, neurons, bin, p-value/modulation/rawavg,odor/
% lick/target/correct/probe/repeat
is_correctMean= any(Tdat(:,:,:,:,2)<(sp.p_cutoff/sum(sp.OdorIndSmall))...
    & Tdat(:,:,:,:,1)>0,4);
is_errorMean= any(Tdat(:,:,:,:,2)<(sp.p_cutoff/sum(sp.OdorIndSmall))...
    & Tdat(:,:,:,:,1)<0,4);
[lookx,lookyz] = find(is_correctMean & is_errorMean);
[looky,lookz] = ind2sub(size(Tdat,2:3),lookyz);
[~,mm] = min(Tdat(lookx,looky,lookz,:,2),[],4);
maxT = arrayfun(@(x) Tdat(lookx(x),looky(x),lookz(x),mm(x),1),1:length(lookx));
is_correctMean(lookx(maxT<0),looky(maxT<0),lookz(maxT<0)) = false;
is_errorMean(lookx(maxT>0),looky(maxT>0),lookz(maxT>0)) = false;


%Probe
Tdat = squeeze(odorselective(:,:,:,sp.OdorIndSmall,:,15:17));
is_probeMean = any(any(Tdat(:,:,:,:,2,:) < ...
    (sp.p_cutoff/(sum(sp.OdorIndSmall)*3)),6),4);
session_probeMean = sum(sum(~isnan(odorselective(:,:,:,:,2,5:7))...
    ,6,'omitnan'),4,'omitnan')>0;

is_probeMean2 = squeeze(any(Tdat(:,:,:,:,2,:) < ...
    (sp.p_cutoff/(sum(sp.OdorIndSmall))),4));
session_probeMean2 = squeeze(sum(~isnan(odorselective(:,:,:,:,2,5:7))...
    ,4,'omitnan'))>0;

%Repeats
Tdat = squeeze(odorselective(:,:,:,sp.OdorIndSmall,:,18:19));
is_repeatMean = any(any(Tdat(:,:,:,:,2,:) < ...
    (sp.p_cutoff/(sum(sp.OdorIndSmall)*2)),6),4);
session_repeatMean = sum(sum(~isnan(odorselective(:,:,:,:,2,8:9))...
    ,6,'omitnan'),4,'omitnan')>0;
        
is_repeatMean2 = squeeze(any(Tdat(:,:,:,:,2,:) < ...
    (sp.p_cutoff/(sum(sp.OdorIndSmall))),4));
session_repeatMean2 = squeeze(sum(~isnan(odorselective(:,:,:,:,2,8:9))...
    ,4,'omitnan'))>0;

%Odor Off

odoroff1 = odorselective(:,:,:,sp.OdorOffInd,:,1);    
odoroff12 = odoroff1(:,:,:,:,1); 
odoroff12(odoroff1(:,:,:,:,2)>=(sp.p_cutoff/sum(sp.OdorOffInd))) = NaN;
[~,maxOffInd] = max(abs(odoroff12),[],4);
maxOff = NaN(size(maxOffInd));
for ix = 1:size(maxOffInd,1)
    for iy = 1:size(maxOffInd,2)
        maxOff(ix,iy,:) = arrayfun(@(x) ...
            abs(meanOdor(ix,iy,x,maxOffInd(ix,iy,x),3)), ...
            1:size(maxOffInd,3));
    end
end

odoroff2 = squeeze(odorselective_all(:,:,:,6,:));   

isOdorOffMean = any(odoroff1(:,:,:,:,2) < ...
    (sp.p_cutoff/sum(sp.OdorOffInd)),4) & ...
    odoroff2(:,:,:,2) < sp.p_cutoff & ...
    ((maxOff>0) == (odoroff2(:,:,:,1)>0)) & ...
    ((maxOff<0) == (odoroff2(:,:,:,1)<0));    

isOdorOff_UpMean = any(odoroff1(:,:,:,:,2) < ...
    (sp.p_cutoff/sum(sp.OdorOffInd)),4) & ...
    odoroff2(:,:,:,2)<sp.p_cutoff & maxOff>0 & odoroff2(:,:,:,1)>0;

%Odor Off target/nontarget           
Tdat = odorselective(:,:,:,sp.OdorOffInd,:,3);
is_target_offMean= any(Tdat(:,:,:,:,2) < ...
    (sp.p_cutoff/sum(sp.OdorOffInd)) & Tdat(:,:,:,:,1)>0,4);
is_nontarget_offMean= any(Tdat(:,:,:,:,2) < ...
    (sp.p_cutoff/sum(sp.OdorOffInd)) & Tdat(:,:,:,:,1)<0,4);
[lookx,lookyz] = find(is_target_offMean & is_nontarget_offMean,3);
[looky,lookz] = ind2sub(size(Tdat,2:3),lookyz);
[~,mm] = min(Tdat(lookx,looky,lookz,:,2),[],4);
maxT = arrayfun(@(x) ...
    Tdat(lookx(x),looky(x),lookz(x),mm(x),1),1:length(lookx));
is_target_offMean(lookx(maxT<0),looky(maxT<0),lookz(maxT<0)) = false;
is_nontarget_offMean(lookx(maxT>0),looky(maxT>0),lookz(maxT>0)) = false;

all_MeanSes = cat(4,isSigOdorMean, islicksigMean, isUpMean, ...
    isDownMean, is_targetMean, is_nontargetMean,  ...
    is_correctMean, is_errorMean, is_probeMean, session_probeMean,...
    is_repeatMean, session_repeatMean, ...
    isOdorOffMean, isOdorOff_UpMean, is_target_offMean, ...
    is_nontarget_offMean, islickUPmean, islickDOWNmean);            



%%
%%%%% BEST SESSION OF EACH MOUSE
% all_BestSes
%mm is the session with the most neurons
% I don't end up using this at all in any of the paper figures or text
[~,mm] = max(sum(~sp.dem_num,3),[],2); 

neuron_info(:,7) = zeros(size(neuron_info,1),1);
all_BestSes = NaN(size(all_MeanSes,1),...
    size(all_MeanSes,3),size(all_MeanSes,4));
for imouse = 1:4
    all_BestSes(imouse,:,:) = squeeze(all_MeanSes(imouse,mm(imouse),:,:));
    neuron_info(neuron_info(:,2) == ...
        imouse & neuron_info(:,3)==mm(imouse),7) = 1;
end


%%
%%%%% EACH NEURON AS INDIVIDUAL (concatonated as if they were all
%%%%% different, independent neurons)
% all_Neurons
     %odorsel = odor/lick/target/correct/probex3/repeatx2

%Odor and lick modulated
odorbins = squeeze(odorsel(:,sp.OdorIndLarge,1,:));
odorbins2 = odorbins(:,:,1); odorbins2(odorbins(:,:,2)>=sp.p_cutoff) = NaN;
[xx,minbinsX1] = max(odorbins2,[],2,'omitnan');
[~,minbinsX2] = min(odorbins2,[],2,'omitnan');
minbins = minbinsX1;
minbins(xx<=0 | isnan(xx)) = minbinsX2(xx<=0 | isnan(xx));

absmod_odor = arrayfun(@(x) odorbins(x,minbins(x),3),1:length(minbins));  
absmod_lickUpDown = arrayfun(@(x) ...
    odorbins(x,minbins(x),1)>0,1:length(minbins));  

lickbins = squeeze(odorsel(:,sp.LickInd,2,:));
lickbins2 = lickbins(:,:,1); lickbins2(lickbins(:,:,2)>=sp.p_cutoff) = NaN;
[~,minbins2] = max(abs(lickbins2),[],2,'omitnan');

absmod_lick = arrayfun(@(x) lickbins(x,minbins2(x),3),1:length(minbins2));    

%played with this cutoff (ratio of 1), this looked good   
islick = false(size(minbins2));    
islick(absmod_lickUpDown) = (absmod_lick(absmod_lickUpDown) > ...
    (absmod_odor(absmod_lickUpDown)*1.0))'; 
islick(~absmod_lickUpDown) = (absmod_lick(~absmod_lickUpDown) < ...
    (absmod_odor(~absmod_lickUpDown)*1.0))';

%%% for checking some
%{ sh = [259 255 420 445 931 746 451 474 492 882 931]; 
%ones that should be
%sum(islick(sh))./length(sh)
%shn = [3 11 85 373 372 632 532 543 654 762 795 1052]; 
%ones that should not be
%sum(islick(shn))./length(shn)

isSigOdor1 = any(odorbins(:,1:5,2)<(sp.p_cutoff/5),2);
isSigOdor = isSigOdor1 & ~islick;
odorbins2 = odorbins(:,1:5,1);
odorbins2(odorbins(:,1:5,2)>=(sp.p_cutoff/5)) = NaN;
isUp = ~islick & isSigOdor & any(odorbins2>0,2);
isDown = ~islick & isSigOdor & any(odorbins2<0,2) & ~isUp;

islicksig = islick & any(lickbins(:,:,2)<(sp.p_cutoff/sum(sp.LickInd)),2);  

is_issig_lick = [islick islicksig];
islickup = islicksig & absmod_lickUpDown';
islickdown = islicksig & ~absmod_lickUpDown';

%target
Tdat = squeeze(odorsel(:,sp.OdorIndSmall,3,:));
is_target= any(Tdat(:,:,2)<(sp.p_cutoff/sum(sp.OdorIndSmall)) & ...
    Tdat(:,:,1)>0,2);
is_nontarget= any(Tdat(:,:,2)<(sp.p_cutoff/sum(sp.OdorIndSmall)) & ...
    Tdat(:,:,1)<0,2);
look = find(is_target & is_nontarget);
[~,mm] = min(Tdat(look,:,2),[],2);
maxT = arrayfun(@(x) Tdat(look(x),mm(x),1),1:length(look));
is_target(look(maxT<0)) = false;
is_nontarget(look(maxT>0)) = false;

%error
Tdat = squeeze(odorsel(:,sp.OdorIndSmall,4,:));
is_correct= any(Tdat(:,:,2)<(sp.p_cutoff/5) & Tdat(:,:,1)>0,2);
is_error= any(Tdat(:,:,2)<(sp.p_cutoff/5) & Tdat(:,:,1)<0,2);
look = find(is_correct & is_error);
[~,mm] = min(Tdat(look,:,2),[],2);
maxT = arrayfun(@(x) Tdat(look(x),mm(x),1),1:length(look));
is_correct(look(maxT<0)) = false;
is_error(look(maxT>0)) = false;


%probe
Tdat = squeeze(odorsel(:,sp.OdorIndSmall,15:17,:));
is_probe= any(any(Tdat(:,:,:,2)<(sp.p_cutoff/(sum(sp.OdorIndSmall)*3)),2),3);
session_probe = neuron_info(:,4)>2;        
is_probe2= squeeze(any(Tdat(:,:,:,2)<(sp.p_cutoff/(sum(sp.OdorIndSmall))),2));
session_probe2 = squeeze(any(~isnan(Tdat(:,:,:,2)),2));

%repeat
Tdat = squeeze(odorsel(:,sp.OdorIndSmall,18:19,:));
is_repeat= any(any(Tdat(:,:,:,2)<(sp.p_cutoff/(sum(sp.OdorIndSmall)*2)),2),3);                
session_repeat = neuron_info(:,5)>2;    
is_repeat2= squeeze(any(Tdat(:,:,:,2)<(sp.p_cutoff/(sum(sp.OdorIndSmall))),2));
session_repeat2 = squeeze(any(~isnan(Tdat(:,:,:,2)),2));

% odor off    
odoroff1 = squeeze(odorsel(:,sp.OdorOffInd,1,:));
odoroff12 = odoroff1(:,:,1); odoroff12(odoroff1(:,:,2) >= ...
    (sp.p_cutoff/sum(sp.OdorOffInd))) = NaN;
[~,maxOffInd] = max(abs(odoroff12),[],2);
maxOff = arrayfun(@(x) odoroff12(x,maxOffInd(x),1),1:length(maxOffInd));
odoroff2 = squeeze(odorsel_all(:,6,:));
isOdorOff = any(odoroff1(:,:,2)<(sp.p_cutoff/5),2) & ...
    odoroff2(:,2) < sp.p_cutoff & ((maxOff'>0) == (odoroff2(:,1)>0)) & ...
    ((maxOff'<0) == (odoroff2(:,1)<0));
isOdorOff_Up = any(odoroff1(:,:,2)<(sp.p_cutoff/5),2) & ...
    odoroff2(:,2) < sp.p_cutoff & maxOff'>0 & odoroff2(:,1)>0;

%odor off target/nontarget
Tdat = squeeze(odorsel(:,sp.OdorOffInd,3,:));
is_target_off= any(Tdat(:,:,2)<(sp.p_cutoff/sum(sp.OdorOffInd)) & ...
    Tdat(:,:,1)>0,2);
is_nontarget_off= any(Tdat(:,:,2)<(sp.p_cutoff/sum(sp.OdorOffInd)) & ...
    Tdat(:,:,1)<0,2);
look = find(is_target_off & is_nontarget_off);
[~,mm] = min(Tdat(look,:,2),[],2);
maxT = arrayfun(@(x) Tdat(look(x),mm(x),1),1:length(look));
is_target_off(look(maxT<0)) = false;
is_nontarget_off(look(maxT>0)) = false;

all_Neurons = cat(2,isSigOdor, islicksig, isUp, isDown, ...
    is_target, is_nontarget, is_correct, is_error, is_probe, ...
    session_probe, is_repeat, session_repeat, isOdorOff, isOdorOff_Up,...
    is_target_off, is_nontarget_off);   

sp.ana_labs = {'AllNeurons';'BestSession';'MeanSessions'};    
%all_Neurons, all_BestSes, all_MeanSes
sp.medcut = 9;    
sp.all_labs = {'SigOdor';'LickMod';'OdorOn-Up';...
    'OdorOn-Down';'TargetPref';'NonTargetPref'; ...
    'CorrectPref';'ErrorPref';'ProbeSelective';...
    'ProbeInclude';'RepeatSelective';'RepeatInclude'; ...
    'OdorOff';'OdorOff-Up';'OdorOff-Target';'OdorOff-NonTarget'};    
         %odorsel = odor/lick/target/correct/probex3/repeatx2
sp.rep_probe = {'Probe';'Repeat'}; 
sp.MM = {'Mean';'Median'};
sp.odorwindowlabel = odorwindowlabel;
selectivity_params = sp;

selectivity_dat.all_Neurons = all_Neurons;
selectivity_dat.all_MeanSes = all_MeanSes;
selectivity_dat.all_BestSes = all_BestSes;
selectivity_dat.is_issig_lick= is_issig_lick;
selectivity_dat.neuron_info = neuron_info;
selectivity_dat.is_probe2 = is_probe2;
selectivity_dat.session_probe2= session_probe2;
selectivity_dat.is_repeat2 = is_repeat2;
selectivity_dat.session_repeat2 = session_repeat2;
selectivity_dat.islickdown = islickdown;
selectivity_dat.islickup = islickup;
selectivity_dat.is_probeMean2 = is_probeMean2;
selectivity_dat.session_probeMean2 = session_probeMean2;
selectivity_dat.is_repeatMean = is_repeatMean;
selectivity_dat.session_repeatMean = session_repeatMean;
selectivity_dat.is_repeatMean2 = is_repeatMean2;
selectivity_dat.session_repeatMean2 = session_repeatMean2;
selectivity_dat.is_targetTrialCorrectedMean = is_targetTrialCorrectedMean;
