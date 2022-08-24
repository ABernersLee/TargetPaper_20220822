function [p,r] = get_shuffle_p_single_neuron_diff(allmod,dat,numshuff,inds)

if size(dat,1)==1
    dat = dat';
elseif size(dat,2)==1
else
    error('wrong size dat')
end

%get random shuffles
[~,shuff] = sort(rand(size(dat,1),numshuff));

%data structure the same size
dat2 = repmat(dat,[1 numshuff]);

%make neural data modulation same size
allmod1 = repmat(allmod,[1 numshuff]); 
allmod2 = allmod1;

% real difference between two types (inds) of labels (dat)
r = mean(allmod(dat==inds(1)))-mean(allmod(dat==inds(2)));

% shuffled labels
allmod1(dat2(shuff)~=inds(1)) = NaN;  
allmod2(dat2(shuff)~=inds(2)) = NaN;

% shuffled differences
s = abs(mean(allmod1,'omitnan')-mean(allmod2,'omitnan'));


% statistical test of likelihood that real is larger than shuffled
% differences
p = (1+sum(s>=abs(r)))./(1+numshuff);