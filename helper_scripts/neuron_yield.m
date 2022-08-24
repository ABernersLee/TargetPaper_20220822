function neuron_yield
% looking at neuron yield across days as a measure of density of the layer
% in relation to Nagappan and Franks 2021
% 
% my thinking: you would expect to hit layer 3, then 2, then 1
% layer 2 is the densest layer, so if anything here we are going from 2 to
% 1
% this means you could expect more semilunar pyrmidal cells, which Nagappan
% and Franks found to be LESS selective. So this goes against the idea that
% moving tetrodes down is why you see more selectivity over time. 

numn = [];
for imouse = 1:4
    sessions = unique(neuron_info(neuron_info(:,2)==imouse,3));
    for isess = 1:length(sessions)
        numn = cat(1,numn,[sum(neuron_info(:,3)==sessions(isess) & neuron_info(:,2)==imouse) sessions(isess) isess imouse]);
    end
    [r1,p1] = corr(numn(numn(:,4)==imouse,2),numn(numn(:,4)==imouse,1));
    [r2,p2] = corr(numn(numn(:,4)==imouse,3),numn(numn(:,4)==imouse,1));
    [imouse r1 p1 r2 p2]
end

% either way, decreases across days

[p,tbl,stats] = anovan(numn(:,1),numn(:,[2 4]),'Display','off','varnames',{'Day','Mouse'},'continuous',1);
[p(1) stats.coeffs(2)]

[p,tbl,stats] = anovan(numn(:,1),numn(:,[3 4]),'Display','off','varnames',{'Day','Mouse'},'continuous',1);
[p(1) stats.coeffs(2)]

