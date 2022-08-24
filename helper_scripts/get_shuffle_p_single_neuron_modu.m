function [p,mod1,postdat] = get_shuffle_p_single_neuron_modu(y,yL,pre,post,ysh,yLsh)

numshuff = size(ysh,2);

if size(pre,1)==1 && size(pre,2)==2
    dv = (abs(diff(pre))+1)/1000;

    if ((abs(diff(post))+1)/1000) ~= dv
        error('not equal times')
    end

    predat = (sum((y(~isnan(yL),:)>=pre(1) & ...
        y(~isnan(yL),:)<pre(2))))./dv; 
    postdat = (sum((y(~isnan(yL),:)>=post(1) & ...
        y(~isnan(yL),:)<post(2))))./dv;     
    preshuff = (sum((ysh(~isnan(yL),:)>=pre(1) & ...
        ysh(~isnan(yL),:)<pre(2))))./dv; 
    postshuff = (sum((ysh(~isnan(yL),:)>=post(1) & ...
        ysh(~isnan(yL),:)<post(2))))./dv; 

    mod1 = (postdat-predat)./(postdat+predat);
    mod1shuff = (postshuff-preshuff)./(postshuff+preshuff);

    p = (1+sum(abs(mod1shuff)>=abs(mod1)))./(1+numshuff);

elseif size(pre,1)==2 && size(pre,2)==2

    dv1 = ((abs(diff(post))+1)/1000);
    dv2 = sum((abs(diff(pre,[],2))))/1000;

    %baseline period
    predat = ...
        sum((y(~isnan(yL)) > pre(1,1) & y(~isnan(yL)) <= pre(1,2))...
        | (y(~isnan(yL)) > pre(2,1) & y(~isnan(yL)) <= pre(2,2)))...
        ./dv2;

    %baseline period shuffled
    preshuff = ...
        sum((ysh(~isnan(yL),:) > pre(1,1) & ysh(~isnan(yL),:) <= pre(1,2))...
        | (ysh(~isnan(yL),:) > pre(2,1) & ysh(~isnan(yL),:) <= pre(2,2)))...
        ./dv2;

    if sum(post>0)==0 && min(y)>0 % lick times instead
        postdat = sum((yL>=post(1) & yL<post(2)))./dv1;     
        postshuff = (sum((yLsh(~isnan(yL),:)>=post(1) & ...
            yLsh(~isnan(yL),:)<post(2))))./dv1; 
    else
        postdat = (sum((y(~isnan(yL),:)>=post(1) & ...
            y(~isnan(yL),:)<post(2))))./dv1;     
        postshuff = (sum((ysh(~isnan(yL),:)>=post(1) & ...
            ysh(~isnan(yL),:)<post(2))))./dv1; 
    end
   
    mod1 = (postdat-predat)./(postdat+predat);
    mod1shuff = (postshuff-preshuff)./(postshuff+preshuff);

    p = (1+sum(abs(mod1shuff)>=abs(mod1)))./(1+numshuff);

    %check
%     figure; histogram(abs(mod1shuff)); hold on; yl = get(gca,'ylim');
%     plot([abs(mod1) abs(mod1)],yl,'r-')
else
    error('wrong size predat')
end