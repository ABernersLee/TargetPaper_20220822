function [decoded_fold_group, decoded_fold, decoded_fold_cat, ...
            decoded_fold_S, decoded_fold_cat_S, decoded_fold_Sprob, ...
            decoded_fold_cat_Sprob] = get_decoding_seperate_trialtypes...
        (ns,trialdatsession,isession,all_labels,itype)


if ~all_labels.todownsample
    all_labels.numfolds2 = 1;
end
%Decoding the trial type with the population        
catgroup = categorical(trialdatsession(:,1));             

%seperating into the 2-3 repeat/probe types to decode seperately
trialty = unique(catgroup);
thesetrials = trialty(3:end);

%Temp variables for both identity decoding and catagorical decoding
decoded_fold_group = ...
    NaN(size(ns,1),3,all_labels.numfolds2,length(thesetrials));
decoded_fold =  ...
    NaN(size(ns,1),3,all_labels.numfolds2,length(thesetrials)); 
decoded_fold_cat = decoded_fold;           
decoded_fold_S =  ...
    NaN(size(ns,1),all_labels.numshuff,all_labels.numfolds2,...
    length(thesetrials)); 
decoded_fold_cat_S = decoded_fold_S; 
decoded_fold_cat_Sprob = decoded_fold_S; 
decoded_fold_Sprob = decoded_fold_S;

tic

for itrialtype = 1:length(thesetrials)
    
    ttt = [trialty(1:2); thesetrials(itrialtype)];
    
    %%%%%% to downsample
    if all_labels.todownsample
        hh = hist(catgroup(ismember(catgroup,ttt)));
        [minh,minw] = min(hh(hh>0)); 

        if minw~=3
            error('Target or nontarget less than expected')
        end     
    else
        minh = NaN;
    end
           
    
    disp(['Min ' num2str(minh) ...
        ', Neuron #: ' num2str(size(ns,2))...
        ', session ' num2str(isession) ', ' all_labels.typelab{itype}])

    if minh<all_labels.minhmin          
        disp('Skipping')
        continue
    end
    
    minh = all_labels.minhmin;      
    
    for ifoldrun = 1:all_labels.numfolds2

        if all_labels.todownsample
            %%%%%%to downsample
            tousenow = false(size(catgroup,1),1);
            tousenow(catgroup==ttt(3)) = 1;                                                
            for itt = 1:2
                jnk = find(catgroup==ttt(itt));
                jnk2 = jnk(randperm(size(jnk,1)));
                tousenow(jnk2(1:minh)) = 1;
            end    
        elseif ~all_labels.todownsample
            %%%%%% to not downsample   
            tousenow = ismember(catgroup,ttt); 
        end
        
        %%%% all
        ns22 = ns(tousenow,:);                
        ns22(:,sum(ns22)<5) = []; 

        %zscore across those trials         
        ns2 = zscore(ns22,1);            
        
        
        catgroup2= categorical(trialdatsession(tousenow,1));
        catgroup_cat2 = categorical(trialdatsession(tousenow,1)==1);                                                 
        ind2 = find(tousenow); 
        
        jnk = repmat(1:all_labels.numfolds,[1 ceil(size(ns2,1)/all_labels.numfolds)]);        
        numfoldschunks = jnk(randperm(size(jnk,2)))';
        numfoldschunks = numfoldschunks(1:size(ns2,1));   
        for ifold = 1:all_labels.numfolds
               X = ns2(numfoldschunks==ifold,:); 
               Y = ns2(numfoldschunks~=ifold,:); 
               group = catgroup2(numfoldschunks~=ifold); 
               group_cat = catgroup_cat2(numfoldschunks~=ifold);

               A = arrayfun(@(x) min(length(unique(Y(:,x)))),1:size(Y,2));     
               X(:,A<all_labels.numexlund) = []; Y(:,A<all_labels.numexlund) = [];

               [R1_c,~,post] = classify(X,Y,group); 
               [R1_cat,~,post_cat] = classify(X,Y,group_cat); 
               if size(post,2)==1 || size(post_cat,2)==1
                   disp('problem')
               end

               %%% getting the correct choice index in order to
               %%% find the probability(correct choice)
               ind1 = find(numfoldschunks==ifold);               
               postind = arrayfun(@(x) ...
                   find(ismember(unique(group), ...
                   catgroup2(ind1(x)))), 1:length(ind1));      
               post2 = arrayfun(@(x) ...
                   post(x,postind(x)),1:length(ind1))'; 
               
               %%% adding to the matrix
               decoded_fold(ind2(numfoldschunks==ifold),...
                   :,ifoldrun,itrialtype) = ...
                   [double(catgroup2(numfoldschunks==ifold,1)) ...
                   catgroup2(numfoldschunks==ifold,1)==R1_c post2];
               decoded_fold_cat(ind2(numfoldschunks==ifold),...
                   :,ifoldrun,itrialtype) = ...
                   [double(catgroup2(numfoldschunks==ifold,1)) ...
                   catgroup_cat2(numfoldschunks==ifold,1)==R1_cat...
                   post_cat(:,1)]; 
               % probability that it is nontarget
               
               
               %%% saving more detailed info
               jnk = (1:size(post,2):size(post,2)*size(post,1))-1;
               post3 = post'; post3(postind+jnk) = NaN;
               decoded_fold_group(ind2(numfoldschunks==ifold),...
                   :,ifoldrun,itrialtype) = [post(:,1:2) ...
                   sum(post3(3:end,:),1,'omitnan')'];
         end

        %Running a shuffle of the trial types to decode to test if we
        %can decode better than chance
        if all_labels.runshuffle
            randd = rand(size(catgroup2,1),all_labels.numshuff);
            [~,randds] = sort(randd,1);
            catgroupS1 = repmat(catgroup2,[1 all_labels.numshuff]);
            catgroup_catS1 = repmat(catgroup_cat2,[1 all_labels.numshuff]);
            catgroupS = catgroupS1(randds);
            catgroup_catS = catgroup_catS1(randds);

           for ifold = 1:all_labels.numfolds                   
               X = ns2(numfoldschunks==ifold,:); Y = ns2(numfoldschunks~=ifold,:); 
               
               A = arrayfun(@(x) min(length(unique(Y(:,x)))),1:size(Y,2));     
               X(:,A<all_labels.numexlund) = []; Y(:,A<all_labels.numexlund) = [];

               if sum(numfoldschunks~=ifold)==1

                   error(['If im going to use this I '...
                       'need to write the code for probability here'])

               elseif sum(numfoldschunks~=ifold)>1    

                   [R1_c,~,postS] = arrayfun(@(x) ...
                       classify(X,Y,...
                       catgroupS(numfoldschunks~=ifold,x)), ...
                       1:all_labels.numshuff,'UniformOutput',false);                            
                   [R1_cat,~,postSc] = arrayfun(@(x) ...
                       classify(X,Y,...
                       catgroup_catS(numfoldschunks~=ifold,x)), ...
                       1:all_labels.numshuff,'UniformOutput',false);
                   
                   %%% getting the correct choice index in order to
                   %%% find the probability(correct choice)                           
                   Ss = reshape(cell2mat(postS),[size(X,1) ...
                       size(postS{1},2) size(postS,2)]);
                   ind1 = find(numfoldschunks==ifold);  
                   post2 = NaN(length(ind1),all_labels.numshuff);
                   for iind1 = 1:length(ind1)
                       postind= cell2mat(arrayfun(@(x) ...
                           find(ismember(unique(group), ...
                           catgroupS(ind1(iind1),x))), ...
                           1:all_labels.numshuff,'UniformOutput',false));                                                           
                       post2(iind1,:) = arrayfun(@(x) ...
                           Ss(iind1,postind(x),x),1:all_labels.numshuff)'; 
                   end
               
                   %%% get the NT out (change to the correct one
                   %%% after)
                   Sc = reshape(cell2mat(postSc),...
                       [size(X,1) 2 size(postSc,2)]);
                   
                                              
                   %%% add the accuracy
                   R1_c2 = cellfun(@(x) x==...
                       catgroup2(numfoldschunks==ifold,1),R1_c,...
                       'UniformOutput',false);  
                   decoded_fold_S(ind2(numfoldschunks==ifold),...
                       :,ifoldrun,itrialtype) = cell2mat(R1_c2);                           
                   R1_cat2 = cellfun(@(x) x==...
                       catgroup_cat2(numfoldschunks==ifold,1),...
                       R1_cat,'UniformOutput',false);  
                   decoded_fold_cat_S(ind2(numfoldschunks==ifold),...
                       :,ifoldrun,itrialtype) = cell2mat(R1_cat2);     
                   
                   %%% add the probability correct
                   decoded_fold_Sprob(ind2(numfoldschunks==ifold)...
                       ,:,ifoldrun,itrialtype) = post2; 
                   decoded_fold_cat_Sprob(ind2(numfoldschunks==ifold)...
                       ,:,ifoldrun,itrialtype) = squeeze(Sc(:,1,:));       
               end
           end   
        end                 
    end        
end
