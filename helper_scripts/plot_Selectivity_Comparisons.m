function plot_Selectivity_Comparisons(singledat,selectivity_dat,...
    selectivity_params,dirs,figs_to_run) 

sp = selectivity_params;
sd = selectivity_dat;

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'Comparisons/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end
     
  %%%%% Comparisons between categories of neurons, between timing too
times = sp.odorwindowlabel(:,1)>=0 & sp.odorwindowlabel(:,2)<=2000;

for iana = [1 3]
   
    is_target = sd.all_Neurons(:,5);                
    is_nontarget = sd.all_Neurons(:,6);    
    isUp = sd.all_Neurons(:,3);            
    isDown = sd.all_Neurons(:,4);

            
    odorbins = squeeze(singledat.odorsel(:,times,1,:));
             
    % Chi-squared test between up and down modulated and modulated by 
    %   target at all
    % Observed data
    n1 = sum((is_target | is_nontarget) & isDown); N1 = sum(isDown);
    n2 = sum((is_target | is_nontarget) & isUp); N2 = sum(isUp);
    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2);
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    % Chi-square test, by hand
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    [~,pp0,PPstats] = chi2gof([1 2 3 4],'freq',observed,'expected',...
        expected,'ctrs',[1 2 3 4],'nparams',2);
    PPstats.Pvalue = pp0;

    D = [sum(~is_target & ~is_nontarget & isDown) ...
        sum(is_target & isDown) sum(is_nontarget & isDown)]./(sum(isDown));
    U = [sum(~is_target & ~is_nontarget & isUp) ...
        sum(is_target & isUp) sum(is_nontarget & isUp)]./(sum(isUp));



    if (sum(ismember(figs_to_run,[0 4]))>0 && iana==1) || ...
            (sum(ismember(figs_to_run,0))>0)
        
        figure; hold on
        bar([D; U],'stacked')
        n1 = sum((is_target | is_nontarget) & isUp);
        n2 = sum(isUp);
        nn1 = sum((is_target | is_nontarget) & isDown);
        nn2 = sum(isDown);
        p01 = 1-binocdf(n1-1,n2,sp.p_cutoff*2);
        p02 = 1-binocdf(nn1-1,nn2,sp.p_cutoff*2);
        set(gca,'xtick',1:2,'xticklabel',{'Down';'Up'})
        text(2.5,.2,{'Chi-squared test of Up vs Down'; ...
            'being modulated by target at all';...
            ['p = ' num2str(round(pp0,2,'significant'))]})
        text(2.5,.85,{'Binomial probability of Up';...
            'cells being modulated,';[num2str(n1) ' of ' num2str(n2) ...
            ', p = ' num2str(round(p01,2,'significant'))]})
        text(2.5,.55,{'Binomial probability of Down';...
            'cells being modulated,';[num2str(nn1) ' of ' ...
            num2str(nn2) ', p = ' num2str(round(p02,2,'significant'))]})        
        xlim([0 7])
        legend({'Not mod.','Target','NonTarget'},'Location','northeast')
        ylabel('Proportion of cells')        
        helper_saveandclosefig([thisdir 'PropSigUpDownTNT_pvalue' ...
            num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_'...
            sp.ana_labs{iana}])
        
        PP(1,:) = fieldnames(PPstats);
        PP(2,:) = struct2cell(PPstats);
        writecell(PP,[thisdir 'PropSigUpDownTNT_pvalue' ...
            num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) ...
            '_Seconds_' sp.ana_labs{iana} '_stats.xlsx'])
    end


end
