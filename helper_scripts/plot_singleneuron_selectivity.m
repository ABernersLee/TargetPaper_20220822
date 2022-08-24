function plot_singleneuron_selectivity(Ls,ineuron,yL,raster_data,dirs,...
    odorsel,sm,toplotOdor,figlablab,mousenames,labeladd,...
    odorwindowlabel,toplotLick,neuron_info,neuronnum,thelicktime,...
    reallabelsTarget,cor,odorsel_all,isProbe,isRepeat,trialtypes)

%notes:
%odorsel_all(ineuron,:,1) = [Target_r Correct_r Target_r_off ...
% Correct_r_off mod1_overall mod1_off];

%odorsel_all(ineuron,:,2) = [Target_p Correct_p Target_p_off ...
% Correct_p_off p_overall p_overall_off];

p_cutoff = 0.05;

if ~isfolder([dirs.figdir '/Neuron_examples/'])
    mkdir([dirs.figdir '/Neuron_examples/'])
end

[x1,y1] = find(raster_data);  
[li,ord] = sort(thelicktime);
li = li+5000;

x = NaN(size(x1));
for ix = 1:max(x1)
    x(x1==ord(ix)) = ix;            
end
y = y1; xsave = x;

for ifig = 1:7
    
    if ifig==1                
        pvals = odorsel(ineuron,:,1,2);                
        pvalsLick = odorsel(ineuron,:,2,2);
        moddir = odorsel(ineuron,:,1,1);
        moddirLick = odorsel(ineuron,:,2,1);
    elseif ifig==2                
        pvals = odorsel(ineuron,:,3,2);
        moddir = odorsel(ineuron,:,3,1);
    elseif ifig>=3 
        pvals = odorsel(ineuron,:,ifig+2,2);            
        moddir = odorsel(ineuron,:,ifig+2,2);
    end
    
    if (ifig>=3 && ifig<=5 && ~isProbe) || (ifig>=6 && ~isRepeat) 
        disp(['Skip early ' num2str(ifig) ' ' num2str(ineuron)])
        continue
    end
    
    if ifig == 2
        
        reallabels2 = reallabelsTarget(ord); 
        cor1 = cor(ord);
        [rlabs,ord2] = sort(reallabels2);
                        
        y = y1;                
        x = NaN(size(x1));
        xsave2 = xsave; 
        li1 = li;
        if size(cor1,1)~=size(reallabels2,1)
            cor1=cor1';
        end

        li1((isnan(reallabels2) | cor1~=1)') = NaN;       
        y(ismember(xsave2,find((isnan(reallabels2) | cor1~=1)')))...
            = NaN;         
        xsave2(ismember(xsave2,find((isnan(reallabels2) | cor1~=1)')))...
            = NaN;
        
        for ix = 1:max(x1)
            x(xsave2==ord2(ix)) = ix;            
        end               
        li2 = li1(ord2);  
        trialtypesFig = reallabelsTarget;
    end
    
    if ifig > 2 && (isProbe || isRepeat)
        trialtypesFig = NaN(size(trialtypes))';                
        trialtypesFig(trialtypes==2) = 2;
        
        if ifig>=3 && ifig<=5 && isProbe
            trialtypesFig(trialtypes==Ls(1,ifig-2)) = 1;
        elseif ifig>=6 && isRepeat
            trialtypesFig(trialtypes==Ls(2,ifig-5)) = 1;
        end
        
        reallabels2 = trialtypesFig(ord);
        cor1 = cor(ord);
        [rlabs,ord2] = sort(reallabels2);
                        
        y = y1;                
        x = NaN(size(x1));
        xsave2 = xsave; 
        li1 = li;
        if size(cor1,1)~=size(reallabels2,1)
            cor1=cor1';
        end
        li1((isnan(reallabels2) | cor1~=1)') = NaN;       
        y(ismember(xsave2,find((isnan(reallabels2) | cor1~=1)')))...
            = NaN;         
        xsave2(ismember(xsave2,find((isnan(reallabels2) | cor1~=1)')))...
            = NaN;
        for ix = 1:max(x1)
            x(xsave2==ord2(ix)) = ix;            
        end               
        li2 = li1(ord2);  
         
        
        if length(unique(trialtypesFig(~isnan(trialtypesFig))))<2                
            disp(['Skip late ' num2str(ifig) ' ' num2str(ineuron)])
            continue
        end
    end
    

    figure; hold on

    if ifig==1
        set(gcf,'Position',[ 1923         222        1074         869])
        if length(y)<15000

            subplot(2,3,1); hold on
            plot([y y]',[x-.5 x+.5]','k','LineWidth',3)            
            plot([li li]',[(1:length(li))-.5; (1:length(li))+.5]...
                ,'c','LineWidth',3)
            plot([5000 5000],[min(x)-.5 max(x)-.5],'b-')            
            plot([7500 7500],[min(x)-.5 max(x)-.5],'g-')            
            plot([7000 7000],[min(x)-.5 max(x)-.5],'r-')
            axis tight
            xt = get(gca,'xtick'); set(gca,'xticklabels',(xt./1000)-5)
            xlabel('Time since odor onset (sec)')
            ylabel('Trial number')
        end


        subplot(2,3,4); hold on
        
        dat = raster_data; 
        tbs = 1:size(dat,2);
        
        m1 = mean(dat); sem1 = std(dat)./sqrt(size(dat,1));               
        m = smoothdata(m1,'movmean',sm)*1000; 
        
        sem = smoothdata(sem1,'movmean',sm)*1000; 
        
        rev = m+sem; 
        fwd = m-sem;
        
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],...
            'black');
        p1.FaceAlpha=.2; 
        p1.EdgeAlpha=0;    
        p1.FaceColor = [0 0 0];
        
        plot(tbs,m,'k')
        
        axis tight                        
        yl = get(gca,'ylim');      
        
        plot([5000 5000],[min(yl) max(yl)],'b-')            
        plot([7500 7500],[min(yl) max(yl)],'g-')            
        plot([7000 7000],[min(yl) max(yl)],'r-')   
        
        xt = get(gca,'xtick'); set(gca,'xticklabels',(xt./1000)-5)
        xlabel('Time since odor onset (sec)')
        ylabel('Firing rate (Hz)')
    else
        set(gcf,'Position',[2090 438 523 684])
    end

    xx = x(y>=4000 & y<=10000);        
    yy = y(y>=4000 & y<=10000);

    if ifig>1
        subplot(2,1,1); hold on
        yc = find(rlabs==2,1,'first');
        plot([yy(xx<yc) yy(xx<yc)]',[xx(xx<yc)-.5 xx(xx<yc)+.5]',...
            'b','LineWidth',3)  
        plot([yy(xx>=yc) yy(xx>=yc)]',[xx(xx>=yc)-.5 xx(xx>=yc)+.5]',...
            'k','LineWidth',3)                  
        plot([li2 li2]',[(1:length(li2))-.5; (1:length(li2))+.5],'c',...
            'LineWidth',3)
    else
        subplot(2,3,2); hold on
        plot([yy yy]',[xx-.5 xx+.5]','k','LineWidth',3)                   
        plot([li li]',[(1:length(li))-.5; (1:length(li))+.5],...
            'c','LineWidth',3)
    end     

    axis tight 
    plot([5000 5000],[min(xx) max(xx)-.5],'r-')            
    plot([7500 7500],[min(xx) max(xx)-.5],'g-')            
    plot([7000 7000],[min(xx) max(xx)-.5],'r-')   
    xlim([4000 10000])
    set(gca,'xtick',1000*(4:10))
    xt = get(gca,'xtick'); 
    set(gca,'xticklabels',(xt./1000)-5)
    xlabel('Time since odor onset (sec)')
    ylabel('Trial number')


    if ifig>1
        subplot(2,1,2); hold on

        cor2 = cor; 
        tt2 = trialtypesFig';
        dat = toplotOdor(cor2==1 & tt2==1,:); 
        tbs = 1:size(dat,2);
        
        m1 = mean(dat); sem1 = std(dat)./sqrt(size(dat,1));               
        m = smoothdata(m1,'movmean',sm)*1000; 
        
        sem = smoothdata(sem1,'movmean',sm)*1000;
        
        rev = m+sem; 
        fwd = m-sem;

        p11 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'b');
        p11.FaceAlpha=.2;
        p11.EdgeAlpha=0;
        p11.FaceColor = [0 0 1];
        
        plot(tbs,m,'b')
        
        dat = toplotOdor(cor2==1 & tt2==2,:);
        m1 = mean(dat); 
        sem1 = std(dat)./sqrt(size(dat,1));             
        
        m = smoothdata(m1,'movmean',sm)*1000; 
        sem = smoothdata(sem1,'movmean',sm)*1000;
        
        rev = m+sem;
        fwd = m-sem;

        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],...
            'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(tbs,m,'k')
        axis tight
    else
        aa = subplot(2,3,5); hold on
        
        dat = toplotOdor; 
        tbs = 1:size(dat,2);
        
        m1 = mean(dat); 
        sem1 = std(dat)./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sm)*1000; 
        sem = smoothdata(sem1,'movmean',sm)*1000;
        rev = m+sem; fwd = m-sem;

        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],...
            'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];            
        plot(tbs,m,'k')
        
        set(gca,'xtick',1000*(0:6))
        xt = get(gca,'xtick'); set(gca,'xticklabels',(xt./1000)-1)
              
    end

       
    OdorIndSmall = (odorwindowlabel(:,1)>=0 & odorwindowlabel(:,2)<=1000)';
    LickInd = (odorwindowlabel(:,1)>=-400 & odorwindowlabel(:,2)<=600)';    
    OdorOffInd = (odorwindowlabel(:,1)>=2000 & odorwindowlabel(:,2)<=3000)';
    
    
    yl = get(gca,'ylim');      
    tp = .05*range(yl);
    plot([5000 5000]-4000,[min(yl) max(yl)],'b-')            
    plot([7500 7500]-4000,[min(yl) max(yl)],'g-')            
    plot([7000 7000]-4000,[min(yl) max(yl)],'r-')                                     

    plot(odorwindowlabel(OdorIndSmall & ...
        pvals>=(p_cutoff/sum(OdorIndSmall)),:)'+5000-4000,...
        (ones(sum(OdorIndSmall & ...
        pvals>=(p_cutoff/sum(OdorIndSmall))),1)*[yl(2)+tp yl(2)+tp])',...
        '.-','color',[.5 .5 .5],'LineWidth',2,'MarkerSize',15);

    plot(odorwindowlabel(OdorIndSmall & ...
        pvals<(p_cutoff/sum(OdorIndSmall)) & moddir<0,:)'+5000-4000,...
        (ones(sum(OdorIndSmall & ...
        pvals<(p_cutoff/sum(OdorIndSmall)) & moddir<0),1)*...
        [yl(2)+tp*2 yl(2)+tp*2])','b.-','LineWidth',2,'MarkerSize',15);  

    plot(odorwindowlabel(OdorIndSmall & ...
        pvals<(p_cutoff/sum(OdorIndSmall)) & moddir>0,:)'+5000-4000,...
        (ones(sum(OdorIndSmall & ...
        pvals<(p_cutoff/sum(OdorIndSmall))  & moddir>0),1)*...
        [yl(2)+tp*3 yl(2)+tp*3])','r.-','LineWidth',2,'MarkerSize',15); 
   
    if ifig==1

        plot(odorwindowlabel(OdorOffInd & ...
            pvals>=(p_cutoff/sum(OdorOffInd)),:)'+5000-4000,...
            (ones(sum(OdorOffInd & ...
            pvals>=(p_cutoff/sum(OdorOffInd))),1)*[yl(2)+tp yl(2)+tp])',...
            '.-','color',[.5 .5 .5],'LineWidth',2,'MarkerSize',15);

        plot(odorwindowlabel(OdorOffInd & ...
            pvals<(p_cutoff/sum(OdorOffInd)) & moddir<0,:)'+5000-4000,...
            (ones(sum(OdorOffInd & ...
            pvals<(p_cutoff/sum(OdorOffInd))  & moddir<0),1)*...
            [yl(2)+tp*2 yl(2)+tp*2])','b.-','LineWidth',2,'MarkerSize',15);          

        plot(odorwindowlabel(OdorOffInd & ...
            pvals<(p_cutoff/sum(OdorOffInd)) & moddir>0,:)'+5000-4000,...
            (ones(sum(OdorOffInd & ...
            pvals<(p_cutoff/sum(OdorOffInd))  & moddir>0),1)*...
            [yl(2)+tp*3 yl(2)+tp*3])','r.-','LineWidth',2,'MarkerSize',15); 
    end

    %overall_p
    if ifig==2 && odorsel_all(ineuron,5,2)<0.05 
        plot([500 1000],[yl(2)+tp*4 yl(2)+tp*4],'k-')
        plot([1100 1600],[yl(2)+tp*4 yl(2)+tp*4],'k*-')
    end

    %Target_p
    if ifig==2 && odorsel_all(ineuron,1,2)<0.05 
        plot([1100 1600],[yl(2)+tp*5 yl(2)+tp*5],'k*-')
    end

    %overall_p_off
    if ifig==2 && odorsel_all(ineuron,6,2)<0.05
        plot([2500 3000],[yl(2)+tp*4 yl(2)+tp*4],'k-')
        plot([3100 3600],[yl(2)+tp*4 yl(2)+tp*4],'k*-')
    end

    %Target_p_off
    if ifig==2 && odorsel_all(ineuron,3,2)<.05 
        plot([3100 3600],[yl(2)+tp*5 yl(2)+tp*5],'k*-')
    end

    if ifig==2
        title(['OdorOn: Mod = ' ...
            num2str(round(odorsel_all(ineuron,1,1),2,'significant')) ...
            ', p = ' num2str(round(odorsel_all(ineuron,1,2),2,...
            'significant')) ', OdorOff: R = ' ...
            num2str(round(odorsel_all(ineuron,3,1),2,'significant')) ...
            ', p = ' num2str(round(odorsel_all(ineuron,3,2),2,...
            'significant'))])
    end

    xlabel('Time since odor onset (sec)')
    ylabel('Firing rate (Hz)')        
    
    ylaa = get(gca,'ylim');
    ylaa(1) = 0;
    
    xt = get(gca,'xtick'); 
    set(gca,'xticklabels',((xt+4000)./1000)-5)   

    if ifig==1
        
        subplot(2,3,3); hold on
        plot([-li -li]'+5000,[(1:length(li))-.5; (1:length(li))+.5]...
            ,'b','LineWidth',3)        
        plot([-li -li]'+5000+2000,[(1:length(li))-.5; (1:length(li))+.5]...
            ,'r','LineWidth',3)        
        plot([-li -li]'+5000+2500,[(1:length(li))-.5; (1:length(li))+.5]...
            ,'g','LineWidth',3)
        plot([yL yL]',[xsave-.5 xsave+.5]','k','LineWidth',3)

        axis tight
        xlim([-2000 4000])         
        plot([0 0],[min(xsave) max(xsave)-.5],'c-')                     
        set(gca,'xtick',1000*(-2:4))        
        xt = get(gca,'xtick'); 
        set(gca,'xticklabels',(xt./1000))

        xlabel('Time since first lick (sec)')
        ylabel('Trial number')

        bb =  subplot(2,3,6); hold on
        dat = toplotLick; 
        tbs = 1:size(dat,2);

        m1 = mean(dat);
        sem1 = std(dat)./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sm)*1000;
        sem = smoothdata(sem1,'movmean',sm)*1000; 
        
        rev = m+sem; fwd = m-sem;
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],...
            'black');
        p1.FaceAlpha=.2;
        p1.EdgeAlpha=0;
        p1.FaceColor = [0 0 0];
        
        plot(tbs,m,'k')
        axis tight                        
        
        yl = get(gca,'ylim');      
        tp = .05*range(yl);
        yl = get(gca,'ylim');
        
        plot([6000 6000],[min(yl) max(yl)],'c-')                   
        xlim([4000 10000])
        set(gca,'xtick',1000*(4:10))
        xt = get(gca,'xtick'); set(gca,'xticklabels',((xt)./1000)-6)   

        plot(odorwindowlabel(LickInd & ...
            pvalsLick>=(p_cutoff/sum(LickInd)),:)'+6000,...
            (ones(sum(LickInd & pvalsLick>=(p_cutoff/sum(LickInd))),1)*...
            [yl(2)+tp yl(2)+tp])','.-','color',[.5 .5 .5],...
            'LineWidth',2,'MarkerSize',15);

        plot(odorwindowlabel(LickInd & ...
            pvalsLick<(p_cutoff/sum(LickInd)) & moddirLick<0,:)'+6000,...
            (ones(sum(LickInd & pvalsLick<(p_cutoff/sum(LickInd)) ...
            & moddirLick<0),1)*[yl(2)+tp*2 yl(2)+tp*2])','b.-',...
            'LineWidth',2,'MarkerSize',15);          
        
        plot(odorwindowlabel(LickInd & ...
            pvalsLick<(p_cutoff/sum(LickInd)) & moddirLick>0,:)'+6000,...
            (ones(sum(LickInd & pvalsLick<(p_cutoff/sum(LickInd)) ...
            & moddirLick>0),1)*[yl(2)+tp*3 yl(2)+tp*3])','r.-',...
            'LineWidth',2,'MarkerSize',15); 
   
        xlabel('Time since first lick (sec)')
        ylabel('Firing rate (Hz)')
        ylbb = get(gca,'ylim');

        set(aa,'ylim',[min([ylaa(1) ylbb(1)]) max([ylaa(2) ylbb(2)])])        
        set(bb,'ylim',[min([ylaa(1) ylbb(1)]) max([ylaa(2) ylbb(2)])])
    end


     if ifig==1

        suptitle(['mouse ' mousenames{neuron_info(ineuron,2)} ...
            ', day ' num2str(neuron_info(ineuron,3)) ', neuron ' ...
            num2str(ineuron) ', sessionneuron ' num2str(neuronnum)])                    

            label = [dirs.figdir '/Neuron_examples/neuron_' ...
                num2str(ineuron) '_mouse' ...
                mousenames{neuron_info(ineuron,2)} '_day' ...
                num2str(neuron_info(ineuron,3)) '_sessionneuron_' ...
                num2str(neuronnum) '_odorselective_' labeladd];            
     else

        suptitle([figlablab{ifig-1} '/Non: mouse ' ...
            mousenames{neuron_info(ineuron,2)} ', day ' ...
            num2str(neuron_info(ineuron,3)) ', neuron ' ...
            num2str(ineuron) ', sessionneuron ' num2str(neuronnum)])
         
        label = [dirs.figdir '/Neuron_examples/neuron_' ...
            num2str(ineuron) '_mouse' mousenames{neuron_info(ineuron,2)}...
            '_day' num2str(neuron_info(ineuron,3)) '_sessionneuron_'...
            num2str(neuronnum) '_' figlablab{ifig-1} 'Modulated_' labeladd];                                                                               
     end
     
     saveas(gcf,[label '.tif'],'tiff') 
     set(gcf,'renderer','Painters')
     saveas(gcf,[label '.pdf'],'pdf')
     saveas(gcf,[label '.eps'],'eps')             
     close gcf
end