function thelicktime = get_licks(raster_labels)

L1 = raster_labels.all_licks_port_one; %right
LL1 = cell2mat(cellfun(@min,L1,'UniformOutput',false));    
fL1 = ~cellfun(@isempty, L1);
L2 = raster_labels.all_licks_port_two;  %left
LL2 = cell2mat(cellfun(@min,L2,'UniformOutput',false));    
fL2 = ~cellfun(@isempty, L2);
rightleft = NaN(size(L1,2),2);
rightleft(fL1,1) = LL1;
rightleft(fL2,2) = LL2;
thelicktime = NaN(size(rightleft,1),1);
thelicktime(strcmp(raster_labels.lick_choice,'right'),1)...
    =rightleft(strcmp(raster_labels.lick_choice,'right'),1);        
thelicktime(strcmp(raster_labels.lick_choice,'left'),1)=...
    rightleft(strcmp(raster_labels.lick_choice,'left'),2);