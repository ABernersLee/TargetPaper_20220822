function trp = helper_text_rp(r,p,num_sigfig)

trp = ['r = ' num2str(round(r,num_sigfig,'significant')) ...
    ', p = ' num2str(round(p,num_sigfig,'significant'))];