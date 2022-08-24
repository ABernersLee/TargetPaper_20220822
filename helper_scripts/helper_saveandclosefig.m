function helper_saveandclosefig(label)


saveas(gcf,[label '.fig'],'fig')
saveas(gcf,[label '.tif'],'tiff')


% illustrator_ready(gcf)

set(gcf,'renderer','Painters')
saveas(gcf,[label '.eps'],'epsc')
% printeps(1,[label '.eps'])
close gcf


%%% use again on other computer

% try
%     export_fig(gcf,label,'-transparent','-eps','-painters') %,'-font_space', '', '-regexprep', 'Arial')
%     
%     close gcf
%     
% catch
%     disp(['Problem, skipped'])
%     
%     set(gcf,'renderer','Painters')
%     saveas(gcf,[label '.eps'],'epsc')
%     close gcf
% end

