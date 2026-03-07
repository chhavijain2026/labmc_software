% plot histogram
% Inputs:
% 1. data : column vector
% 2. nbins : no. of bins for histogram
% 3. axislabel : string containing ylabel
% 4. flag = if 1: mark mean, if 0: do not mark mean
% 5. splot = subplot parameters

function [maxbinc] = plothist(data,nbins,axislabel,flag,splot)

subplot(splot(1),splot(2),splot(3))
%hist(data,nbins)
h=histogram(data,nbins); %,'Normalization','probability');
maxbinc=round(max(h.BinCounts)*10)/10;
% h=findobj(gca,'Type','patch');
set(h,'FaceColor',[0.75 0.75 0.75],'Edgecolor','k')
ylabel(axislabel,'FontSize',14)
if flag == 1
    md = mean(data); sd=std(data);
    hold on
    %plot([mean(data) mean(data)],[0 maxbinc],'r--')
    h=plot([md md],[0 maxbinc],'r--',[md-sd md-sd],[0 maxbinc],'r:',...
        [md+sd md+sd],[0 maxbinc],'r:');
    set(h,'LineWidth',1.5)
    hold off
end
return;
