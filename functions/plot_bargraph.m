function plot_bargraph(Value, group_idx, condition_names, group_names, Ylabel, Ylim, Fontsize, Title)

h = dabarplot(Value,'groups',group_idx,...
    'xtlabels', condition_names,'errorbars',0,...
    'scatter',1,'scattersize',15,'scatteralpha',0.5,...
    'barspacing',0.8,'legend',group_names); 
ylabel(Ylabel);
yl = ylim; ylim(Ylim);  % make more space for the legend
set(gca,'FontSize',Fontsize); title(Title);


