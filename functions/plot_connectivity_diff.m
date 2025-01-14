function plot_connectivity_diff(target_conn, numROI, ROI_names, CAxis, task)

imagesc(target_conn);
ylabel('Connection');
xlabel('Connection');
title(sprintf('Dynamic functional connectivity in %s', task));
caxis(CAxis);
textStrings = num2str(target_conn(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:numROI);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
	'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
% Choose white or black for the text color of the strings so they can be easily seen over the background color
textColors = repmat(target_conn(:) > midValue, 1, 3);  
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors

set(gca, 'XTick', 1:numROI, ...                             % Change the axes tick marks
         'XTickLabel', ROI_names, ...  %   and tick labels
         'YTick', 1:numROI, ...
         'YTickLabel', ROI_names, ...
         'TickLength', [0 0]);
         
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
% mycolormap = customcolormap([0 0.5 1], {'#9d0142','#f66e45','#ffffbb'});
% mycolormap = customcolormap(linspace(0,1,5), {'#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3'});
colorbar('southoutside');
colormap(mycolormap);


