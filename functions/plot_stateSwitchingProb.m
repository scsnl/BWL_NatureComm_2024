function StatTransProb_Task= plot_stateSwitchingProb(TemEvol_SubjStat, taskList, group_idx, num, group_name, reorder)

%% Plot state transition probability of each task conditions
figure;
for task=1:num.Task
    subplot(1,3,task);
% 	imagesc(TemEvol_GroupStat.StatTransProb{1,task});
    StatTransProb=squeeze(nanmean(TemEvol_SubjStat.StatTransProb{1,task}(group_idx,:,:)));

    StatTransProb=StatTransProb(reorder, reorder);

    StatTransProb_Task{1,task}=StatTransProb;
    
    imagesc(StatTransProb);
    ylabel('state(t)');
    xlabel('state(t+1)');
    title(sprintf('State transition prob. of %s Condition of %s', taskList{task}, group_name));
    caxis([0 1]);
    textStrings = num2str(StatTransProb(:), '%0.2f');       % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    [x, y] = meshgrid(1:5);  % Create x and y coordinates for the strings
    hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                    'HorizontalAlignment', 'center');
    midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
    % Choose white or black for the text color of the strings so they can be easily seen over the background color
    textColors = repmat(StatTransProb(:) > midValue, 1, 3);  
    set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors

    set(gca, 'XTick', 1:5, ...                             % Change the axes tick marks
             'XTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...  %   and tick labels
             'YTick', 1:5, ...
             'YTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...
             'TickLength', [0 0]);
         
    mycolormap = customcolormap([0 0.5 1], {'#9d0142','#f66e45','#ffffbb'});
    colorbar('southoutside');
    colormap(mycolormap);
end
