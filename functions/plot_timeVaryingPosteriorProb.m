function plot_timeVaryingPosteriorProb(StatePP_GroupStat, StatePP_SubjStat, num, groupList, group_idx, targetGroup)

colorstrings = {'#515564','#C25B56','#BEB9B5', '#74828F', '#96C0CE' };
xpts = 0:length(StatePP_GroupStat.Avg.wholeTB{1,1})-1; 
xrange = [0, (length(StatePP_GroupStat.Avg.wholeTB{1,1})-1).*0.49]; 
yrange = [0, 0.8]; 
linewidth = 4; fontsize=15;

targetGroupIdx=group_idx;

PlotList={'LL', 'HL', 'DL'};
figure;
sgtitle(sprintf('Posterior probability changes of the States of %s',groupList{targetGroup}), 'FontSize', fontsize);
for task = 1:num.Task
    subplot(1, num.Task, task)
	for state = 1:num.State
        y=mean(StatePP_SubjStat.Avg.wholeTB{1,task}(targetGroupIdx,:,state)); % your mean vector;            

%     %------Add shade (standard deviation)--------%         
%             x = 1:numel(y);            
%             std_dev = std(StatePP_SubjStat.Avg.wholeTB{1,task}(:,:,state));
%             curve1 = y + std_dev;
%             curve2 = y - std_dev;
%             x2 = [x, fliplr(x)];
%             inBetween = [curve1, fliplr(curve2)];            
%             fill(x2*0.49, inBetween, 'g');
%     %--------------------------------------------%
        hold on;
        plot(xpts*0.49, y, 'Color', colorstrings{state}, 'LineWidth', linewidth);
	end
    ylim(yrange); xlim(xrange); xlabel('Time(s)', 'FontSize', 12); ylabel('posterior probability', 'FontSize', 12);
    title(sprintf('%s', PlotList{task}), 'FontSize', 15);
end