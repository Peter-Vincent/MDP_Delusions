% Peter Vincent
% 04/10/2020
% This function generates a bar chart showing the mean performance across
% different habit parameters

function mean_perf(habit_range,fin_perf)

num_params = length(habit_range);
if num_params ~= size(fin_perf,1)
    error("dimensions of input data do not match");
end
mean_perf = mean(fin_perf,2);
dev = std(fin_perf,0,2);
num_simulations = size(fin_perf,2);
figure('color','w')
ax = gca;
range = 1:num_params;
hold(ax,'on')

chart = bar(ax,range,mean_perf,'FaceColor','#1C1BFF','EdgeColor','k');
for i = 1:num_params
    scatter(ax,ones(1,num_simulations)*range(i),fin_perf(i,:),50,...
    'MarkerFaceColor','#8133FF','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',0.5);
end
err = errorbar(ax,range,mean_perf,dev,dev,'.','LineWidth',2,'Color','r');
hold(ax,'off')
set(ax,'XTickLabels',{'','1','5','50','100','300','600',''});
xlabel(ax,'Habit Parameter');
ylabel(ax,'Performance')