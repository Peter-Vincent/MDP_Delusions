% 01/08/2020 Peter Vincent

function vary_habit_harness(habit_range,num_repeats,a_conting,initial_states)

% A harness to run the batch vary habit, requiring the same arguments.
% This function processes the output and uses these to generate figures

[fin_perf,fin_inf_correct,fin_habits] = batch_vary_habit(habit_range,num_repeats,a_conting,initial_states);

mean_performance = mean(fin_perf,2);
std_performance = std(fin_perf,0,2);

mean_inf_cor = mean(fin_inf_correct,2);
std_inf_cor  = std(fin_inf_correct,0,2);

mean_habits  = squeeze(mean(fin_habits,2));
std_habits   = squeeze(std(fin_habits,0,2));

% Make performance figure
figure('color','w')
perf = gca;
hold(perf,'on')
plot(perf,habit_range,mean_performance,'color','k','LineWidth',4);
plot(perf,habit_range,mean_performance+std_performance,'color','r','LineWidth',4);
plot(perf,habit_range,mean_performance-std_performance,'color','r','LineWidth',4);
set(perf,"XDir","reverse")
plot(perf,habit_range,fin_perf,'--k')
perf.Title.String = "Performance across trials";
perf.XLabel.String = "Habit parameter";
perf.YLabel.String = "Average performance";

% Make inference figure
figure('color','w')
infer = gca;
hold(infer,'on')
plot(infer,habit_range,mean_inf_cor,'color','k','LineWidth',4);
plot(infer,habit_range,mean_inf_cor+std_inf_cor,'color','r','LineWidth',4);
plot(infer,habit_range,mean_inf_cor-std_inf_cor,'color','r','LineWidth',4);
set(infer,"XDir","reverse")
plot(infer,habit_range,fin_inf_correct,'--k')
infer.Title.String = "Accurate inference across trials";
infer.XLabel.String = "Habit parameter";
infer.YLabel.String = "Average performance";

% Make habits figure
figure('color','w')
habits = gca;
hold(habits,'on');
plot(habits,habit_range,mean_habits,'color','k','LineWidth',2);
plot(habits,habit_range,mean_habits+std_habits,'color','r','LineWidth',4);
plot(habits,habit_range,mean_habits-std_habits,'color','r','LineWidth',4);
set(habits,"XDir","reverse")
habits.Title.String = "Habits across trials";
habits.XLabel.String = "Habit parameter";
habits.YLabel.String = "Trajectory weights";