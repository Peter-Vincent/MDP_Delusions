% 28/06/2020 Peter Vincent

function plot_behaviour(mdp)

[decision_traj,feedback_traj,card_traj] = behaviour(mdp);

num_trials = length(mdp);
figure('color','w')

hold on
max_val = max(decision_traj)+5;
min_val = min(decision_traj)-5;
xlim([0.5 (num_trials+0.5)]);
ylim([min_val max_val])
for trial = 1:num_trials
    
    cur_bounds = decision_traj(1,trial);
    if feedback_traj(1,trial) == 1
        feedback_colour = 'b';
    else
        feedback_colour = 'g';
    end
    rectangle('Position',[trial-0.5,min_val,1,cur_bounds-min_val],...
        'FaceColor',feedback_colour,'EdgeColor','none');
    if card_traj(1,trial) == 1
        card_colour = 'r';
    else
        card_colour = 'k';
    end
    rectangle('Position',[trial-0.5,cur_bounds,1,max_val-cur_bounds],...
        'FaceColor',card_colour,'EdgeColor','none');
end

plot(1:num_trials,decision_traj,'LineWidth',2,'Color','k');
    
    
