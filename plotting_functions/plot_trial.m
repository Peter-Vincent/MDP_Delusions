% 28/06/2020 Peter Vincent

function plot_trial(mdp)

[decision_traj,feedback_traj,card_traj] = behaviour(mdp);
[advised_card,trustworthy] = advisor(mdp);

num_trials = length(mdp);
figure('color','w')
max_pos = length(mdp);
max_total = round(max_pos * 1.5);
hold on
xlim([0.5 (num_trials+0.5)]);
ylim([-max_total max_total])
sep = 10;
for trial = 1:num_trials
    
    cur_bounds = decision_traj(1,trial);
    %% Plot the behavior of the agent
    if feedback_traj(1,trial) == 1
        feedback_colour = 'b';
    else
        feedback_colour = 'g';
    end
    rectangle('Position',[trial-0.5,-max_pos,1,cur_bounds+max_pos],...
        'FaceColor',feedback_colour,'EdgeColor','none');
    if card_traj(1,trial) == 1
        card_colour = 'r';
    else
        card_colour = 'k';
    end
    rectangle('Position',[trial-0.5,cur_bounds,1,max_pos-cur_bounds],...
        'FaceColor',card_colour,'EdgeColor','none');
    %% Plot the behaviour of the advisor
    if trustworthy(1,trial) == 1
        trustworthy_colour = 'b';
    else
        trustworthy_colour = 'g';
    end
    rectangle('Position',[trial-0.5,-max_total,1,max_total-max_pos-sep],...
        'FaceColor',trustworthy_colour,'EdgeColor','none');
    if advised_card(1,trial) == 1
        card_colour = 'r';
    else
        card_colour = 'k';
    end
    rectangle('Position',[trial-0.5,max_pos + sep,1,max_total-max_pos-sep],...
        'FaceColor',card_colour,'EdgeColor','none');
    
end

plot(1:num_trials,decision_traj,'LineWidth',2,'Color','k');
ax = gca;
set(ax,'Box','off');
ax.YAxis.Visible = 'off';    
