% 28/06/2020 Peter Vincent

function [decision_traj,feedback_traj,card_traj] = behaviour(mdp)

% Takes in and MDP output tile, plots the deicision of the agent to trust
% the advisor, the feedback that was recieved and which card was
% selected

num_trials = length(mdp);
decision_traj = zeros(1,num_trials);
card_traj     = zeros(1,num_trials);
feedback_traj    = zeros(1,num_trials);

for trial = 1:num_trials
    cur_trial = mdp(trial);
    cur_outcomes=cur_trial.o;
    feedback_traj(1,trial) = cur_outcomes(2,3);
    card_traj(1,trial)     = cur_outcomes(4,3);
    if cur_outcomes(4,3) == cur_outcomes(1,2)
        decision_traj(1,trial:end) = decision_traj(1,trial) + 1;
    else
        decision_traj(1,trial:end) = decision_traj(1,trial) - 1;
    end
end

