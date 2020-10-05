% 05/10/2020 Peter Vincent

function [advised_card,trustworthy] = advisor(mdp)
% This function parses the MDP to extract whether the advisor was
% trustworthy when the advice was given, and the card that was advised.

num_trials   = length(mdp);
advised_card = zeros(1,num_trials);
trustworthy  = zeros(1,num_trials);

for trial = 1:num_trials
    cur_trial    = mdp(trial);
    cur_states   = cur_trial.s;
    cur_outcomes = cur_trial.o;
    trustworthy(1,trial) = cur_states(1,2);
    advised_card(1,trial)= cur_outcomes(1,2);
end

