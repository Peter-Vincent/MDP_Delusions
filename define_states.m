% 14/06/2020 Peter Vincent

% function to generate a predetermined set of states

function state_matrix =  define_states(num_trials,card,advisor)
% Assumes there are 5 states, and will set the states that aren't card or
% advisor as some default state

% num_trials = no trials at start in which events are fixed
% card: 1 = A or 2 = B
% advisor: 1 = trustworthy, 2 = untrustworthy
% NB card and advisor should be scalars (repeated for num_trials) or vector
% sequences of cards/trustworthiness of length(num_trials)

state_matrix = zeros(5,num_trials);
for trial = 1:num_trials
    if length(advisor) == 1
        state_matrix(1,trial) = advisor;
    elseif length(advisor) == num_trials
        state_matrix(1,:) = advisor;
    end
    if length(card) == 1
        state_matrix(2,trial) = card;
    elseif length(advisor) == num_trials
        state_matrix(2,:) = card;
    end
    state_matrix(3,trial) = 3;
    state_matrix(4,trial) = 3;
    state_matrix(5,trial) = 1;
end
