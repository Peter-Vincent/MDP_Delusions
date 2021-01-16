% 25/10/2020 Peter Vincent

function summary = extract_data(cur_mdp)

% This function goes through the MDP and extracts 
%       - Trust or not
%       - Feedback given
%       - Card Selected
%       - Advised Card
%       - Trustworthy state
%       - Prior belief over advisor
%       - Posterior belief over advisor
%       - Prior belief over card
%       - Posterior belief over card
%       - Autonomic State
%       - If the agent had a delusion about the advisor state
%       - The decision trajectory
% For each trial and for each repeat
num_trials     = length(cur_mdp);
trust_advice   = zeros(1,num_trials);
feedback       = zeros(1,num_trials);
card_chosen    = zeros(1,num_trials);
advised_card   = zeros(1,num_trials);
advisor_state  = zeros(1,num_trials);
prior_advisor  = zeros(1,num_trials);
post_advisor   = zeros(1,num_trials);
autonomic_state= zeros(1,num_trials);
decision_traj  = zeros(1,num_trials);
delusion_ad    = zeros(1,num_trials);
trustworthy_val=2;

%%% Now we can extract the information we need from the MDP
for trial = 1:num_trials
    cur_states = cur_mdp(trial).s;
    cur_outcomes=cur_mdp(trial).o;
    inf_ad = cur_mdp(trial).xn{1,1};
    advisor_state(1,trial) = cur_states(1,2);
    advised_card(1,trial)= cur_outcomes(1,2);
    feedback(1,trial) = cur_outcomes(2,3);
    card_chosen(1,trial)     = cur_outcomes(4,3);
    if advised_card(1,trial) == card_chosen(1,trial)
        trust_advice(1,trial) = 1;
        decision_traj(1,trial:end) = decision_traj(1,trial) + 1;
    else
        decision_traj(1,trial:end) = decision_traj(1,trial) - 1;
    end
    prior_advisor(1,trial) = inf_ad(end,1,2,2);
    post_advisor(1,trial)  = inf_ad(end,1,2,3);
    if (cur_outcomes(1,2) == cur_outcomes(4,3)) && (cur_outcomes(2,3) == 1)
        trustworthy_val = 1;
    elseif (cur_outcomes(1,2) ~= cur_outcomes(4,3)) && (cur_outcomes(2,3) == 1)
        trustworthy_val = 0;
    elseif (cur_outcomes(1,2) == cur_outcomes(4,3)) && (cur_outcomes(2,3) == 2)
        trustworthy_val = 0;
    elseif (cur_outcomes(1,2) ~= cur_outcomes(4,3)) && (cur_outcomes(2,3) == 2)
        trustworthy_val = 1;
    end
    if (trustworthy_val == 1) &&  (post_advisor(1,trial) > 0.5)
        delusion_ad(1,trial) = 0;
    elseif (trustworthy_val == 0) && (post_advisor(1,trial) < 0.5)
        delusion_ad(1,trial) = 0;
    else
        delusion_ad(1,trial) = 1;
    end
    autonomic_state(1,trial) = cur_states(4,2);
end
summary.trust_advice = trust_advice;
summary.feedback = feedback;
summary.card_chosen = card_chosen;
summary.advised_card = advised_card;
summary.advisor_state = advisor_state;
summary.prior_advisor = prior_advisor;
summary.post_advisor  = post_advisor;
summary.autonomic_state = autonomic_state;
summary.decision_traj = decision_traj;
summary.delusions_ad  = delusion_ad;

