% 25/10/2020 Peter Vincent

function extract_data(precision,affect,habit_param,sequence_path,num_simulations,save_dir)

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
setup
randseeds = linspace(0,100*num_simulations,num_simulations+1);
randseeds(2:end);
sequence_dir = load(sequence_path);
num_trials = size(sequence_dir.sequence,2);
trust_advice   = zeros(num_simulations,num_trials);
feedback       = zeros(num_simulations,num_trials);
card_chosen    = zeros(num_simulations,num_trials);
advised_card   = zeros(num_simulations,num_trials);
advisor_state  = zeros(num_simulations,num_trials);
prior_advisor  = zeros(num_simulations,num_trials);
post_advisor   = zeros(num_simulations,num_trials);
autonomic_state= zeros(num_simulations,num_trials);
decision_traj  = zeros(num_simulations,num_trials);
delusion_ad    = zeros(num_simulations,num_trials);
trustworthy_val=2;
for cur_sim = 1:num_simulations
    cur_mdp = MDP_Delusions_harness(num_trials,sequence_dir.sequence,randseeds(cur_sim),...
         precision,habit_param,affect);
    %%% Now we can extract the information we need from the MDP
    for trial = 1:num_trials
       cur_states = cur_mdp(trial).s;
       cur_outcomes=cur_mdp(trial).o;
       inf_ad = cur_mdp(trial).xn{1,1};
       advisor_state(cur_sim,trial) = cur_states(1,2);
       advised_card(cur_sim,trial)= cur_outcomes(1,2);
       feedback(cur_sim,trial) = cur_outcomes(2,3);
       card_chosen(cur_sim,trial)     = cur_outcomes(4,3);
       if advised_card(cur_sim,trial) == card_chosen(cur_sim,trial)
           trust_advice(cur_sim,trial) = 1;
           decision_traj(cur_sim,trial:end) = decision_traj(cur_sim,trial) + 1;
       else
           decision_traj(cur_sim,trial:end) = decision_traj(cur_sim,trial) - 1;
       end
       prior_advisor(cur_sim,trial) = inf_ad(end,1,2,2);
       post_advisor(cur_sim,trial)  = inf_ad(end,1,2,3);
       if (cur_outcomes(1,2) == cur_outcomes(4,3)) && (cur_outcomes(2,3) == 1)
           trustworthy_val = 1;
       elseif (cur_outcomes(1,2) ~= cur_outcomes(4,3)) && (cur_outcomes(2,3) == 1)
           trustworthy_val = 0;
       elseif (cur_outcomes(1,2) == cur_outcomes(4,3)) && (cur_outcomes(2,3) == 2)
           trustworthy_val = 0;
       elseif (cur_outcomes(1,2) ~= cur_outcomes(4,3)) && (cur_outcomes(2,3) == 2)
           trustworthy_val = 1;
       end
       if (trustworthy_val == 1) &&  (post_advisor(cur_sim,trial) > 0.5)
           delusion_ad(cur_sim,trial) = 0;
       elseif (trustworthy_val == 0) && (post_advisor(cur_sim,trial) < 0.5)
           delusion_ad(cur_sim,trial) = 0;
       else
           delusion_ad(cur_sim,trial) = 1;
       end
       autonomic_state(cur_sim,trial) = cur_states(4,2);
    end
end
break_down = strsplit(sequence_path,filesep);
file_title = strcat("BatchSim_",num2str(precision),"_",num2str(affect),"_",num2str(habit_param),"_",break_down(end));
save_path  = fullfile(save_dir,file_title);
save(save_path,"trust_advice","feedback","card_chosen","advised_card","advisor_state","prior_advisor","post_advisor","autonomic_state","decision_traj","delusion_ad");



