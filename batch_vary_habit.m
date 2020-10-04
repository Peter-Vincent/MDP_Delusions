% 31/07/2020 Peter Vincent

function [fin_perf,fin_inf_correct,fin_habits] = batch_vary_habit(habit_range,num_repeats,a_conting,initial_states)

% Function to run many iterations of the MDP_Delusions script.  Will store
% the % of correct inference, the final habits count, the % performance

num_sims = length(habit_range);
num_trials  = size(initial_states,2);
performance = zeros(num_sims,num_repeats);
inf_correct = zeros(num_sims,num_repeats);
habits      = zeros(num_sims,num_repeats,4); % There are 4 trajectories 

rand_seeds  = randi(1000000,num_sims,num_repeats);


parfor habit_val = 1:num_sims
    disp(strcat("Processing habit param ",num2str(habit_val)));
    cur_habit = habit_range(habit_val);
    for repeat = 1:num_repeats
        cur_seed = rand_seeds(habit_val,repeat);
        warning("off")
        cur_mdp  = MDP_Delusions(num_trials,initial_states,cur_seed,...
                                 a_conting,cur_habit);
        warning("on")
                             %%%%%% Now we parse the output to extract the information %%%%%%%
        final_habit = cur_mdp(end).e;
        norm_habit  = final_habit ./ sum(final_habit);
        habits(habit_val,repeat,:) = norm_habit;
        % Find the performance and inferences
        perform_vector = zeros(1,num_trials);
        correct_post_inf = zeros(1,num_trials);
        for trial = 1:num_trials
            
            outcomes = cur_mdp(trial).o;
            perform_vector(1,trial) = outcomes(2,end);
            if perform_vector(1,trial) == 2
                perform_vector(1,trial) = 0;
            end
            inf_over_ad=cur_mdp(trial).xn{1,1};
            post_feedback_inf = inf_over_ad(end,1,2,3);

            if (outcomes(1,2) == outcomes(4,3)) && (outcomes(2,3) == 1)
                trustworthy = 1;
            elseif (outcomes(1,2) ~= outcomes(4,3)) && (outcomes(2,3) == 1)
                trustworthy = 0;
            elseif (outcomes(1,2) == outcomes(4,3)) && (outcomes(2,3) == 2)
                trustworthy = 0;
            elseif (outcomes(1,2) ~= outcomes(4,3)) && (outcomes(2,3) == 2)
                trustworthy = 1;
            end
            if (trustworthy == 1) &&  (post_feedback_inf > 0.5)
                correct_post_inf(1,trial) = 1;
            elseif (trustworthy == 0) && (post_feedback_inf < 0.5)
                correct_post_inf(1,trial) = 1;
            else
                correct_post_inf(1,trial) = 0;
            end
        end
        inf_correct(habit_val,repeat) = sum(correct_post_inf)/num_trials;
        performance(habit_val,repeat) = sum(perform_vector)/num_trials;
    end
end
        
fin_perf = performance;
fin_inf_correct = inf_correct;
fin_habits = habits;

