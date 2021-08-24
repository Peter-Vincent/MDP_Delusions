% Script: vary all parameters and save results
clear
dbstop if error

UseMood   = false;   % use model with/without mood?
TreatBeta = false;  % if UseMood, use model with 'treatment' of raising beta (lowering precision)
TrustorSwitch = 1;  % which sequence, stable (1) or with changes in trustworthiness (2)?

likelihoods = [0.65 0.75 0.85 0.95];    likel_rand = -0.05:0.01:0.04; % add
habits      = [5  32  200];             habit_rand = 0.5:0.25:3;      % multiply
if UseMood
    moods       = [-3  0  3];           mood_rand  = -1.5:0.3:1.5;    % add
end
alphas      = [0.75  1.5  2.25];        alpha_rand = -0.25:0.075:0.425; % add
inv_betas   = [0.5  1  1.5];            beta_rand = -0.25:0.05:0.25;    % add
consistency = 0;                      % 1:3 (see below) for const_rand; 0 if trust_sequence

% Macbook paths
% save_path  = '/Users/rickadams/Data/MDP_sim_symm_a1_beta/'; % Macbook
save_path  = '/Data/MDP_simulations/MDP_sim_trust_beta_NoAf/'; % iMac
code_path  = '/Users/rickadams/Dropbox/Rick/Academic/Delusions model/working_model/';
spm_path   = '/Users/rickadams/Code/SPM/spm12_v7771/';

addpath(save_path, code_path, spm_path)
addpath([spm_path 'toolbox/DEM/'])

cd(code_path)

load('switch_sequence_250.mat') % 250 trials, 125 untrustworthy at start
load('trust_sequence_250.mat')  % 250 trials, mostly trustworthy (stable)

% load('switch_sequence.mat') - NB 256 trials, 156 trustworthy at start
% make advisor untrustworthy first? i.e. swap 1s for 2s & vice versa
% switch_sequence(1,:) = abs(switch_sequence(1,:)-3);

if TrustorSwitch == 1
    sequence = trust_sequence;
else
    sequence = switch_sequence;
end

if UseMood
    total_sims = length(likelihoods)*length(habits)*length(alphas)*length(inv_betas)*length(consistency)*length(moods);
else
    total_sims = length(likelihoods)*length(habits)*length(alphas)*length(inv_betas)*length(consistency);
end
t = 1;      % Initialise count

% set rand_seed
keyboard % IF RESETTING SEED, STOP SCRIPT AND RESTART
rand_seed=864;  % 4*3^4 = 0, 324, 648, 972, 1296, 1620, 1944
rng(rand_seed,'twister');

for l = 1:length(likelihoods)
    for h = 1:length(habits)
%         for m = 1:length(moods)      % COMMENT LOOP OUT IF NOT USING MOOD
            for a = 1:length(alphas)
                for b = 1:length(inv_betas)
                    for c = consistency
                        tic
                        % rearrange trials to vary consistent no. at start
                        switch c
                            case 1
                                const_trials = randi([0 9],1);
                            case 2
                                const_trials = randi([1 11],1)*2+10;
                            case 3
                                const_trials = randi([0 9],1)*10+35;
                            case 0
                                const_trials = 250;
                        end
                        % const_trials = consistent(c) * const_rand(randi(length(const_rand)));
                        % const_trials = round(const_trials); % ensure integer
                        
                        if const_trials < 125 % max number of consistent trials
                            states = sequence(:,[1:const_trials ...
                                126:length(sequence) const_trials+1:125]);
                        else
                            states = sequence;
                        end
                        
                        n         = length(states);
                        rand_seed = rand_seed+1;
                        % set parameters with some random variation
                        a_conting = likelihoods(l)  + likel_rand(randi(length(likel_rand)));
                        e_dir     = floor(habits(h) * habit_rand(randi(length(habit_rand))));
                        if UseMood
                            c_3   = moods(m)        + mood_rand(randi(length(mood_rand)));
                        else 
                            c_3   = NaN;
                        end
                        alpha     = alphas(a)       + alpha_rand(randi(length(alpha_rand)));
                        inv_beta  = inv_betas(b)    + beta_rand(randi(length(beta_rand)));
                        % NB inv_beta (precision) is transformed to beta in Run model command
                        
                        % Run model
                        if UseMood
                            if TreatBeta
                                % Set beta to decrease by some fraction, 
                                % after a certain number of false inferences?  
                                invbeta_drop     = 0.5; 
                                invbeta_floor    = min(inv_betas)+min(beta_rand);
                                OPTIONS.new_beta = 1/((inv_beta-invbeta_floor)*invbeta_drop+invbeta_floor);
                                OPTIONS.false_inf_count = 0;  % initialise
                                OPTIONS.false_inf_max   = 1; % # false infs before treatment starts
                                % Set new beta to start at a given trial?
                                % OPTIONS.new_beta_trial  = 135; 
                                mdp = MDP_Delusions_Affect_VaryAll_TreatBeta(n,states,rand_seed,a_conting,e_dir,c_3,alpha,1/inv_beta,OPTIONS);
                            else
                                mdp = MDP_Delusions_Affect_VaryAll(n,states,rand_seed,a_conting,e_dir,c_3,alpha,1/inv_beta);
                            end
                        else
                            mdp = MDP_Delusions_VaryAll(n,states,rand_seed,a_conting,e_dir,alpha,1/inv_beta);
                        end
                        
                        savename = ['MDP_l' num2str(a_conting)   '_e' num2str(e_dir) ...
                            '_c' num2str(c_3) '_a' num2str(alpha) '_b' num2str(inv_beta) ...
                            '_seq' num2str(const_trials) '_r' num2str(rand_seed) '.mat'];
                        
                        cd(save_path)
                        save(savename,'mdp');
                        disp(['Completed ' num2str(t) '/' num2str(total_sims) ' simulations']);
                        t = t+1; toc
                        clear mdp states a_conting e_dir c_3 alpha const_trials
                    end
                end
            end
%         end
    end
end

