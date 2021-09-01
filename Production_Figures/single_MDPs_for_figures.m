% MDP single simulations for figures

clear
close all
dbstop if error

% save_path  = '/Users/rickadams/Dropbox/Rick/Academic/Delusions model/MDP_sims_for_Figures/';  % Macbook
save_path  = '/Users/rickadams/Dropbox/Rick/Academic/Delusions model/MDP_sims_for_Figures/';    % iMac
code_path  = '/Users/rickadams/Dropbox/Rick/Academic/Delusions model/working_model/';
plot_path  = [code_path 'plotting_functions/'];
spm_path   = '/Users/rickadams/Code/SPM/spm12_v7771/';
addpath(save_path, code_path, spm_path, plot_path)
addpath([spm_path 'toolbox/DEM/'])
addpath(genpath(fun_path))

cd(save_path)

load('switch_sequence_250.mat') % 250 trials, 125 untrustworthy at start
load('trust_sequence_250.mat')  % 250 trials, mostly trustworthy (stable)

% Figure 2 - Trust sequences, no mood, without/with habits

states    = trust_sequence; % sequence of states
n         = length(states); % number of trials
rand_seed = 1;              % set random seed
a_conting = 0.9;            % likelihood
e_dir     = 600;            % Dirichlet over habits (high = no habits)
alpha     = 1.5;              % Action precision (softmax inv temp)
inv_beta  = 1;              % Policy precision
c_3       = NaN;            % Mood (if used)
with_mood = 0;              % Plot figure with/without mood?

try
    load('mdp2a.mat')
catch
    mdp2a     = MDP_Delusions_VaryAll(n,states,rand_seed,a_conting,e_dir,alpha,1/inv_beta);
    save('mdp2a','mdp2a')
end
plot_trial(mdp2a,with_mood)
suptitle('Figure 2A')

e_dir     = 2;
try
    load('mdp2b.mat')
catch
    mdp2b     = MDP_Delusions_VaryAll(n,states,rand_seed,a_conting,e_dir,alpha,1/inv_beta);
    save('mdp2b','mdp2b')
end
plot_trial(mdp2b,with_mood)
suptitle('Figure 2B')

% Figure 3 - Switch sequences, no mood, without/with habits, changing likelihood

states    = switch_sequence;
n         = length(states);
rand_seed = 1;
a_conting = 0.9;
e_dir     = 2;
alpha     = 1.5;
inv_beta  = 1;
c_3       = NaN;
with_mood = 0;      

try
    load('mdp3a.mat')
catch
    mdp3a     = MDP_Delusions_VaryAll(n,states,rand_seed,a_conting,e_dir,alpha,1/inv_beta);
    save('mdp3a','mdp3a')
end
plot_trial(mdp3a,with_mood)
suptitle('Figure 3A')

a_conting = 0.75;
try
    load('mdp3b.mat')
catch
    mdp3b     = MDP_Delusions_VaryAll(n,states,rand_seed,a_conting,e_dir,alpha,1/inv_beta);
    save('mdp3b','mdp3b')
end
plot_trial(mdp3b,with_mood)
suptitle('Figure 3B')

a_conting = 0.6;
e_dir     = 600;
try
    load('mdp3c.mat')
catch
    mdp3c     = MDP_Delusions_VaryAll(n,states,rand_seed,a_conting,e_dir,alpha,1/inv_beta);
    save('mdp3c','mdp3c')
end
plot_trial(mdp3c,with_mood)
suptitle('Figure 3C')

% Figure 4 - Switch sequence, with mood, varying habits & likelihood

states    = switch_sequence;
n         = length(states);
rand_seed = 1;
a_conting = 0.9;
e_dir     = 600;
alpha     = 1.5;
inv_beta  = 1;
c_3       = -1;
with_mood = 1; 

try
    load('mdp4a.mat')
catch
    mdp4a     = MDP_Delusions_Affect_VaryAll(n,states,rand_seed,a_conting,e_dir,c_3,alpha,1/inv_beta);
    save('mdp4a','mdp4a')
end
plot_trial(mdp4a,with_mood)
suptitle('Figure 4A')

e_dir     = 2;
try
    load('mdp4b.mat')
catch
    mdp4b     = MDP_Delusions_Affect_VaryAll(n,states,rand_seed,a_conting,e_dir,c_3,alpha,1/inv_beta);
    save('mdp4b','mdp4b')
end
plot_trial(mdp4b,with_mood)
suptitle('Figure 4B')

a_conting = 0.75;
try
    load('mdp4c.mat')
catch
    mdp4c     = MDP_Delusions_Affect_VaryAll(n,states,rand_seed,a_conting,e_dir,c_3,alpha,1/inv_beta);
    save('mdp4c','mdp4c')
end
plot_trial(mdp4c,with_mood)
suptitle('Figure 4C')

% To plot 4 trials from start of Fig 4A
states    = switch_sequence;
n         = length(states);
rand_seed = 1;
a_conting = 0.9;
e_dir     = 600;
alpha     = 1.5;
inv_beta  = 1;
c_3       = -1;
with_mood = 1; 
OPTIONS.plot = 1; 

mdp = MDP_Delusions_Affect_VaryAll_TreatBeta(4,states(:,1:4),rand_seed,a_conting,e_dir,c_3,alpha,1/inv_beta,OPTIONS);
