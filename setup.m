%% Initialise paths
paths = {};
paths{1} = 'C:\Users\Murri\Documents\MDP_AI\SPM12';
paths{2} = 'C:\Users\Murri\Documents\MDP_AI\MDP_Delusions\plotting_functions';
paths{3} = 'C:\Users\Murri\Documents\MDP_AI\MDP_Delusions';

for cur_path = 1:length(paths)
    addpath(genpath(paths{cur_path}))
end
