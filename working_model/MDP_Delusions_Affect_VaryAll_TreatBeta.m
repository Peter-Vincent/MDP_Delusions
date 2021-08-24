function mdp = MDP_Delusions_Affect_VaryAll_TreatBeta(num_trials,initial_states,rand_seed,a_conting,...
                             e_dir,c_3,alpha,beta,OPTIONS)
% num_trials        = total no. trials 
% initial_states    = output from define_states
%                    (?use whole sequences from Adams_seq: 120/240 trials)

% SET PARAMETERS
A_conting = 0.9;       % for process model
a_dir     = 600;       % Dirichlet parameter for A matrix
b_1_1     = [1 1; 
             0 0];  % transitions for trustworthy policy
b_1_2     = [0 0; 
             1 1]; % transitions for untrustworthy policy
b_2       = [1 0;
             0 1];     % transitions for deck (inference +/- process)
b_dir     = 600;       % Dirichlet parameter for B matrices
% c_3       = 1;       % Predicted affective outcomes 
% beta     = 4;         % Precision over policies (alpha = over actions)

try rand_seed; catch rand_seed = 1; end   
% Initialisation
%--------------------------------------------------------------------------
% spm_path   = 'D:\Other_Science\Halucinations_SPM_paper\spm12';
% model_path = '';
% spm_path   = '/Users/rickadams/Code/SPM/spm12_v7771/';
% model_path = '/Users/rickadams/Dropbox/Rick/Academic/Delusions model/working_model';
% addpath(genpath(spm_path))
% addpath(genpath(model_path))
rng default
dbstop if error
tic
if num_trials < size(initial_states,2)
    num_trials = size(initial_states,2);
end
% Specify generative model
%__________________________________________________________________________

% Prior over initial states
%--------------------------------------------------------------------------

D{1} = [1;1];   % Advisor trustworthy or not
D{2} = [1;1];   % Correct option (A or B)
D{3} = [0;0;1]; % Choose A or B or null
D{4} = [1;1;0]; % Aroused or relaxed or null
D{5} = [1;0;0]; % Null, advice or choice/feedback

% Likelihood
%--------------------------------------------------------------------------
%% Construct A matrices
for g = 1:4
    A{g} = zeros(1,length(D{1}),length(D{2}),length(D{3}),length(D{4}),length(D{5}));
end

for f1 = 1:length(D{1})
    for f2 = 1:length(D{2})
        for f3 = 1:length(D{3})
            for f4 = 1:length(D{4})
                for f5 = 1:length(D{5})
                    if f5 == 1 % Null stage
                        A{1}(3,f1,f2,f3,f4,f5) = 1; % No advice given
                        A{2}(3,f1,f2,f3,f4,f5) = 1; % during advice giving stage, always give no feedback (3)
                    elseif f5 == 2 % Advice stage
                        if f1 == 1 % If advisor is trustworthy...
                            %... advice is consistent with correct option...
                            % MDP.o therefore reports the advice received
                            A{1}(f2,f1,f2,f3,f4,f5) = 1;
                        else % ...if advisor is not trustworthy...
                            %... advice is inconsistent with correct option.
                            A{1}(-(f2-3),f1,f2,f3,f4,f5) = 1;
                        end
                        A{2}(3,f1,f2,f3,f4,f5) = 1; % Actions are greeted with "incorrect"
                    elseif f5 == 3 % Feedback stage
                        A{1}(3,f1,f2,f3,f4,f5) = 1; % No advice given
                        if f2 == f3
                            A{2}(1,f1,f2,f3,f4,f5) = A_conting; % Correct
                            A{2}(2,f1,f2,f3,f4,f5) = 1-A_conting;
                        else
                            A{2}(1,f1,f2,f3,f4,f5) = 1-A_conting;
                            A{2}(2,f1,f2,f3,f4,f5) = A_conting; % Incorrect
                        end
                    end
                    % Interoception
                    %------------------------------------------------------
                    % Autonomic state gives rise to cardiac outcome
                    % (tachy/bradycardia)
                    A{3}(f4,f1,f2,f3,f4,f5) = 1;
                    
                    % Proprioception
                    %------------------------------------------------------
                    % Can observe the choice being made
                    A{4}(f3,f1,f2,f3,f4,f5) = 1;
                end
            end
        end
    end
end
%% Construct A matrices
for g = 1:4
    a{g} = zeros(1,length(D{1}),length(D{2}),length(D{3}),length(D{4}),length(D{5}));
end

for f1 = 1:length(D{1})
    for f2 = 1:length(D{2})
        for f3 = 1:length(D{3})
            for f4 = 1:length(D{4})
                for f5 = 1:length(D{5})
                    if f5 == 1 % Null stage
                        a{1}(3,f1,f2,f3,f4,f5) = 1; % No advice given
                        a{2}(3,f1,f2,f3,f4,f5) = 1; % during advice giving stage, always give no feedback (3)
                    elseif f5 == 2 % Advice stage
                        if f1 == 1 % If advisor is trustworthy...
                            %... advice is consistent with correct option...
                            % MDP.o therefore reports the advice received
                            a{1}(f2,f1,f2,f3,f4,f5) = 1;
                        else % ...if advisor is not trustworthy...
                            %... advice is inconsistent with correct option.
                            a{1}(-(f2-3),f1,f2,f3,f4,f5) = 1;
                        end
                        a{2}(3,f1,f2,f3,f4,f5) = 1; % Actions are greeted with "incorrect"
                    elseif f5 == 3 % Feedback stage
                        a{1}(3,f1,f2,f3,f4,f5) = 1; % No advice given
                        if f2 == f3
                            a{2}(1,f1,f2,f3,f4,f5) = a_conting; % Correct
                            a{2}(2,f1,f2,f3,f4,f5) = 1-a_conting;
                        else
                            a{2}(1,f1,f2,f3,f4,f5) = 1-a_conting;
                            a{2}(2,f1,f2,f3,f4,f5) = a_conting; % Incorrect
                        end
                    end
                    % Interoception
                    %------------------------------------------------------
                    % Autonomic state gives rise to cardiac outcome
                    % (tachy/bradycardia)
                    a{3}(f4,f1,f2,f3,f4,f5) = 1;
                    
                    % Proprioception
                    %------------------------------------------------------
                    % Can observe the choice being made
                    a{4}(f3,f1,f2,f3,f4,f5) = 1;
                end
            end
        end
    end
end

% Preclude learning in A
a{1} = a{1}*a_dir;
a{2} = a{2}*a_dir;
a{3} = a{3}*a_dir;
a{4} = a{4}*a_dir;
%%
% Transition probabilities
%--------------------------------------------------------------------------

% Generative process:
B{1} = repmat([0.9 0.1; 0.1 0.9],1,1,2); % Advisor trustworthiness invariant to action.
B{2} = eye(2); % or b_2;             
B{3} = zeros(3,3,3);                     % Control over choice states
for k = 1:3
    B{3}(k,:,k) = 1;
end


% B{4} = zeros(3,3,3);
% for k = 1:3
%     B{4}(k,k,1) = 1;
%     B{4}(k,k,2) = 1;
%     B{4}(k,k,3) = 1;
% end
B{4} = zeros(3,3,3);                    % Autonomic state depends (stochastically) on choices
for k = 1:2
    B{4}(k,:,k) = 2/3;
end
for k = 1:2
    B{4}(abs(k-3),:,k) = 1/3;
end
B{4}(3,:,3) = 1;
% B{4} = zeros(2,2,2);
% for k = 1:2
%     B{4}(k,:,k) = 1;
% end
B{5} = [0 0 1;1 0 0;0 1 0];                       % Progress through null, advice, choice

% Generative model:
b{1}(:,:,1) = b_1_1.^(1)*b_dir;  % Under action 1, transition to more trustworthy
b{1}(:,:,2) = b_1_2.^(1)*b_dir;  % Under action 2, transition to less trustworthy
b{2} = b_2*b_dir;
b{3} = B{3}*b_dir;
b{4} = B{4}*b_dir;
b{5} = B{5}*b_dir;

% (128 scales the Dirichlet parameters to preclude learning)

%% Preferences
%--------------------------------------------------------------------------
for c = 1:numel(A)
    C{c} = zeros(size(A{c},1),1);
end
C{2} = [3 -3 0]';
C{3} = [c_3 -c_3 0]';
% Policies
%--------------------------------------------------------------------------
V = zeros(2,4,5);
V(:,:,1) = [ones(2,2) 2*ones(2,2)];
V(:,:,2) = ones(2,4);
V(:,:,3) = vertcat(3*ones(1,4),[1 2 1 2]);
%V(:,:,4) = [ones(2,2) 2*ones(2,2)];
V(:,:,4) = [ones(2,2) 2*ones(2,2)];
V(:,:,5) = ones(2,4);
E = ones(4,1) .* e_dir;

% Compile MDP
%--------------------------------------------------------------------------
MDP.A = A;
MDP.a = a;
MDP.B = B;
MDP.b = b;
MDP.C = C;
MDP.D = D;
MDP.e = E;
MDP.V = V;
MDP.zeta = inf;
MDP.alpha = alpha;
MDP.beta  = beta;

% Simulate and plot
%__________________________________________________________________________
MDP(1:num_trials) = deal(MDP); % Simulate 16 worth of trials
%% Now add in the prededined states
for decided = 1:size(initial_states,2)
    cur_state = initial_states(:,decided);
    MDP(decided).s = cur_state;
end
%%
mdp = spm_MDP_VB_X_rand(MDP,rand_seed,OPTIONS);

end