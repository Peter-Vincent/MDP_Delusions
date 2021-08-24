% plot state space trajectories - Rick Adams

function [all_traj] = plot_state_space_traj(Trajectories, Results, raw, width, moving, MDP_sim_symm_beta)

% Takes Trajectories structure with fields (from analyse_varied_MDPs):
% Trajectories.e_traj          = e dirichlet parameters over trials          
% Trajectories.post_trust_traj = posterior inference re advisor over trials 
% Trajectories.pol_prec_traj   = posterior policy precision (1/beta) "  "

% Results = vector of outcome by which colours are determined, e.g.
%           number of false inferences

% raw    = Using raw data that needs processing (1) or mean summary data (0)
% width  = LineWidth for plot
% moving = 0 or 1 for animated plotting
% MDP_sim_symm_beta = contains axes limits to make all plots same

% Produces plot and also record of all plotted trajectories in all_traj


% Remap the colours according to outcome (e.g. number of false inferences)
cmap        = jet(256);
v           = rescale(Results, 1, 256);
numSubj     = length(Results);
colours     = zeros(numSubj, 3);

% Assign colours according to the value of the outcome
for s = 1 : numSubj
    row = round(v(s));
    colours(s, :) = cmap(row, :);
end

numTrials = length(Trajectories(1).e_traj);
all_traj  = [];

% Plot
for s = 1: numSubj
    
    if raw  % if processing required
        
        e_traj_trans = NaN(numTrials,1);
        
        for t = 1:numTrials
            % get relative strengths of habits of trusting vs not trusting
            e_traj_trans(t,1) = log(sum(Trajectories(s).e_traj(t,1:2))/...
                sum(Trajectories(s).e_traj(t,3:4)));
        end
        e_traj_trans    = abs(smooth(e_traj_trans));
        e_traj_trans(1) = []; % make same length for plotting
        
        % get cumulative measure of trial to trial absolute shifts in posteriors
        cum_post_shift = smooth(cumsum(abs(diff(Trajectories(s).post_trust_traj))));
        
        % smooth posterior policy precision trajectory
        smooth_pol_prec    = smooth(Trajectories(s).pol_prec_traj,50);
        smooth_pol_prec(1) = []; % make same length for plotting
        
    else        
        e_traj_trans    = Trajectories(s).e_traj;
        cum_post_shift  = Trajectories(s).post_trust_traj;
        smooth_pol_prec = Trajectories(s).pol_prec_traj;
    end
        
    % plot (omitting 1st trial in non-cumulative ones so lengths all match)
    plot3(e_traj_trans,cum_post_shift,smooth_pol_prec,'color',colours(s,:),'LineWidth',width)
    hold on
    
    if moving 
        M(s) = getframe;
    end
    
    if raw
        % record plotted trajectories for averaging
        t = numTrials-1;
        
        all_traj(s,1)             = Results(s);
        all_traj(s,2:t+1)         = e_traj_trans;
        all_traj(s,t+2:t*2+1)     = cum_post_shift;
        all_traj(s,t*2+2:t*3+1)   = smooth_pol_prec;
    end
    
end

xlabel('Habit strength')
ylabel('Cumulative changes in trust posteriors')
zlabel('Posterior policy precision')
grid on
try % to make all plots with same axes limits
    if raw
        xlim([0 MDP_sim_symm_beta.xlimits(2)]);
        ylim([0 MDP_sim_symm_beta.ylimits(2)]);
        zlim([0 MDP_sim_symm_beta.zlimits(2)]);
    else % mean limits
        xlim([0 MDP_sim_symm_beta.mean_xlimits(2)]);
        ylim([0 MDP_sim_symm_beta.mean_ylimits(2)]);
        zlim([0 MDP_sim_symm_beta.mean_zlimits(2)]);
    end
catch
    axis tight
end

% figure 
% movie(M)