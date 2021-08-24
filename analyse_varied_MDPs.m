% Analyse results of varying parameters
clearvars -except MDP_sim_trust_beta_NoAf MDP_sim_trust_beta MDP_sim_symm_beta_NoAf MDP_sim_symm_beta MDP_sim_symm_beta_treat05_10
% close all
dbstop if error

doResults = 0; % Make Results and Trajectories? If not, load them
doAnalysis= 1;
doSave    = 0; % If running for first time

% save_root  = '/Users/rickadams/Dropbox/Rick/Academic/Delusions model/MDP_sims_for_Figures/'; % PETER - CHANGE PATH TO YOUR DROPBOX
% save_root  = '/Users/rickadams/Data/MDP_simulations/';  % Macbook
save_root  = '/Data/MDP_simulations/';                  % iMac
save_fold  = 'MDP_sim_symm_beta/';
code_path  = '/Users/rickadams/Dropbox/Rick/Academic/Delusions model/working_model/';
plot_path  = [code_path 'plotting_functions/'];
spm_path   = '/Users/rickadams/Code/SPM/spm12_v7771/';
fun_path   = '/Users/rickadams/Dropbox/Downloaded_functions/';
save_path  = [save_root save_fold];
addpath(save_path, code_path, spm_path, plot_path)
addpath([spm_path 'toolbox/DEM/'])
addpath(genpath(fun_path))

cd(save_path)

if doResults
    
    MDP_files = dir;
    % Remove '.' and '..' files and folders
    for f = length(MDP_files):-1:1
        if startsWith(MDP_files(f).name,'.')
            MDP_files(f) = [];
        elseif MDP_files(f).isdir == 1
            MDP_files(f) = [];
        end
    end
    
    Results = NaN(length(MDP_files),15);
    
    for f = 1:length(MDP_files)
        
        disp([num2str(f) ' of ' num2str(length(MDP_files))])
        load(MDP_files(f).name)
        
        % find & store parameters from filename
        l_ind = strfind(MDP_files(f).name,'_l');
        e_ind = strfind(MDP_files(f).name,'_e');
        c_ind = strfind(MDP_files(f).name,'_c');
        a_ind = strfind(MDP_files(f).name,'_a');
        b_ind = strfind(MDP_files(f).name,'_b');
        s_ind = strfind(MDP_files(f).name,'_s');
        m_ind = strfind(MDP_files(f).name,'.m');
        r_ind = strfind(MDP_files(f).name,'_r');
        
        Results(f,1) = str2num(MDP_files(f).name(l_ind+2:e_ind-1));
        Results(f,2) = log(str2num(MDP_files(f).name(e_ind+2:c_ind-1)));
        if isempty(strfind(MDP_files(f).name,'_b'))
            Results(f,4) = str2num(MDP_files(f).name(c_ind+2:a_ind-1));
            Results(f,5) = str2num(MDP_files(f).name(a_ind+2:s_ind-1));
        elseif isempty(strfind(MDP_files(f).name,'_a'))
            Results(f,4) = str2num(MDP_files(f).name(c_ind+2:b_ind-1));
            Results(f,6) = str2num(MDP_files(f).name(b_ind+2:s_ind-1));
        else % if both alpha and beta estimated
            Results(f,4) = str2num(MDP_files(f).name(c_ind+2:a_ind-1));
            Results(f,5) = str2num(MDP_files(f).name(a_ind+2:b_ind-1));
            Results(f,6) = str2num(MDP_files(f).name(b_ind+2:s_ind-1));
        end
        if isempty(r_ind) % first versions done without randseed in name
            Results(f,3) = sqrt(str2num(MDP_files(f).name(s_ind+4:m_ind-1)));
        else
            Results(f,3) = sqrt(str2num(MDP_files(f).name(s_ind+4:r_ind-1)));
        end
        
        % compute & store performance metrics
        num_trials       = length(mdp);
        Scnd_half        = num_trials/2+1:num_trials;
        perform_vector   = zeros(1,num_trials);
        correct_post_inf = zeros(1,num_trials);
        trustOrNot_inf   = NaN(1,num_trials);
        false_Inf_conf   = NaN(1,num_trials);
        
        for trial = 1:num_trials
            
            outcomes = mdp(trial).o;
            choices_vector(1,trial) = outcomes(4,3);
            perform_vector(1,trial) = outcomes(2,end);
            if perform_vector(1,trial) == 2
                perform_vector(1,trial) = 0;
            end
            inf_over_ad       = mdp(trial).xn{1,1}; % same as mdp(trial).un(1)
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
            else % if inference re advisor is incorrect
                correct_post_inf(1,trial) = 0;
                trustOrNot_inf(1,trial)   = post_feedback_inf;
                if abs(round(post_feedback_inf)-post_feedback_inf) < 0.1 % if false inf >80% confident (NB 0.5 is 0% so 0.1 is 80%)
                    false_Inf_conf(1,trial) = 1;
                else
                    false_Inf_conf(1,trial) = 0;
                end
            end
            
            % State/habit trajectories
            e_traj(trial,:)          = mdp(trial).e;
            post_trust_traj(trial,1) = post_feedback_inf;
            pol_prec_traj(trial,1)   = mdp(trial).w(3);
            
        end
        false_infr = 1-sum(correct_post_inf)/num_trials;
        falInf_conf= nansum(false_Inf_conf)/(num_trials-sum(correct_post_inf));
        fI_mean_con= post_trust_traj(~correct_post_inf);
        fI_mean_con= mean(abs(0.5-fI_mean_con)*2);
        falInf_incn= sum(diff(correct_post_inf)==1)/(num_trials-sum(correct_post_inf)); % diff=1 means 0 (false) followed by 1 (correct)
        trustOrNot = nanmean(trustOrNot_inf-0.5); % centre at 0
        correct_ch = sum(perform_vector)/num_trials;
        choice_hab = abs(mean(choices_vector-1.5))*2;
        
        final_habit = mdp(end).e;
        norm_habit  = final_habit ./ sum(final_habit);
        hab_entropy = -sum(norm_habit.*log2(norm_habit));
        blue_habit  = sum(norm_habit([1 3]));
        green_habit = sum(norm_habit([2 4]));
        trust_habit = sum(norm_habit([1 2]));
        distr_habit = sum(norm_habit([3 4]));
        crd_entropy = -sum([blue_habit green_habit].*log2([blue_habit green_habit]));
        tru_entropy = -sum([trust_habit distr_habit].*log2([trust_habit distr_habit]));
        
        % Record results
        Results(f,7)  = false_infr;      % fraction false inferences out of total
        Results(f,8)  = trustOrNot;      % strength of false inferences, trust/not trust on average? 
        Results(f,9)  = correct_ch;      % fraction correct choices
        Results(f,10) = choice_hab;      % 50/50 vs habitual card choice (0 to 1 scale)
        Results(f,11) = hab_entropy;     % entropy of habits
        Results(f,12) = falInf_conf;     % proportion of false inferences with >80% confidence 
        Results(f,13) = 1-falInf_incn;   % incorrigibility/false inference consistency (what proportion are followed by false inferences?)
        Results(f,14:17) = norm_habit;   % B/trust, G/trust, B/distrust, G/distrust
        Results(f,18) = crd_entropy;     % Entropy over card habits
        Results(f,19) = tru_entropy;     % Entropy over trust habits
        Results(f,20) = fI_mean_con;     % Mean confidence (0-100%) in false inferences (independent of trust/not trust)
        Results(f,21) = nansum(Results(f,[7 20 13])); % Delusion 'score'
        
        % Record state space trajectories for plotting
        Trajectories(f).e_traj          = e_traj;          clear e_traj
        Trajectories(f).post_trust_traj = post_trust_traj; clear post_trust_traj
        Trajectories(f).pol_prec_traj   = pol_prec_traj;   clear pol_prec_traj
        Trajectories(f).results         = Results(f,:);
        clear mdp
    end
    
else
    cd('Results')
    load('Results.mat')
    load('Trajectories.mat')
end

if doAnalysis
    
    % Transform beta to precision (1/beta, or inv_beta) if stored as beta
    if min(Results(:,6)) > 0.4 % lowest inv_beta is 0.25, lowest beta is 0.5
        Results(:,6) = 1./Results(:,6);
    end
    
    % Remove subjects with extreme parameter values?
    % Results_orig  = Results;            % save original version
    % ind_to_remove = Results(:,6)<0.75;
    % Results(ind_to_remove,:) = [];
    
    % Sort Results by variables: best at top, worst at bottom
    Results_Del = sortrows(Results,[7 12],{'ascend' 'ascend'});
    Results_Cor = sortrows(Results,[9 7],{'descend' 'ascend'});
    Results_ChH = sortrows(Results,[10 8],{'ascend' 'ascend'});
    Results_Ent = sortrows(Results,[11 9],{'descend' 'descend'});
    
    % Correlations between outcomes and parameters
    FIn_lik = corr(Results(:,7),Results(:,1),'type','Spearman');
    FIn_hab = corr(Results(:,7),Results(:,2),'type','Spearman');
    FIn_con = corr(Results(:,7),Results(:,3),'type','Spearman');
    FIn_moo = corr(Results(:,7),Results(:,4),'type','Spearman');
    FIn_alp = corr(Results(:,7),Results(:,5),'type','Spearman');
    FIn_bet = corr(Results(:,7),Results(:,6),'type','Spearman');
    
    Del_lik = corr(Results(:,21),Results(:,1),'type','Spearman');
    Del_hab = corr(Results(:,21),Results(:,2),'type','Spearman');
    Del_con = corr(Results(:,21),Results(:,3),'type','Spearman');
    Del_moo = corr(Results(:,21),Results(:,4),'type','Spearman');
%     Del_moo = corr(Results(:,21),abs(Results(:,4)),'type','Spearman');
    Del_alp = corr(Results(:,21),Results(:,5),'type','Spearman');
    Del_bet = corr(Results(:,21),Results(:,6),'type','Spearman');
    
    Cor_lik = corr(Results(:,9),Results(:,1),'type','Spearman');
    Cor_hab = corr(Results(:,9),Results(:,2),'type','Spearman');
    Cor_con = corr(Results(:,9),Results(:,3),'type','Spearman');
    Cor_moo = corr(Results(:,9),Results(:,4),'type','Spearman');
    Cor_alp = corr(Results(:,9),Results(:,5),'type','Spearman');
    Cor_bet = corr(Results(:,9),Results(:,6),'type','Spearman');
    
    Crd_lik = corr(Results(:,10),Results(:,1),'type','Spearman');
    Crd_hab = corr(Results(:,10),Results(:,2),'type','Spearman');
    Crd_con = corr(Results(:,10),Results(:,3),'type','Spearman');
    Crd_moo = corr(Results(:,10),Results(:,4),'type','Spearman');
    Crd_alp = corr(Results(:,10),Results(:,5),'type','Spearman');
    Crd_bet = corr(Results(:,10),Results(:,6),'type','Spearman');
    
    Ent_lik = corr(Results(:,11),Results(:,1),'type','Spearman');
    Ent_hab = corr(Results(:,11),Results(:,2),'type','Spearman');
    Ent_con = corr(Results(:,11),Results(:,3),'type','Spearman');
    Ent_moo = corr(Results(:,11),Results(:,4),'type','Spearman');
    Ent_alp = corr(Results(:,11),Results(:,5),'type','Spearman');
    Ent_bet = corr(Results(:,11),Results(:,6),'type','Spearman');
    
    disp('Parameter correlations with false inferences re trustworthiness:')
    disp(['Likelihood rho = ' num2str(FIn_lik)])
    disp(['Habit rho = ' num2str(FIn_hab)])
    disp(['Initial consistency rho = ' num2str(FIn_con)])
    disp(['Mood (not transformed) rho = ' num2str(FIn_moo)])
    disp(['Alpha rho = ' num2str(FIn_alp)])
    disp(['Beta rho = ' num2str(FIn_bet)])
    
    disp('Parameter correlations with delusion-like score:')
    disp(['Likelihood rho = ' num2str(Del_lik)])
    disp(['Habit rho = ' num2str(Del_hab)])
    disp(['Initial consistency rho = ' num2str(Del_con)])
    disp(['Mood (not transformed) rho = ' num2str(Del_moo)])
    disp(['Alpha rho = ' num2str(Del_alp)])
    disp(['Beta rho = ' num2str(Del_bet)])
    
    disp('Parameter correlations with correct initial choices:')
    disp(['Likelihood rho = ' num2str(Cor_lik)])
    disp(['Habit rho = ' num2str(Cor_hab)])
    disp(['Initial consistency rho = ' num2str(Cor_con)])
    disp(['Mood (not transformed) rho = ' num2str(Cor_moo)])
    disp(['Alpha rho = ' num2str(Cor_alp)])
    disp(['Beta rho = ' num2str(Cor_bet)])
    
    disp('Parameter correlations with habits over cards:')
    disp(['Likelihood rho = ' num2str(Crd_lik)])
    disp(['Habit rho = ' num2str(Crd_hab)])
    disp(['Initial consistency rho = ' num2str(Crd_con)])
    disp(['Mood (not transformed) rho = ' num2str(Crd_moo)])
    disp(['Alpha rho = ' num2str(Crd_alp)])
    disp(['Beta rho = ' num2str(Crd_bet)])
    
    % Transform mood to find correlation (if mood used) w false inf
    if ~isnan(Results(1,4))
        % Compute correlations without false inference-free subjects:
        % Replace 0 with NaN in fraction false inferences?
        % Results(Results(:,7)==0,7) = NaN;
        no_falseInf = Results(:,7);
        no_falseInf(no_falseInf(:,1)==0,1) = NaN;
        
        figure
        mood_tr = abs(Results(:,4));
        rho1 = corr(mood_tr,no_falseInf,'type','Spearman','rows','complete');
        subplot(1,3,1)
        scatter(mood_tr,no_falseInf,8,'k','filled'); lsline
        xlabel('Absolute mood, |c_3|'); ylabel('False/Total inferences')
        text(2,0.65,['\rho = ' num2str(rho1)])
        % Mood in subset
        %     imprecise_lik_ind = Results(:,1)<0.65 & Results(:,4)<0;
        %     rho2 = corr(mood_tr(imprecise_lik_ind),Results(imprecise_lik_ind,7),'type','Spearman');
        %     subplot(1,3,2)
        %     scatter(mood_tr(imprecise_lik_ind),Results(imprecise_lik_ind,7),'kx'); lsline
        %     xlabel('Absolute mood'); ylabel('False inferences in subset w low lik precision')
        subplot(1,3,2)
        rho2 = corr(mood_tr,abs(Results(:,8)),'type','Spearman','rows','complete');
%         rho2 = corr(mood_tr,Results(:,20),'type','Spearman','rows','complete');
        scatter(mood_tr,abs(Results(:,8).*2),8,'k','filled'); lsline
%         scatter(mood_tr,Results(:,20),8,'k','filled'); lsline
        xlabel('Absolute mood, |c_3|'); ylabel('Strength of false inferences')
        text(2,0.65,['\rho = ' num2str(rho2)])
        subplot(1,3,3)
        rho3 = corr(Results(:,4),Results(:,8),'type','Spearman','rows','complete');
        scatter(Results(:,4),Results(:,8)+0.5,8,'k','filled');
        xlabel('Mood, c_3'); ylabel('Mean posteriors of false inferences')
        
        % Compute correlations with false inference-free subjects:
        rho4 = corr(mood_tr,Results(:,7),'type','Spearman','rows','complete');
        % Replace NaN with 0 in mean false inferences
        trustOrNot_wFalInf                        = Results(:,8);
        trustOrNot_wFalInf(isnan(Results(:,8)),1) = 0;
        rho5 = corr(mood_tr,abs(trustOrNot_wFalInf),'type','Spearman','rows','complete');
        
        disp(['Absolute mood vs False inf (incl 0s), rho = ' num2str(rho4)])
        %     disp(['Absolute mood vs False inf if lik prec<0.65, rho = ' num2str(rho2)])
        disp(['Absolute mood vs Strength of false inferences (incl 0s), rho = ' num2str(rho5)])
        disp(['Absolute mood vs False inf direction, rho = ' num2str(rho3)])
    end
    
    %% Analyse relationships with results - choose outcome type
    for outcome = [21] % 7 = false trust inf, 9 = correct, 10 = card habits, 21 = delusion scores
        
        transf_mood = 0; % transform mood to absolute values?
        disp(['Transform mood = ' num2str(transf_mood)])
        
        % Multiple regression w non-NaN parameters ± transformed mood
        param_ind = ~isnan(Results(1,1:6));
        X = [Results(:,param_ind)];
        if ~isnan(Results(1,4)) && transf_mood == 1
            X(:,4) = mood_tr; % if mood used, transform it
        end
        if Results(1,3) > 15 % if all consistency is identical, i.e. exp(15.8) (=250), remove it
            X(:,3) = [];
        end
        X = (X - mean(X))./std(X,0,1);
        X = [ones(length(Results),1) X];
        Y = Results(:,outcome); 
        Y = (Y - mean(Y))./std(Y,0,1);
        disp(['For outcome ' num2str(outcome)])
        [b,b_CI] = regress(Y,X)
        
        % Look at interactions and plot betas
        %  a    Dir(e)    |c|   \alpha  1/\beta 
%         p1 = 1; p2 = 2; p3 = 4; p4 = 5; p5 = 6;
% %         X = [X X(:,p2+1).*X(:,p4+1).*X(:,p5+1)];
%         X = [X X(:,p1+1).*X(:,p2+1) X(:,p1+1).*X(:,p3+1) X(:,p1+1).*X(:,p4+1) X(:,p1+1).*X(:,p5+1) ...
%             X(:,p2+1).*X(:,p5+1) X(:,p1+1).*X(:,p2+1).*X(:,p5+1)];
%         [b,b_CI] = regress(Y,X)
%         
%         x = categorical({'Likelihood, a','Dir(e)','Choice precision, \alpha','Policy precision, 1/\beta','a*Dir(e)','a*|c|','a*\alpha','a*1/\beta',...
%             'Dir(e)*1/\beta'}); % 'Mood, |c|','a*Dir(e)*1/\beta'
%         x = reordercats(x,{'Likelihood, a','Dir(e)','Choice precision, \alpha','Policy precision, 1/\beta','a*Dir(e)','a*|c|','a*\alpha','a*1/\beta',...
%             'Dir(e)*1/\beta'}); % 'Mood, |c|','a*Dir(e)*1/\beta'
%         figure; b([1 4 5 13]) = []; bar(x,b) % remove intercept and nonsig parameters/interactions
%         xlabel('Parameters and interactions')
%         ylabel('Effect size (beta weight)')
%         title('Effect sizes of parameters and their interactions on delusion scores') 
        
        clear X Y
        
        % Create the scatter plots - first remap the colours
        max_scale = 2.6; % CHOOSE MAX VALUE FOR ALL COLORBARS
        cmap = jet(256);
        v = rescale([Results(:,outcome); max_scale], 1, 256);
        v = v(1:end-1);  % remove max_scale 
        numValues = length(Results(:,outcome));
        markerColors = zeros(numValues, 3);
        % Now assign marker colors according to the value of the data.
        for k = 1 : numValues
            row = round(v(k));
            markerColors(k, :) = cmap(row, :);
        end
        % Create the scatter plots - first decide which parameters to plot
        xlabels = {'Likelihood, a','Habit resistance, Dir(e)','Initial consistency','Mood, c','Choice precision, \alpha','Policy precision, 1/\beta'};
        a = 1; b = 2; c = 4; d = 5; e = 6; % choose 5 out of 6 parameters
        
        % Use transformed mood with some offset for plotting purposes?
        Results_tr = Results;
        if ~isnan(Results(1,4)) && transf_mood == 1
            Results_tr(Results(:,4)<0,4) = abs(Results(Results(:,4)<0,4))-0.15;
        end
        
        figure
        for pl = 1:4
            p = [b c d e];
            subplot(4,4,pl)
            scatter(Results_tr(:,a), Results_tr(:,p(pl)), [], markerColors,'filled');
            colormap(jet(256)); %cb = colorbar;
            xlabel(xlabels{a})
            ylabel(xlabels{p(pl)})
            axis tight
        end
        for pl = 1:4
            p = [a c d e];
            subplot(4,4,pl+4)
            scatter(Results_tr(:,b), Results_tr(:,p(pl)), [], markerColors,'filled');
            colormap(jet(256)); %cb = colorbar;
            xlabel(xlabels{b})
            ylabel(xlabels{p(pl)})
            axis tight
        end
        for pl = 1:4
            p = [a b d e];
            subplot(4,4,pl+8)
            scatter(Results_tr(:,c), Results_tr(:,p(pl)), [], markerColors,'filled');
            colormap(jet(256)); %cb = colorbar;
            xlabel(xlabels{c})
            ylabel(xlabels{p(pl)})
            axis tight
        end
        for pl = 1:4
            p = [a b c e];
            subplot(4,4,pl+12)
            scatter(Results_tr(:,d), Results_tr(:,p(pl)), [], markerColors,'filled');
            colormap(jet(256)); %cb = colorbar;
            xlabel(xlabels{d})
            ylabel(xlabels{p(pl)})
            axis tight
        end
        suptitle(['Parameter relationships for outcome ' num2str(outcome)])
        % c.Limits = [0 max(Results(:,6))];
        % grid on
        
        % Figures for paper
        if outcome == 9
%             figure % FIG 2
%             for pl = 1:3
%                 subplot(1,3,pl)
%                 p = [a e d];
%                 scatter(Results_tr(:,2), Results_tr(:,p(pl)), [], markerColors,'filled');
%                 colormap(jet(256)); cb = colorbar;
%                 xlabel(xlabels{2})
%                 ylabel(xlabels{p(pl)})
%                 axis tight
%             end
%             suptitle(['Colourbar scale from ' num2str(min(Results(:,9))) ' to ' num2str(max(Results(:,9)))])
            
                figure % FIG 4
                for pl = 1:2
                    subplot(1,2,pl)
                    p = [b e];
                    scatter(Results_tr(:,4), Results_tr(:,p(pl)), [], markerColors,'filled');
                    colormap(jet(256)); cb = colorbar;
                    xlabel(xlabels{4})
                    ylabel(xlabels{p(pl)})
                    axis tight
                end
                suptitle(['Colourbar scale from ' num2str(min(Results(:,9))) ' to ' num2str(max(Results(:,9))) ...
                    '; label= Proportion of initial decisions that are correct'])
            
        elseif outcome == 21
            figure % FIG 5
            for pl = 1:3
                subplot(1,3,pl)
                p = [b c e];
                scatter(Results_tr(:,1), Results_tr(:,p(pl)), [], markerColors,'filled');
                colormap(jet(256)); cb = colorbar;
                xlabel(xlabels{1})
                ylabel(xlabels{p(pl)})
                axis tight
            end
            suptitle(['Colourbar scale from ' num2str(min(Results(:,21))) ' to ' num2str(max_scale) ...
                '; label= Delusion score'])
            
        end
    end
    
    %% Get proportions with false inferences/delusions
    Threshold            = 0.66;
    Delusn_thresh_trust  = Threshold; 
    Delusn_thresh_symm   = Threshold/2; % threshold halved because delusion will always result in 50% correct
    Delusn_conf          = Threshold;
    Delusn_const         = Threshold;
    Card_habit_threshold = Threshold/2;    % = >0.66 prob of choosing one/other
    Entropy_threshold    = 0.8113;  % -sum([0.25 0.75].*log2([0.25 0.75]))
    cd(save_path) 
    [~,curr_dir,~]       = fileparts(pwd);
    
    if contains(curr_dir,'trust')
        Delusn_FI_thresh = Delusn_thresh_trust;
    elseif contains(curr_dir,'symm')
        Delusn_FI_thresh = Delusn_thresh_symm;
    end
        
    Delusn_subj          = Results(:,7)>Delusn_FI_thresh & ...
                           Results(:,12)>Delusn_conf & ...
                           Results(:,13)>Delusn_const;
    
    % Plot delusions within false inferences
    figure
    plot3(Results(~Delusn_subj,7),Results(~Delusn_subj,12),Results(~Delusn_subj,13),'ko')
    hold on
    plot3(Results(Delusn_subj,7),Results(Delusn_subj,12),Results(Delusn_subj,13),'ro')
    xlim([0 1]); ylim([0 1]); zlim([0 1])
    xlabel('False as proportion of total inferences')
    ylabel('Proportion of false infs of >90% confidence')
    zlabel('Consistency of false inferences')
    grid on
    suptitle(strrep(curr_dir,'_',' '))
    
    FI = sum(Results(:,7)>0);     % false inf
    TI = sum(Results(:,7)==0);    % true inf
    Del = sum(Delusn_subj);       % delusions
    NDe = sum(~Delusn_subj);      % non-delusions
    FIt = sum(Results(:,7)>=Delusn_FI_thresh);      % delusions
    TIt = sum(Results(:,7)<Delusn_FI_thresh);       % non-delusions
    CHb = sum(Results(:,10)>=Card_habit_threshold); % card habit
    nCH = sum(Results(:,10)<Card_habit_threshold);  % no card habit
    try
    if contains(curr_dir,'trust')
        FIt(2) = sum(sum(Results(:,16:17),2)>=0.5);      % delusions
        TIt(2) = sum(sum(Results(:,16:17),2)<0.5);       % non-delusions
    else
        FIt(2) = sum(Results(:,19)<=Entropy_threshold);      % delusions
        TIt(2) = sum(Results(:,19)>Entropy_threshold);       % non-delusions
    end
    CHb(2) = sum(Results(:,18)<=Entropy_threshold);  % card habit
    nCH(2) = sum(Results(:,18)>Entropy_threshold);   % no card habit
    end
    
    switch curr_dir
        case 'MDP_sim_trust_beta_NoAf'
            MDP_sim_trust_beta_NoAf.FI = FI;
            MDP_sim_trust_beta_NoAf.TI = TI;
            MDP_sim_trust_beta_NoAf.Del = Del;
            MDP_sim_trust_beta_NoAf.NDe = NDe;
            MDP_sim_trust_beta_NoAf.FIt = FIt;
            MDP_sim_trust_beta_NoAf.TIt = TIt;
            MDP_sim_trust_beta_NoAf.CHb = CHb;
            MDP_sim_trust_beta_NoAf.nCH = nCH
        case 'MDP_sim_trust_beta'
            MDP_sim_trust_beta.FI = FI;
            MDP_sim_trust_beta.TI = TI;
            MDP_sim_trust_beta.Del = Del;
            MDP_sim_trust_beta.NDe = NDe;
            MDP_sim_trust_beta.FIt = FIt;
            MDP_sim_trust_beta.TIt = TIt;
            MDP_sim_trust_beta.CHb = CHb;
            MDP_sim_trust_beta.nCH = nCH
        case 'MDP_sim_symm_beta_NoAf'
            MDP_sim_symm_beta_NoAf.FI = FI;
            MDP_sim_symm_beta_NoAf.TI = TI;
            MDP_sim_symm_beta_NoAf.Del = Del;
            MDP_sim_symm_beta_NoAf.NDe = NDe;
            MDP_sim_symm_beta_NoAf.FIt = FIt;
            MDP_sim_symm_beta_NoAf.TIt = TIt;
            MDP_sim_symm_beta_NoAf.CHb = CHb;
            MDP_sim_symm_beta_NoAf.nCH = nCH
        case 'MDP_sim_symm_beta'
            MDP_sim_symm_beta.FI = FI;
            MDP_sim_symm_beta.TI = TI;
            MDP_sim_symm_beta.Del = Del;
            MDP_sim_symm_beta.NDe = NDe;
            MDP_sim_symm_beta.FIt = FIt;
            MDP_sim_symm_beta.TIt = TIt;
            MDP_sim_symm_beta.CHb = CHb;
            MDP_sim_symm_beta.nCH = nCH
        case 'MDP_sim_symm_beta_treat05_10' % = 0.5
            MDP_sim_symm_beta_treat05_10.FI = FI;
            MDP_sim_symm_beta_treat05_10.TI = TI;
            MDP_sim_symm_beta_treat05_10.Del = Del;
            MDP_sim_symm_beta_treat05_10.NDe = NDe;
            MDP_sim_symm_beta_treat05_10.FIt = FIt;
            MDP_sim_symm_beta_treat05_10.TIt = TIt;
            MDP_sim_symm_beta_treat05_10.CHb = CHb;
            MDP_sim_symm_beta_treat05_10.nCH = nCH
        case 'MDP_sim_symm_beta_treat05_01' % = 0.5
            MDP_sim_symm_beta_treat05_01.FI = FI;
            MDP_sim_symm_beta_treat05_01.TI = TI;
            MDP_sim_symm_beta_treat05_01.Del = Del;
            MDP_sim_symm_beta_treat05_01.NDe = NDe;
            MDP_sim_symm_beta_treat05_01.FIt = FIt;
            MDP_sim_symm_beta_treat05_01.TIt = TIt;
            MDP_sim_symm_beta_treat05_01.CHb = CHb;
            MDP_sim_symm_beta_treat05_01.nCH = nCH
    end
    % Do chi squared tests
    % Z1 = [MDP_sim_symm_beta_NoAf.FI   MDP_sim_symm_beta_NoAf.TI;
    %       MDP_sim_symm_beta.FI        MDP_sim_symm_beta.TI         ];
    % [~,p1,X1] = chi2cont(Z1);
    %
    % Z2 = [MDP_sim_symm_beta_NoAf.Del  MDP_sim_symm_beta_NoAf.NDe;
    %       MDP_sim_symm_beta.Del       MDP_sim_symm_beta.NDe         ];
    % [~,p2,X2] = chi2cont(Z2);
    %
    % Z3 = [MDP_sim_symm_beta_NoAf.CHb  MDP_sim_symm_beta_NoAf.nCH;
    %       MDP_sim_symm_beta.CHb       MDP_sim_symm_beta.nCH         ];
    % [~,p3,X3] = chi2cont(Z3);
    
    % Plot state space trajectories
    % Load previous analysis to plot on same axes
    load([save_root 'MDP_sim_symm_beta/Results/MDP_sim_symm_beta.mat'])
    
    figure
    subplot(1,2,1)
    plot_state_space_traj(Trajectories,Results(:,7),1,1,0,MDP_sim_symm_beta); view(-67,45);
    subplot(1,2,2)
    all_traj = plot_state_space_traj(Trajectories,Results(:,7),1,1,0,MDP_sim_symm_beta); view(-200,10);
    try
        if strcmp(curr_dir,'MDP_sim_symm_beta') == 1
            MDP_sim_symm_beta.xlimits = xlim;
            MDP_sim_symm_beta.ylimits = ylim;
            MDP_sim_symm_beta.zlimits = zlim;
        end
    end
    % Plot mean trajectories
    all_traj     = sortrows(all_traj,'ascend');
    num_per_mean = floor(length(Results)/10);
    num_trials   = length(Trajectories(1).e_traj);
    t            = num_trials-1;
    for m = 1:10
        subj_in_mean = [(m-1)*num_per_mean+1 m*num_per_mean];
        mean_traj(m).e_traj          = mean(all_traj(subj_in_mean(1):subj_in_mean(2),2:t+1),1);
        mean_traj(m).post_trust_traj = mean(all_traj(subj_in_mean(1):subj_in_mean(2),t+2:t*2+1),1);
        mean_traj(m).pol_prec_traj   = mean(all_traj(subj_in_mean(1):subj_in_mean(2),t*2+2:t*3+1),1);
        mean_Results(m,1)            = mean(all_traj(subj_in_mean(1):subj_in_mean(2),1),1);
    end
    figure
    subplot(1,2,1)
    dummy = plot_state_space_traj(mean_traj,mean_Results,0,5,0,MDP_sim_symm_beta); view(-67,45);
    subplot(1,2,2)
    dummy = plot_state_space_traj(mean_traj,mean_Results,0,5,0,MDP_sim_symm_beta); view(-200,10);
    try
        if strcmp(curr_dir,'MDP_sim_symm_beta') == 1
            MDP_sim_symm_beta.mean_xlimits = xlim;
            MDP_sim_symm_beta.mean_ylimits = ylim;
            MDP_sim_symm_beta.mean_zlimits = zlim;
        end
    end
    % Adjust view to be consistent?
    % [a,e] = view(gca) % get azimuth and elevation of repositioned plot
    % view(a,e)         % apply to different plot (select it first)
    % Animate and save?
    % cd([save_path 'Results'])
    % Options.FrameRate=30;Options.Duration=12;Options.Periodic=false;
    % CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],[curr_dir '_movie'],Options)


% Plot relative amounts of false inferences re advisor/cards
v = 2; % which version (1/2) of delusions re trust/cards to use? 
try
    Data(1,1) = MDP_sim_trust_beta_NoAf.FIt(v)/(MDP_sim_trust_beta_NoAf.FIt(v)+MDP_sim_trust_beta_NoAf.TIt(v));
    Data(1,2) = MDP_sim_trust_beta_NoAf.CHb(v)/(MDP_sim_trust_beta_NoAf.CHb(v)+MDP_sim_trust_beta_NoAf.nCH(v));
    Data(2,1) = MDP_sim_trust_beta.FIt(v)/(MDP_sim_trust_beta.FIt(v)+MDP_sim_trust_beta.TIt(v));
    Data(2,2) = MDP_sim_trust_beta.CHb(v)/(MDP_sim_trust_beta.CHb(v)+MDP_sim_trust_beta.nCH(v));
    Data(3,1) = MDP_sim_symm_beta_NoAf.FIt(v)/(MDP_sim_symm_beta_NoAf.FIt(v)+MDP_sim_symm_beta_NoAf.TIt(v));
    Data(3,2) = MDP_sim_symm_beta_NoAf.CHb(v)/(MDP_sim_symm_beta_NoAf.CHb(v)+MDP_sim_symm_beta_NoAf.nCH(v));
    Data(4,1) = MDP_sim_symm_beta.FIt(v)/(MDP_sim_symm_beta.FIt(v)+MDP_sim_symm_beta.TIt(v));
    Data(4,2) = MDP_sim_symm_beta.CHb(v)/(MDP_sim_symm_beta.CHb(v)+MDP_sim_symm_beta.nCH(v));
    Data(6,1) = MDP_sim_symm_beta_treat05_10.FIt(v)/(MDP_sim_symm_beta.FIt(v)+MDP_sim_symm_beta.TIt(v));
    Data(6,2) = MDP_sim_symm_beta_treat05_10.CHb(v)/(MDP_sim_symm_beta.CHb(v)+MDP_sim_symm_beta.nCH(v));
    
    h = findobj('type','figure');
    n = length(h);
    figure(n+1)
    
    bar(Data);
    row1 = {'Stable advice' 'Stable advice' 'Unstable advice' 'Unstable advice' '' 'Unstable advice'};
    row2 = {'No mood' 'With mood' 'No mood' 'With mood' '' 'With mood'};
    row3 = {'' '' '' '' '' 'TREATED'};
    labelArray = [row1; row2; row3];
    labelArray = strjust(pad(labelArray),'center'); % 'left'(default)|'right'|'center
    tickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\n', labelArray{:}));
    % Assign ticks and labels
    ax = gca(figure(n+1));
    ax.XTickLabel = tickLabels;
    legend('Trusting advisor','Choosing cards')
    ylabel('Frequency of delusion-like behaviour')
end

end

if doSave
    % Save Results structure?
    cd([save_path 'Results'])
    save('Results','Results')
    save('Trajectories','Trajectories')
    try
        save(curr_dir,curr_dir)
    end
end

% Check plots
% s = 1; % number of subject's row in Results
% plot_trial(mdp)
% suptitle(['Fraction False Inf = ' num2str(Results(s,7)) ', Proportion >90% confidence = ' ...
%     num2str(Results(s,12)) ', Proportion consistent = ' num2str(Results(s,13))])

% figure
% for p = 1:5
% subplot(1,5,p)
% scatter(Results(:,p),Results(:,7),'kx')
% lsline
% xlabel(xlabels{p}); ylabel('False inferences')
% end

% Transform habit resistance to find interaction with initial consistency
% habitRes_tr = abs(Results(:,2)-max(Results(:,2)));
% interaction = habitRes_tr.*Results(:,3);
% figure; subplot(1,3,1)
% scatter(interaction,Results(:,7),'kx')
% lsline
% corr(interaction,Results(:,7),'type','Spearman')
% xlabel('habit*consistency'); ylabel('False inferences')