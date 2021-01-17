% 28/06/2020 Peter Vincent

function plot_trial(mdp)

[decision_traj,feedback_traj,card_traj] = behaviour(mdp);
[advised_card,trustworthy] = advisor(mdp);

num_trials = length(mdp);

prior_belief    =  zeros(1,num_trials);
posterior_belief=  zeros(1,num_trials);
prior_card      =  zeros(1,num_trials);
posterior_card  =  zeros(1,num_trials);

correct_post_inf=  zeros(1,num_trials);
correct_card_inf=  zeros(1,num_trials);
affect_state    =  zeros(1,num_trials);
fig1 = figure('color','w');
max_pos = length(mdp);
max_total = round(max_pos * 1.5);
hold on
xlim([0.5 (num_trials+0.5)]);
ylim([-max_total max_total])
sep = 10;
%% Define the colours
card_1_col        = '#0100FF';
card_2_col        = '#70CF32';
trust_col         = '#E00000';
not_trust_col     = '#EDB9B9';
corr_feedback_col = '#6D0BBF';
false_feedback_col= '#D1B4E9';
posterior_col     = '#8B988C';
poor_inf_col      = '#36F3CB';

for trial = 1:num_trials
    
    cur_bounds = decision_traj(1,trial);
    %% Plot the behavior of the agent
    if feedback_traj(1,trial) == 1
        feedback_colour = corr_feedback_col;
    else
        feedback_colour = false_feedback_col;
    end
    rectangle('Position',[trial-0.5,-max_pos,1,cur_bounds+max_pos],...
        'FaceColor',feedback_colour,'EdgeColor','none');
    if card_traj(1,trial) == 1
        card_colour = card_1_col;
    else
        card_colour = card_2_col;
    end
    rectangle('Position',[trial-0.5,cur_bounds,1,max_pos-cur_bounds],...
        'FaceColor',card_colour,'EdgeColor','none');
    %% Plot the behaviour of the advisor
    if trustworthy(1,trial) == 1
        trustworthy_colour = trust_col;
    else
        trustworthy_colour = not_trust_col;
    end
    rectangle('Position',[trial-0.5,-max_total,1,max_total-max_pos-sep],...
        'FaceColor',trustworthy_colour,'EdgeColor','none');
    if advised_card(1,trial) == 1
        card_colour = card_1_col;
    else
        card_colour = card_2_col;
    end
    rectangle('Position',[trial-0.5,max_pos + sep,1,max_total-max_pos-sep],...
        'FaceColor',card_colour,'EdgeColor','none');
    
    
    %% Get the beliefs of the agent
    cur_trial = mdp(trial);
    outcomes   = cur_trial.o;
    inf_over_ad=cur_trial.xn{1,1};
    inf_over_card = cur_trial.xn{1,2};
    
    prior_card(1,trial) = inf_over_card(1,1,3,2);
    posterior_card(1,trial) = inf_over_card(1,1,3,3);
    
    card_choice = prior_card(1,trial) >= 0.5;
    final_card_choice = posterior_card(1,trial) >= 0.5;
    if (outcomes(2,3) == 1) && (card_choice == final_card_choice)
        correct_card_inf(1,trial) = 1;
    elseif (outcomes(2,3) == 2) && (card_choice ~= final_card_choice)
        correct_card_inf(1,trial) = 1;
    else
        correct_card_inf(1,trial) = 0;
    end
        
    
    
    prior_belief(1,trial) = inf_over_ad(end,1,2,2);
    posterior_belief(1,trial) = inf_over_ad(end,1,2,3);
    
    if (outcomes(1,2) == outcomes(4,3)) && (outcomes(2,3) == 1)
        trustworthy_val = 1;
    elseif (outcomes(1,2) ~= outcomes(4,3)) && (outcomes(2,3) == 1)
        trustworthy_val = 0;
    elseif (outcomes(1,2) == outcomes(4,3)) && (outcomes(2,3) == 2)
        trustworthy_val = 0;
    elseif (outcomes(1,2) ~= outcomes(4,3)) && (outcomes(2,3) == 2)
        trustworthy_val = 1;
    end
    if (trustworthy_val == 1) &&  (posterior_belief(1,trial) > 0.5)
        correct_post_inf(1,trial) = 1;
    elseif (trustworthy_val == 0) && (posterior_belief(1,trial) < 0.5)
        correct_post_inf(1,trial) = 1;
    else
        correct_post_inf(1,trial) = 0;
    end
    
    states = cur_trial.s;
    affect_state(1,trial) = states(4,2);
    
end
% Calculate performance metrics
per_correct = (sum(feedback_traj == 1)/length(feedback_traj)) * 100;
per_good_inf=(1-(sum(correct_post_inf == 1)/length(correct_post_inf)))*100;
per_good_card_inf = (1 - (sum(correct_card_inf == 1)/length(correct_card_inf)))*100;

plot(1:num_trials,decision_traj,'LineWidth',2,'Color','k');
ax = gca;
set(ax,'Box','off');
ax.XLabel.String = "Trial number";
ax.YAxis.Visible = 'off'; 
%% Add legends
traj_message = sprintf("Decision\n trajectory");
card_advice_message = sprintf("Advised\n card");
card_choice_message = sprintf("Chosen\n card");
feedback_message    = sprintf("Feedback\n recieved");
advisor_message     = sprintf("Advisor\n trustworthiness");
text(ax,0,0,traj_message,'HorizontalAlignment','right');
text(ax,0,max_pos+sep + 0.5*(max_total-max_pos-sep),card_advice_message,'HorizontalAlignment','right');
text(ax,0,max_pos/2,card_choice_message,'HorizontalAlignment','right');
text(ax,0,-max_pos/2,feedback_message,'HorizontalAlignment','right');
text(ax,0,-max_total+0.5*(max_total-max_pos-sep),advisor_message,'HorizontalAlignment','right');

text(ax,num_trials+1,-max_total+0.5*(max_total-max_pos-sep),'Trustworthy','Color',trust_col);
text(ax,num_trials+1,-max_pos/2,'Correct','Color',corr_feedback_col);
text(ax,num_trials+1,-max_pos,strcat(num2str(per_correct),'%'),'Color',corr_feedback_col);


%%%% Make the second figure to illustrate the beliefs of the agent with
%%%% respect to the trial (prior, posterior, ...
fig2 = figure('color','w');
hold on
xlim([0.5 (num_trials+0.5)]);
ylim([0 1])
ax2 = gca;
for trial = 1:num_trials
    if trustworthy(1,trial) == 1
        trustworthy_colour = trust_col;
    else
        trustworthy_colour = not_trust_col;
    end
    rectangle(ax2,'Position',[trial-0.5,0,1,2],...
        'FaceColor',trustworthy_colour,'EdgeColor','none');
    if correct_post_inf(1,trial) == 0
        rectangle(ax2,'Position',[trial-0.5,0.5,1,0.5],...
            'FaceColor',poor_inf_col,'EdgeColor','none');
    end
end
plot(ax2,1:num_trials,prior_belief,'color','k','LineWidth',2);
plot(ax2,1:num_trials,posterior_belief,'.','color',posterior_col,'MarkerSize',15);
posterior_message = sprintf('Posterior\nBelief');
prior_message     = sprintf('Prior\nBelief');
false_message     = sprintf('False\nInference');
text(ax2,num_trials+1,0.5,prior_message,'Color','k');
text(ax2,num_trials+1,0.25,strcat(num2str(per_good_inf),'%'),'Color',poor_inf_col);

text(ax2,num_trials+1,1,posterior_message,'Color',posterior_col);
text(ax2,num_trials+1,0,false_message,'Color',poor_inf_col);
ax2.XLabel.String = "Trial number";
ax2.YLabel.String = "Probability";
hold('off');
%%
%%%% Make the third figure to illustrate the affective state of the agent
%%%% with smoothing across some number of trials (10 in this case)
smoothing_kernel = 10;
plot_range = smoothing_kernel:(length(affect_state)-smoothing_kernel);
fig3 = figure('color','w');
xlim([0.5 (num_trials+0.5)]);
ylim([0.8 2.2]);
ax3 = gca;
hold(ax3,'on');
smooth = movmean(affect_state,smoothing_kernel);
affect_state_low = find(affect_state == 1);
affect_state_high= find(affect_state == 2);
scatter(ax3,affect_state_low,affect_state(affect_state_low),25,'filled','MarkerFaceColor',trust_col);
scatter(ax3,affect_state_high,affect_state(affect_state_high),25,'filled','MarkerFaceColor',not_trust_col);
plot(ax3,plot_range,1.5*ones(1,length(plot_range)),'LineWidth',1,'color','k');
plot(ax3,plot_range,smooth(plot_range),'LineWidth',2,'color','k');
ax3.XLabel.String = "Trial number";
ax3.YAxis.Visible = 'off'; 
text(ax3,0,1.5,"Arousal state",'HorizontalAlignment','right');
text(ax3,num_trials+1,2,"High arousal",'Color','k');
text(ax3,num_trials+1,1,"Low arousal",'Color','k');
hold(ax3,'off');

%% Make figure to plot the inferences over the card

fig4 = figure('color','w');
hold on
xlim([0.5 (num_trials+0.5)]);
ylim([0 1])
ax4 = gca;
for trial = 1:num_trials
    if feedback_traj(1,trial) == 1
        feedback_colour = corr_feedback_col;    
    else
        feedback_colour = false_feedback_col;
    end
    rectangle(ax4,'Position',[trial-0.5,0,1,2],...
        'FaceColor',feedback_colour,'EdgeColor','none');
    if correct_card_inf(1,trial) == 0
        rectangle(ax4,'Position',[trial-0.5,0.5,1,0.5],...
            'FaceColor',poor_inf_col,'EdgeColor','none');
    end
end

plot(ax4,1:num_trials,prior_card,'color','k','LineWidth',2);
plot(ax4,1:num_trials,posterior_card,'.','color',posterior_col,'MarkerSize',15);
%posterior_message = sprintf('Posterior\nBelief');
%prior_message     = sprintf('Prior\nBelief');
%text(ax4,num_trials+1,0.5,prior_message,'Color','k');
text(ax4,num_trials+1,0.25,strcat(num2str(per_good_card_inf),'%'),'Color',poor_inf_col);

%text(ax4,num_trials+1,1,posterior_message,'Color',posterior_col);
ax4.XLabel.String = "Trial number";
ax4.YLabel.String = "Probability";
hold('off');

%% Make the combined figure without arousal
fig1_axes = findobj('Parent',fig1,'Type','axes');
fig2_axes = findobj('Parent',fig2,'Type','axes');
fig3_axes = findobj('Parent',fig4,'Type','axes');
h1_ax = fig1_axes(1);
h2_ax = fig2_axes(1);
h3_ax = fig3_axes(1);
figure('Color','w');
sub_ax1 = subplot(5,1,1:3);
set(sub_ax1,'Box','off');
sub_ax1.YAxis.Visible = 'off'; 
sub_ax1.XAxis.Visible = 'off';
sub_ax1.XLim = [0.5 (num_trials+0.5)];
sub_ax1.YLim = [-max_total max_total];

sub_ax2 = subplot(5,1,4:5);
sub_ax2.XLim = [0.5 (num_trials+0.5)];
sub_ax2.YLim = [0 1];
sub_ax2.YLabel.String = "Probability";


copyobj(get(h1_ax,'Children'),sub_ax1);
copyobj(get(h2_ax,'Children'),sub_ax2);
%% Make the combined figure with arousal
fig3_axes = findobj('Parent',fig3,'Type','axes');
h4_ax = fig3_axes(1);
figure('Color','w');

sub_ax1 = subplot(6,1,1:3);
set(sub_ax1,'Box','off');
sub_ax1.YAxis.Visible = 'off'; 
sub_ax1.XAxis.Visible = 'off';
sub_ax1.XLim = [0.5 (num_trials+0.5)];
sub_ax1.YLim = [-max_total max_total];

sub_ax2 = subplot(6,1,4);
set(sub_ax2,'Box','off');
sub_ax2.YAxis.Visible = 'off'; 
sub_ax2.XAxis.Visible = 'off';
sub_ax2.XLim = [0.5 (num_trials+0.5)];
sub_ax2.YLim = [0.8 2.2];

sub_ax3 = subplot(6,1,5:6);
sub_ax3.XLim = [0.5 (num_trials+0.5)];
sub_ax3.YLim = [0 1];
sub_ax3.XAxis.Visible = 'off';
sub_ax3.YLabel.String = "Probability";

% sub_ax4 = subplot(8,1,7:8);
% sub_ax4.XLim = [0.5 (num_trials + 0.5)];
% sub_ax4.YLim = [0 1];
% sub_ax4.XLabel.String = "Trial number";
% sub_ax4.YLabel.String = "Probability";

copyobj(get(h1_ax,'Children'),sub_ax1);
copyobj(get(h4_ax,'Children'),sub_ax2);
copyobj(get(h2_ax,'Children'),sub_ax3);
% copyobj(get(h3_ax,'Children'),sub_ax4);