% Peter Vincent 01/11/2020

function plot_single_trial(mdp_trial)

plot_width = 300;
box_width  = 80;
action     = 20;
gaps       = 5;

card_1_col        = '#0100FF';
card_2_col        = '#70CF32';
trust_col         = '#E00000';
not_trust_col     = '#EDB9B9';
corr_feedback_col = '#6D0BBF';
false_feedback_col= '#D1B4E9';
posterior_col     = '#8B988C';
poor_inf_col      = '#36F3CB';
null_state        = [0.1 0.1 0.1 0.3];
inf_colour        = [1 1 1];


inf_on_ad = mdp_trial.xn{1,1};
outcomes  = mdp_trial.o;
states = mdp_trial.s;
card_advice = outcomes(1,2);
card_selected = outcomes(4,3);
feedback = outcomes(2,3);
trustworthy = states(1,2);
pre_trial_belief = inf_on_ad(:,1,2,1);
prior_belief     = inf_on_ad(:,1,2,2);
posterior_belief = inf_on_ad(:,1,2,3);
if (outcomes(1,2) == outcomes(4,3)) && (outcomes(2,3) == 1)
    trustworthy_val = 1;
elseif (outcomes(1,2) ~= outcomes(4,3)) && (outcomes(2,3) == 1)
    trustworthy_val = 0;
elseif (outcomes(1,2) == outcomes(4,3)) && (outcomes(2,3) == 2)
    trustworthy_val = 0;
elseif (outcomes(1,2) ~= outcomes(4,3)) && (outcomes(2,3) == 2)
    trustworthy_val = 1;
end
if (trustworthy_val == 1) &&  (posterior_belief(end) > 0.5)
    correct_post_inf = 1;
elseif (trustworthy_val == 0) && (posterior_belief(end) < 0.5)
    correct_post_inf = 1;
else
    correct_post_inf = 0;
end



fig1 = figure('color','w');
ax = gca;
hold(ax,'on');
cur_ref = 1;
rectangle('Position',[cur_ref,0,box_width-1,1],...
    'FaceColor',null_state,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
plot(cur_ref:(cur_ref+box_width-1),interp1(1:length(pre_trial_belief),...
    pre_trial_belief,linspace(1,length(pre_trial_belief),box_width)),...
    'color',inf_colour,'LineWidth',2);

ylim([0 1]);
xlim([1 plot_width+1]);
set(ax,'Box','off');
ax.XAxis.Visible = 'off';
ax.YLabel.String = "Probability";

text(ax,cur_ref,-0.03,"Initial Stage",'HorizontalAlignment','left','FontSize',20);
%%%%%%%%%%%%%%%%
cur_ref = cur_ref + box_width + gaps + action + gaps;

if card_advice == 1
rectangle('Position',[cur_ref,0.5,box_width-1,0.5],...
    'FaceColor',card_1_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
else
    rectangle('Position',[cur_ref,0.5,box_width-1,0.5],...
    'FaceColor',card_2_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
end

if trustworthy == 1
rectangle('Position',[cur_ref,0,box_width-1,0.5],...
    'FaceColor',trust_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
else
    rectangle('Position',[cur_ref,0,box_width-1,0.5],...
    'FaceColor',not_trust_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
end

plot(cur_ref:(cur_ref+box_width-1),interp1(1:length(prior_belief),...
    prior_belief,linspace(1,length(prior_belief),box_width)),...
    'color',inf_colour,'LineWidth',2);

text(ax,cur_ref,-0.03,"Advice Stage",'HorizontalAlignment','left','FontSize',20);

%%%%%%%%%%%%%%
cur_ref = cur_ref + box_width + gaps;


if card_selected == 1
rectangle('Position',[cur_ref,0,action,1],...
    'FaceColor',card_1_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
else
    rectangle('Position',[cur_ref,0,action,1],...
    'FaceColor',card_2_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
end
text(ax,cur_ref,-0.03,"Action",'HorizontalAlignment','left','FontSize',20);


%%%%%%%%%%%
cur_ref = cur_ref + action + gaps;

if feedback == 1
rectangle('Position',[cur_ref,0,box_width-1,1],...
    'FaceColor',corr_feedback_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
else
    rectangle('Position',[cur_ref,0,box_width-1,1],...
    'FaceColor',false_feedback_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
end

if correct_post_inf == 0
    rectangle('Position',[cur_ref,0.5,box_width-1,0.5],...
    'FaceColor',poor_inf_col,'EdgeColor','k','LineWidth',2,...
    'Curvature',0.2);
end

plot(cur_ref:(cur_ref+box_width-1),interp1(1:length(posterior_belief),...
    posterior_belief,linspace(1,length(posterior_belief),box_width)),...
    'color',posterior_col,'LineWidth',2);

text(ax,cur_ref,-0.03,"Feedback Stage",'HorizontalAlignment','left','FontSize',20);









