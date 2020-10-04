% 26/04/2020 Peter Vincent
% This function parses the MDP structure to generate summary plots

function plot_comparison(mdp)
feature('DefaultCharacterSet','UTF-8');
timestep   = 3; % CHANGE THIS TO CHANGE TIMESTEP BEING PLOTTED
num_trials = length(mdp);
num_actions = num_trials;
process_blue_correct = zeros(1,num_actions);
infer_blue_correct   = zeros(1,num_actions);
process_advisor_trustworthy = zeros(1,num_actions);
infer_advisor_trustworthy = zeros(1,num_actions);
chosen_colour = zeros(1,num_actions);
feedback = zeros(1,num_actions);
arousal  = zeros(1,num_actions);
advice_followed = zeros(1,num_actions);
advice_given    = zeros(1,num_actions);
advice_correct  = zeros(1,num_actions);
for trial = 1:num_trials
    cur_trial = mdp(trial);
    outcomes  = cur_trial.o;
    states    = cur_trial.s;
    posteriors= cur_trial.xn;
    action1   = trial;
    feedback(1,action1) = outcomes(2,3);
    arousal(1,action1)  = states(4,2);
    chosen_colour(1,action1) = states(3,3);
    % Store the advice given
    advice_given(1,action1) = outcomes(1,2);
    % Store whether the advice is followed
    if outcomes(1,2) == states(3,3)
        advice_followed(1,action1) = 1;
    else
        advice_followed(1,action1) = 2; 
    end
    % Store whether the advice was correct (from subject's POV)
    if advice_followed(1,action1) == feedback(1,action1) 
        advice_correct(1,action1) = 1;
    else
        advice_correct(1,action1) = 0;
    end
    % Store whether the advice was correct (from advisor's POV)
%     if advice_given(1,action1) == states(2,1) 
%         advice_correct(1,action1) = 1;
%     else
%         advice_correct(1,action1) = 0;
%     end
    
    % Store probability of advisor being trustworthy
    A_mat = cur_trial.A{1,1};
    advice_mat_trust_state = A_mat(:,:,1,1,1,2);
    if states(1,2) == 1
        process_advisor_trustworthy(1,action1) = advice_mat_trust_state(1,1);
    else
        process_advisor_trustworthy(1,action1) = advice_mat_trust_state(2,1);
    end
    % Store inference of advisor being trustworthy
    infer_advisor_trustworthy(1,action1) = posteriors{1,1}(end,1,2,timestep); %%%% Change timestep at which inferences are made with final index
    % Store which card is good
    A_mat = cur_trial.A{1,2};
    correct_card_blue_state = A_mat(:,:,1,1,1,3);
    if states(2,1) == 1
        process_blue_correct(1,action1) = correct_card_blue_state(1,1);
    else
        process_blue_correct(1,action1) = correct_card_blue_state(2,1);
    end
    % Store inferences of blue being the correct card
    infer_blue_correct(1,action1) = posteriors{1,2}(end,1,2,timestep); %%%% Change timestep at which inferences are made with final index
    % Store whether outcome matched contingency for red (1) or blue (0)
    % (i.e. if choice == outcome (i.e. red1=correct1, blue2=incorrect2)
    if states(3,3) == outcomes(2,3)
        outcomeVcontingency(1,action1) = 1;
    % (or if choice ~= outcome (i.e. red1=incorrect2, blue2=correct1)
    elseif states(3,3) ~= outcomes(2,3)
        outcomeVcontingency(1,action1) = 0;
    end
end

index_array = 1:num_trials;

%% Top Figure

figure('color','w')
ax1 = subplot(3,1,1);
xlim([0 (num_trials + 1)])
ylim([-0.5 1.5])
hold(ax1,'on')
plot(index_array,process_blue_correct,'g','LineWidth',2);
plot(index_array,infer_blue_correct,'color',[0.9 0.4 0.2],'LineWidth',2);
choiceA = chosen_colour == 1;
choiceB = chosen_colour == 2;
choiceA_ind = index_array(choiceA);
choiceB_ind = index_array(choiceB);
scatter(choiceA_ind,ones(1,sum(choiceA))+0.3,20,'r','filled','s');
scatter(choiceB_ind,zeros(1,sum(choiceB))-0.3,20,'b','filled','s');
for choice_ind = 1:length(index_array)
    if feedback(choice_ind) == 1
        text_type = char(hex2dec('2713'));%'?';
    else
        text_type = 'X';
    end
    if outcomeVcontingency(choice_ind) == 1 % chosen_colour(choice_ind) == 1
        text_loc_y = 1;
    else
        text_loc_y = 0;
    end
    text_loc_x = index_array(choice_ind)-0.5; % shift text to align it correctly
    text(text_loc_x,text_loc_y,text_type);
end
set(gca,'XColor','none')
set(gca,'ytick',[0 0.5 1])
arrayfun( @(i) line([i i],[0 1],'Color','k','LineStyle','-'),10:10:num_trials,'UniformOutput',0)
ylabel('p(Red deck = correct)')
title(['Choice (R/B), Correct? (' char(hex2dec('2713')) '/X), With/Vs contingency (on/off green), Process (green), Inference (orange)'])
hold(ax1,'off')   
%% Middle Figure
ax2 = subplot(3,1,2);
xlim([0 (num_trials + 1)])
ylim([-0.5 1.5])
hold(ax2,'on')
plot(index_array,process_advisor_trustworthy,'g','LineWidth',2);
plot(index_array,infer_advisor_trustworthy,'color',[0.9 0.4 0.2],'LineWidth',2)
for advice_ind = 1:length(index_array)
    text_loc_x = index_array(advice_ind);
    if advice_followed(advice_ind) == 1
        text_type = char(hex2dec('2713'));%'?';
        text_loc_y= advice_correct(advice_ind);%1;
        advice_loc = 1.3;
    else
        text_type = 'X';
        text_loc_y= advice_correct(advice_ind);%0;
        advice_loc = 1.3;%-0.5;
    end
    
    text(text_loc_x-0.5,text_loc_y,text_type); % shift text to align it correctly
    if advice_given(advice_ind) == 1
        scatter(text_loc_x,advice_loc,20,'r','filled','s')
    else
        scatter(text_loc_x,advice_loc,20,'b','filled','s')
    end
end
set(gca,'XColor','none')
set(gca,'ytick',[0 0.5 1])
arrayfun( @(i) line([i i],[0 1],'Color','k','LineStyle','-'),10:10:num_trials,'UniformOutput',0)
ylabel('p(Advisor = trustworthy)')
title(['Advice (R/B), Trusted? (' char(hex2dec('2713')) '/X), With/Vs contingency (on/off green), Process (green), Inference (orange)'])
hold(ax2,'off') 
%% Bottom plot
ax3 = subplot(3,1,3);
xlim([0 (num_trials + 1)])
plot(index_array,arousal,'k','LineWidth',2); box off
set(gca,'ytick',[1 2],'yticklabels',{'Relaxed','Aroused'}); xlabel('Trial number')
arrayfun( @(i) line([i i],[1 2],'Color','k','LineStyle','-'),10:10:num_trials,'UniformOutput',0)
xlim([0 (num_trials + 1)])
% adjust advice_followed to allow comparison in title below
advice_followed(advice_followed==2) = 0;
title(['Overall: ' num2str(round(100*sum(feedback==1)/num_trials)) '% choices correct; ' ...
    num2str(round(100*sum(advice_followed==process_advisor_trustworthy)/num_trials)) '% advice appraised correctly; ' ...
    'timestep ' num2str(timestep)])
ylim([0.5 2.5])
hold(ax3,'off')
end
