% 12/07/2020 Peter Vincent

function check_inferences(mdp)

num_trials = length(mdp);

pre_action_inf = zeros(1,num_trials);
post_feedback_inf=zeros(1,num_trials);
correct_post_inf = zeros(1,num_trials);


for trial = 1:num_trials
    cur_trial = mdp(trial);
    outcomes   = cur_trial.o;
    inf_over_ad=cur_trial.xn{1,1};
    pre_action_inf(1,trial) = inf_over_ad(end,1,2,2);
    post_feedback_inf(1,trial) = inf_over_ad(end,1,2,3);
    if (outcomes(1,2) == outcomes(4,3)) && (outcomes(2,3) == 1)
        trustworthy = 1;
    elseif (outcomes(1,2) ~= outcomes(4,3)) && (outcomes(2,3) == 1)
        trustworthy = 0;
    elseif (outcomes(1,2) == outcomes(4,3)) && (outcomes(2,3) == 2)
        trustworthy = 0;
    elseif (outcomes(1,2) ~= outcomes(4,3)) && (outcomes(2,3) == 2)
        trustworthy = 1;
    end
    if (trustworthy == 1) &&  (post_feedback_inf(1,trial) > 0.5)
        correct_post_inf(1,trial) = 1;
    elseif (trustworthy == 0) && (post_feedback_inf(1,trial) < 0.5)
        correct_post_inf(1,trial) = 1;
    else
        correct_post_inf(1,trial) = 0;
    end   
end


% Now generate plot

figure('color','w')
ax = gca;
hold(ax,'on')
pre_index = 1:2:(2*num_trials);
post_index= 2:2:(2*num_trials+1);
for trial = 1:num_trials
    plot(ax,[pre_index(trial) post_index(trial)],...
        [pre_action_inf(trial) post_feedback_inf(trial)],...
        '-k','Linewidth',2)
    if correct_post_inf(1,trial) == 1
        col = 'g';
    else
        col = 'r';
    end
    plot(ax,[(pre_index(trial)-0.5) (post_index(trial)+0.5)],[1 1],...
        'color',col,'Linewidth',2)
end
scatter(ax,pre_index,pre_action_inf,30,'filled','MarkerFaceColor','m');
scatter(ax,post_index,post_feedback_inf,30,'filled','MarkerFaceColor','b');
num_false_inf = sum(correct_post_inf)/length(correct_post_inf);
ax.Title.String = strcat("% correct inference = ",num2str(num_false_inf*100));
ax.YLim = [0 1];
    

    
    