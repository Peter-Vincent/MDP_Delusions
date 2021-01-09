# Documentation for figure 2

## Figure 2 summary
Panels in this figure illustrate the benefits and drawbacks of habits in behaviour, depending on the statistics of the environment.  Titles to the figures are provided in the file names, and these match the .mat data files, allowing for re-construction of the figure from the same data.  

### General interpretation
These figures show the most basic results of a set of simulations.  THe black line through the middle of the figure shows a "trajectory" of whether the agent did or did not trust the advisor.  If they did trust the advisor, the value represented by the black line increases by 1.  If they did not trust the advisor, the value decreases by 1.  The top row (Advised card) shows the advise given in each trial, using the colours given in the text (blue and green).  THe section of the middle row, on top of the black trajectory, shows the choices made, again in the green and blue of the cards, as described in the text.  The middle row below the black line shows the feedback recieved by the agent.  Violet representes negative feedback whilst dark purple represents positive feedback.  The bottom panel represents the state of the advisor, where light red represents un-trustworthy and dark red represents trustworthy.  THerefore, the central strip represents the adctions are observations of the agent, whilst the top and bottom strips represent the true state of the world.

### trust_sequence_no_habits.png
This panel shows that no learning/forming of habits in a given environment leads to poor performance, even when sufficient information is given for performance to be high.  THe advisor is consistently trustworthy (with only the occasional fallacy), but because the agent does not learn this, they make random guesses about the card.
`output = MDP_Delusions(256,trust_sequence,1,0.9,600);` 

### trust_sequence_high_habits.png
This panel shows that learning/forming of habits in a given environment leads to good performance, since the environment is broadly stationary.  The agent easily forms a "habit" of trusting the advisor, and since the advisor is normally trustworthy (bottom row), the agent is normally correct.
`output = MDP_Delusions(256,trust_sequence,1,0.9,2);`

### switch_sequence_no_habits.png
This panel shows that no learning/forming of habits in a given environment leads to poor performance, even when sufficient information is given for performance to be high.  THe advisor is consistently trustworthy (with only the occasional fallacy), but because the agent does not learn this, they make random guesses about the card.  
This panel is superflous given "trust_sequence_no_habits.png".
`output = MDP_Delusions(256,switch_sequence,1,0.9,600);`

### switch_sequence_high_habits.png
This panel shows that learning/forming of habits when the environment undergoes a sudden change can lead to plenty of errors, but that if the habit is re-learnt performance can pick back up again.  See this in the change in advisor trustworthiness in the bottom strip, and the corresponding changes in light and dark purple in the middle panel under the black line.
`output = MDP_Delusions(256,switch_sequence,1,0.9,2);`

## Construction of these figures
These figures are made with either "trust_sequence.mat" or with "switch_sequence.mat".  The code to execute these is given at the bottom of each figure.  `plot_trial(output)` is used to generate the figure.  This code produces 6 seperate figures, the first is used in these reports.