# Documentation for figure 3

## Figure 3 summary
Panels in this figure illustrate the impact of low precision in the agents beliefs about the environment.  Titles to the figures are provided in the file names, and these match the .mat data files, allowing for re-construction of the figure from the same data.

### General interpretation
The interpretation of the top panels of each figure is the same as in figure 2 (see docs in Fig2 folder).  These figures, however, contain an additional panel at the bottom.  This additional panel shows the prior precision over advisor trustworthiness in the black line.  The grey dots give the posterior precision over advisor trustworthiness.  The colour of the bars (red or pale pink) give the trsutworthiness of the advisor.  In some trials the top half of a given trial is coloured in cyan.  This indicates that the agent experienced a "delusion" in that particular trial.  A delusion is defined here as a posterior belief about the trustworthiness of the advisor that is inconsistent with the observed evidence.

### switch_sequence_habits_06_prec.png
This figure shows that low precision coupled with habit forming behaviours leads to delusions about the advisor.
`output = MDP_Delusions(256,switch_sequence,1,0.6,2);`

### switch_sequence_no_habits_06_prec.png
This figure shows that even with low precision, if there is no habit forming behaviour there are no delusions about the advisor.
`output = MDP_Delusions(256,switch_sequence,1,0.6,600);`

### switch_sequence_habits_075_prec.png
This figure shows delusions occur even when low precision beleifs aren't at the extreme end of the spectrum.  One can also compare the habit forming behaviour between this simulation and when beliefs about precision are extrmemly low.
`output = MDP_Delusions(256,switch_sequence,1,0.75,2);`

### switch_sequence_habits_09_prec.png
This figure shows that low precision is neccessary for delusion forming behaviours, if with strong habit forming behaviour.  This figure is a duplicate of **switch_sequence_high_habits.png** in figure 2.
`output = MDP_Delusions(256,switch_sequence,1,0.9,2);`

## Construction of these figures
These figures are made with "switch_sequence.mat".  The code to execute these is given at the bottom of each figure.  `plot_trial(output)` is used to generate the figure.  This code produces 6 seperate figures, the fifth is used in these reports.