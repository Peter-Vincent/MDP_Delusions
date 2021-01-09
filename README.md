# MDP_Delusions

This repository contains data and code relating to the MDP_Delusions project, a collaboration between Rick Adams, Thomas Parr, Peter Vincent, David Benrimoh and Karl Friston

## Project summary

The aim of the project is to provide a framework for simulating certain types of delusions in patients with paranoid schizophrenia.  The main contribution of this work is to provide an active inference based model to understand why delusions might occur in certain circumstances.  In particular, we focus on the following observations
- Habits can be beneficial in environments with regular statistics, and are flexible to the exact nature of those statistics
- Uncertainty over the reliability of the environment can lead imprecise posteriors about the environment
- Imprecise posteriors coupled with very precise priors can lead to delusions.  The subject of these delusions is driven by the subject of the strong priors
- Affect can be beneficial in environments with regular statistics, but are in-flexible to the exact nature of those statistics
- Uncertainty over posteriors is essential to form delusions

## Figure summary
### Figure 1
Figure 1 will show a schematic of the overall experimental structure, generative model and generative process, and provide a (clear) example of a single trial in a couple of different scenarios
- [ ] Basic figure draft
- [ ] Final figure complete

### Figure 2
Figure 2 will show how habits can be beneficial in environments with regular statistics, but can lead to inflexibility if those statistics suddenly change and habits are over-learnt
- [ ] Basic figure draft
- [ ] Final figure complete

### Figure 3
Figure 3 will show how imprecision in beliefs about the environment can couple with very strong priors to give rise to delusions, and will explore the nature of these delusions
- [ ] Basic figure draft
- [ ] Final figure complete

### Figure 4
Figure 4 will show a parameters sweep of precision in beleifs about feedback and habit learning to show how these two parameters interplay to form delusions.  Habit parameter will be on a logarithmic scale
- [ ] Basic figure draft
- [ ] Final figure complete

### Figure 5
Figure 5 will introduce the concept of affect, and show how it can be beneficial, but is in general very inflexible and will lead to poor performance, but not delusions
- [ ] Basic figure draft
- [ ] Final figure complete

### Figure 6
Figure 6 will show a parameter sweep of precision in beliefs about feedback, with a set of different parameters over the habit forming and the strength of the affect
- [ ] Basic figure draft
- [ ] Final figure complete

## Key experimental details
For all experiments we enforce habits to be formed over the advisor, not the cards.  We do this by beginning every experiment with a sequence of 50 trials where the advisor is always trustworthy OR untrustworthy and the cards are entirely random.
