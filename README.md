# MDP_Delusions

This repository contains data and code relating to the MDP_Delusions project, a collaboration between Rick Adams, Thomas Parr, Peter Vincent, David Benrimoh and Karl Friston

## Project summary

Delusions are, by popular definition, false beliefs that are held with certainty and resistant to contradictory evidence. They seem at odds with the notion that the brain at least approximates Bayesian inference. This is especially the case in schizophrenia, a disorder thought to relate to decreased – rather than increased – certainty in the brain’s model of the world. We use an active inference Markov decision process model (a Bayes-optimal decision-making agent) to perform a simple task involving social and non-social inferences. We show that even moderate changes in some model parameters – decreasing confidence in sensory input and increasing confidence in states implied by its own (especially habitual) actions – can lead to delusions as defined above. Incorporating affect in the model increases delusions, specifically in the social domain. The model also reproduces some classic psychological effects, including choice-induced preference change, and an optimism bias in inferences about oneself. A key observation is that no change in a single parameter is both necessary and sufficient for delusions; rather, delusions arise due to conditional dependencies that create ‘basins of attraction’ which trap Bayesian beliefs. Simulating the effects of antidopaminergic antipsychotics – by reducing the model’s confidence in its actions – demonstrates that the model can escape from these attractors, through this synthetic pharmacotherapy.



