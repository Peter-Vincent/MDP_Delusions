# MDP_Delusions

This repository contains data and code relating to the MDP_Delusions project, a collaboration between Rick Adams, Thomas Parr, Peter Vincent, David Benrimoh and Karl Friston

## Project summary

Delusions are, by popular definition, false beliefs that are held with certainty and resistant to contradictory evidence. They seem at odds with the notion that the brain at least approximates Bayesian inference. This is especially the case in schizophrenia, a disorder thought to relate to decreased – rather than increased – certainty in the brain’s model of the world. We use an active inference Markov decision process model (a Bayes-optimal decision-making agent) to perform a simple task involving social and non-social inferences. We show that even moderate changes in some model parameters – decreasing confidence in sensory input and increasing confidence in states implied by its own (especially habitual) actions – can lead to delusions as defined above. Incorporating affect in the model increases delusions, specifically in the social domain. The model also reproduces some classic psychological effects, including choice-induced preference change, and an optimism bias in inferences about oneself. A key observation is that no change in a single parameter is both necessary and sufficient for delusions; rather, delusions arise due to conditional dependencies that create ‘basins of attraction’ which trap Bayesian beliefs. Simulating the effects of antidopaminergic antipsychotics – by reducing the model’s confidence in its actions – demonstrates that the model can escape from these attractors, through this synthetic pharmacotherapy.

## Using the code

All the code (almost) needed to run the simulations included in this work can be found in this repo.  You will also need a version of Matlab (this project was developed using MATLAB R2019b) and SPM12, which you can download [here](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

Figures **2a**, **2b**, **3a**, **3b**, **3c**, **4a**, **4b** and **4c** can be run quickly by running `Production_Figures/single_MDPs_for_figures.m` which will load in the relevant simulations used for each graphic - note that paths at the beginning of the file (lines 8, 9, 11, 12) will need to be updated for your setup.  Alternatively, this script also gives the code and parameters required to re-create each simulation.  These can be easily changed within the script.

To re-create the remaining figures (which show many simulations to demonstrate the effect changing parameters on the results), use `working_model/analyse_varied_MDPs.m`.  Unfortunately, we are not able to upload all of the simulation results here, since some of the simulations are >5GB.  We are very happy to share these large files, so please email [Dr. Rick Adams](rick.adams@ucl.ac.uk) and we'll find a way to share these files.  However, you can still generate some the figures already using the results in `Production_Figures/MDP_sim_symm_beta*` and `Production_Figures/MDP_sim_trust_beta*`.  Again, you will need to change the paths at the beginning of the code.  

You can also email Rick if you have any further questions about the code, or feel free to raise an issue here.

The online manuscript can be found [here](https://www.sciencedirect.com/science/article/pii/S0920996421003054?via%3Dihub)











