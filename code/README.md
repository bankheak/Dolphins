# Code

This file contains all the code for the Dolphins repository

## Data analysis process
<img src="https://github.com/user-attachments/assets/910b3c00-0a76-4735-b39b-6cb383488591" align="middle" width="500px"/>

## PART 1: Data Wrangling

I start by fixing and combining data from 1993-2014 and seperating 6-year periods between 1995-2012. I then only include individuals that have been seen at least 10 times in all three study periods.

## PART 2: Social Association Matrix

I calculated the dyadic social associations in each study period using the Simple Ratio-Index: $SRI_{AB} = x/(x + y_{AB} + y_A + y_B)$

Where:
- $x$ = [the number of sampling periods with A and B observed associated]
- $y_{AB}$ = [the number of sampling periods with just A identified]
- $y_A$ = [the number of sampling periods with just B identified]
- $y_B$ = [the number of sampling periods with A and B identified but not associated]

## PART 3: CV and Modularity

I create permutated null models to calculate random distributions of the coefficient of variation ($CV$) and modularity ($Q$). For the $CV$ I sent the cv_null.R script to OSU's high performance cluster to run 1000 permutation of each study period's social network to get the distribution of random $CV$'s compared to each observed $CV$. For $Q$ I calculated the cluster edge betweeness of the 1000 permutations in each study period and compared them to each observed $Q$.

## PART 4: Create Individual-level Variables (ILV) and Human-centric Foraging Behavior Predictors

### Sex and Age Matrices

I fixed sex and age data of the 117 dolphins and created similarity matrices of both sex and age. For sex, I assigned pairs that had the same sex a 1 and those of opposite sex a 0. For age I calculated the absolute difference in the birth year of each pair and subtracted that scaled value by 1 to get a similarity in age. I then calculated some demographics of age and sex.

### Home-range Overlap Matrix

I calculated the home-range overlap of each pair of dolphins using the kernel density overlap: $𝐻𝑅𝑂_{𝑖𝑗} = (𝑅_{𝑖𝑗}/𝑅_𝑖)*(𝑅_{𝑖𝑗}/𝑅_𝑗)$

Where:
- $𝑅_{𝑖𝑗}$ = [the kernel density overlap of individual $i$ and $j$]
- $𝑅_𝑖$ = [the kernel density estimate of individual $i$]
- $𝑅_𝑗$ = [the kernel density estimate of individual $j$]

### Genetic relatedness Matrix

I used the paternity data of each individual to calculate the pedidree genetic relateness of each pair of individuals.

### Human-centric Foraging Behavior Similarity Matrix

I calculated the proportion of time each indidvidual spent engaging in human-centric foraging behaviors compared to the total number of observations they were seen. I then calculated the difference in these proportions between individuals and subtracted their scaled values by one to get a similarity in human-centric foraging behavior engagement.

## PART 5: Run MCMC GLMM

I built a multi-membership generalized linear mixed model (GLMM) using a Markov chain Monte Carlo (MCMC) sampler under a Bayesian statistical framework. I compared leave-one-out cross-validation information criteria (LOOIC) between a null model, an additive model, and an interactive model. I selected the model with the lowest significant change in the expected log-likelihood predictive density (ELPD) of the approximated LOOIC. The most parsimonious model was the interaction model, which described each covariates’ effect on the dyadic association index:

Process model:
$SRI_{i,j,p} ~ (u_i+u_j)+β1_{Sex_{i,j}}+β2_{Age_{i,j}}+β3_{HRO_{i,j,p}}+β4_{HC_{i,j,p}}+β5_{During HAB_{i,j,p}}+β6_{After HAB_{i,j,p}}+β7_{HC_{i,j,p}}*During HAB+β8_{HC_{i,j,p}}*After HAB$

Observation model:
$SRI_{i,j,p} ~ Beta[μ_{SRI_{i,j,p}},ϕ_{SRI_{i,j,p}}]

I then create the figures of the effect sizes of these predictors.

## PART 6: Display Networks

I create the figure for the network containing clusters on and off the Sarasota map.

