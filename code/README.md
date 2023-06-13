# Code

This files contains all the code for the Dolphins repository

## PART 1: Data Wrangling

In this [script](https://github.com/bankheak/Dolphins/blob/main/code/social_associations.R) I organize 30 years of data into lists and correct sample sizes for each step of analysis.

## PART 2: Social Associations

In this [script](https://github.com/bankheak/Dolphins/blob/main/code/social_associations.R) I calculate dyadic social associations using the [simple-ratio index](https://github.com/bankheak/Dolphins/blob/main/code/functions.R) to create an association matrix for each year of data.

## PART 3: Network Structure

In this [script](https://github.com/bankheak/Dolphins/blob/main/code/network_structure.R) I create undirected weighted social networks using the association matrix and caluculate local and global network metrics.

## PART 4: Regression Models

In this [script](https://github.com/bankheak/Dolphins/blob/main/code/mantel_tests.R) I run multiple linear regression models using dyadic HI simmularity, individual-level variables, homerange overlap and genetic relatedness as predictors for the assciation index.

## PART 5: Social Diffusion Analysis

In this [script]() I fit an OADA to a network that uses different combinations of vertical learning, horizontal learning, genetic relatedness and environmental simularity networks.
