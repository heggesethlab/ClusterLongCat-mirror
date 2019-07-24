
## Concept of a Mixture of Markov Models
This mixture model assumes that there are $K$ distinct clusters in the data, each of which can be modeled by a Markov process. This assumes that the probability of an observation moving from one state to the next state depends only on the current state. That is that $P(S_{t+1}|S_{t},S_{t-1},\dots,S_0) = P(S_{t+1}|S_{t})$ where $S_t$ is the state at position $t$ in the sequence of that observation.

This model naturally fits longitudinal categorical data, where each state in a category is a state in the Markov chain and each time is a position. The model is easily represented either through a set of graphs of the chains displaying transition probabilities, or matrices of transition probabilities, where one graph or matrix corresponds to one cluster. However, this model doesn't fit all situations. Often times the Markov assumption is too strict, and previous states influence the probability of transitioning to the next state.

We provide three methods of modeling data this way: A maximum likelihood estimation approach, a Bayesian approach, and an extension relying on Bayesian estimation allowing for the Markov chains within each cluster to vary. All of these models are computationally complex and grow as the number of states, number of observations, and length of observations grows.

## Using Maximum Likelihood Estimation

### Mixture of Markov Models

Using maximum likelihood estimation (MLE) calculated with the expectation maximization (EM) algorithm, the parameters for the models and the cluster assignments of observations are computed. (Helske and Helske 2019) This method works for sequences of variable length, but computational complexity grows quickly for longer sequences and for more states. In order to prevent getting stuck at local minimums on the likelihood surface, we run the EM algorithm with 10 restarts for each model.

To support interpreting the results, we create several visualizations for MLE mixtures of Markov models

#### Plots for interpretation
* Posterior probability of observation belonging to most probable cluster
* Estimated Markov chains
* Estimated transition matrices probabilities


## Using Bayesian Estimation


In light of the rising field of Bayesian statistics, Bayesian clustering methods were discovered to yield similar results as frequentist approaches when the sample size is large. Some advantages of Bayesian analysis include a natural way of combing prior knowledge with data and interpretable inferences without relying on asymptotic approximation.

### Mixture of Markov Models

Similar to the Mixture of Markov Models using MLE, the Bayesian inference method also assumes that each member of a cluster can be modeled by a cluster specific transition matrix. It is assumed that the rows of this matrix are independent are that each row is modeled by a Dirichlet distribution, and that this underlying distribution is identical for all clusters. Notably, this method's complexity grows much more slowly for longer length sequences, and is recommended when performing clusterings on the uncompressed data.  Our app supports setting the 2 values for priors. The parameters of the Dirichlet distributions corresponding to on diagonal entries in the transition matrices and those corresponding to off diagonal entries. Parameters of this distribution are estimated using Markov chain Monte Carlo (MCMC). (Pamminger and Frühwirth-Schnatter 2010)

To support interpreting the results, we create several visualizations for Bayesian estimated mixture of Markov models:

#### Plots for Model Analysis
* Posterior probability of each observation belonging to its most probable cluster (calculated over the last 50 MCMC draws)
* Estimated Markov chains
* Average of estimated transition probabilities across MCMC runs after burn-in
* Standard deviations of transition probabilities across MCMC runs after burn-in

#### Trace plots for MCMC evaluation
For trace plots, if there is a visible trend or if the chain is getting stuck at one value for a significant number of iterations, the MCMC isn't estimating the values properly. Either it should be run for longer to reach a stable position, or the parameters should be changed to encourage it to explore the probability space more. We create:

* Trace plot of estimated cluster sizes ($eta_h$)
* Trace plot of estimated diagonal entries of transition matrices


### Dirichlet Multinomial Clustering

As opposed to Mixture of Markov Models clustering, where each member of a cluster is assumed to be modeled by a cluster specific transition matrix, Dirichlet Multinomial clustering method assumes that each member of a cluster has their own transition matrix. Instead of all clusters having one Dirichlet distribution as priors, each cluster has its own Dirichlet distribution, whose parameters are modeled as a negative multinomial distribution with hyperparameters that can be set in the app interface. This allows each observation to be modeled by its own transition matrix, while the cluster specific transition matrix is modeled by the expectation of the Dirichlet distributions for that cluster. Parameters of the Distribution are estimated using MCMC. (Pamminger and Frühwirth-Schnatter 2010)

To support interpreting the results, we create several visualizations for Dirichlet Multinomial clustering:

#### Plots for model analysis
* Posterior probability of each observation belonging to its most probable cluster (calculated over the last 50 MCMC draws)
* Estimated Markov chains from the expectation of cluster specific Dirichlet distributions averaged across MCMC iterations after burn-in (*equation 3 from paper*)
* Average of Estimated Transition probabilities from the expectation of cluster specific Dirichlet distributions averaged across MCMC iterations after burn-in (*equation 3 from paper*)
* Standard Deviations of Transition Probabilities across MCMC runs after burn-in
* Posterior expectation of within cluster variation for each transition ($e_{h,ij}$ *from paper*)
* Posterior expectation and standard deviation of row specific unobserved heterogeneity in each cluster (*equation 5 from paper*)

#### Trace plots for MCMC evaluation
For trace plots, if there is a visible trend or if the chain is getting stuck at one value for a significant number of iterations, the MCMC isn't estimating the values properly. Either it should be run for longer to reach a stable position, or the parameters should be changed to encourage it to explore the probability space more. We create:

* Trace plot of estimated cluster sizes ($eta_h$ *from paper*)
* Trace plot of estimated diagonal entries of transition matrices
* Trace plot of diagonal entries of estimated within group heterogeneity ($e_{h,ii}$ *from paper*)
