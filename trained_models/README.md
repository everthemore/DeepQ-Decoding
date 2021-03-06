This folder contains all trained models and final results.

In particular, for both bitflip and depolarising noise, the relevant folder contains:

1. A fixed_config.p file, dictionary containing the fixed hyper-parameters used for all simulations.
2. A subfolder for each error rate at which agents were trained


Each error-rate subfolder then contains the results for the best performing agent trained at that error rate. Specifically, each error-rate subfolder contains:

1. all results.p - the average decoded logical qubit lifetime for each error rate at which the agent was tested.
2. final_dqn_weights.h5f - the final weights of the agent. Can be used for reloading this agent.
3. training_history.json - a complete record of the training process through which the agent was obtained. A dictionary containing all relevant metrics.
4. variable_config_xx.p - a dictionary containing the values of the variable hyper-parameters at which this particular agent was obtained.
5. detailed results - a folder containing  all relevant metrics for every episode, at all error rates, at which this agent was evaluated.
