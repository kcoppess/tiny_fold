# tiny_fold

## parameters.py
### holds parameter values that are used across multiple .py files:
### free energy parameters, priors used in cost function, learning parameters for Adam method, training/testing data for testing code

## partition.py
### holds function definitions for calculating the partition function and gradient (and individual derivatives) for both linear and circular RNA sequences
### recursion modified from NUPACK n^4 algorithm in Dirks and Pierce 2003

## log_partition.py
### holds function definitions for calculating log(partition) and its derivative for both linear and circular RNA sequences

## random_sequence_generator.py
### generates training and testing set of random RNA sequences and calculates their partition functions

## sequences_train.txt & sequences_test.txt
### store the RNA sequences (and partition function) that were used in that run to train/test the model

## stochastic_gradient_descent.py
### implementation of stochastic gradient descent using Adam method

## train_toy.py
### holds the cost function that toy.py minimizes

## toy.py
### trains model using methods in scipy.optimize

## convergence.py
### graphs the value of the cost function after each iteration of SGD to see how it converges
