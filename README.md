# tiny_fold

### partition.py
    holds function definitions for calculating the partition function and its derivative for both linear and circular RNA sequences
    recursion modified from NUPACK n^4 algorithm in Dirks and Pierce 2003

### training_sequence_generator.py
    generates training set of random RNA sequences and calculates their partition functions

### test_sequence_generator.py
    generates testing set of random RNA sequences and their partition functions to test how trained model predicts partition function for given sequence

### sequences_train.txt & sequences_test.txt
    store the RNA sequences (and partition function) that were used in that run to train/test the model

### train_toy.py
    trains the model on the training data generated by minimizing a modified least squares

### toy.py
    tests the trained model using the test data generated

### toy_analytical.py
    has analytical derivative of partition function (just for linear sequences)
