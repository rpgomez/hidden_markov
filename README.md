## Introduction
This package implements order 1 Hidden Markov Models
for discrete categorical variables. It also provides a jupyter
notebook documenting how the Hidden Markov Model's parameters are
computed and updated.

## Installation:

1. Clone this package.
2. pip install .

## Usage:

1. You have a sequence of observed symbols y_n.
2. Instantiate a hmm object to fit your data:

```python
import hmm

# instantiate a model with 5 hidden states
my_model = hmm.hmm(num_hiddenstates=5)

# fit the model to the data
my_model.fit(y_n)

# Compute the log likelihood of the observed data:
my_model.compute_loglikelihood()

# Compute the gamma sequence for the data:
gamma = my_model.transform()

# Compute the most probable hidden state at time t:
most_likely = gamma[t].argmax()

# Return the transition matrix:
my_model.A

# Recover the observable probability matrix:
my_model.B

# Recover the most likely initial hidden state:
my_model.pi.argmax()

# Map symbol 'foo' to observable integer 'k':
my_model.map['foo']

# Convert observable integer 'k' to the corresponding symbol 'foo':
my_model.inv_map[k]
```
