# Stochastic Gradient Descent estimation for gaussian mixtures models with fisher information

Some reproduction and analysis for the algorithm described in this paper : https://arxiv.org/abs/2306.12841. 
It optimizes the step of descent with a finite positive fisher estimation and 
an evolutive learning rate (pre-heating, heating, decreasing step).
We use this algorithm for estimation in gaussian mixtures models and compare it to classic stochastic gradient descent methods.
