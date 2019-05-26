###########################
# Stan Model specification for
#
# Fitting dependent and independent causal models for the language dataset from
# 'Human sound systems are shaped by post-Neolithic changes in bite configuration'
#
# Tarasov, Uyeda 2019
##########################

#--- Interpretation of the model parameters and variables, see also Supplementary text

# N = N_i number of languages per Area i
# N_hg = H_i number of languages with the hunter-gathering society per Area i
# L_hg = HL_i number of languages where the society is hunter-gathering and labiodental phonemes are present per Area i
# L_agr = AL_i number of languages where the society is agricultural and labiodental phonemes are present per Area i

# theta_S = theta_Hi probability of languages with hunter-gathering societies per Area i
# theta_L1 = theta_L1i probability of languages with labiodental phonemes present for hunter-gatherers and per Area i
# theta_L2 = theta_L2i probability of languages with labiodental phonemes present for agricultural societies and per Area i
#--------------

# Dependent model for the parameter inference

M_G1_dep =
  "
data {
int N; // number of languages per Area
int N_hg; // number of languages where the society is HG per Area
int L_hg; // number of languages where the society is HG and LB sounds are present per Area
int L_agr; // number of languages where the society is AGR and LB sounds are present per Area
}

parameters {
real<lower=0,upper=1> theta_S; 
real<lower=0,upper=1> theta_L1;
real<lower=0,upper=1> theta_L2;
}

model {
theta_S ~ uniform(0, 1);
theta_L1 ~ uniform(0, 1);
theta_L2 ~ uniform(0, 1);
N_hg ~ binomial(N, theta_S);
L_hg ~ binomial(N_hg, theta_L1);
L_agr ~ binomial(N - N_hg, theta_L2);
}
"

# Log scale dependent model for the marginal Ln estimation

M_G1_dep_Log =
  "
data {
int N;
int N_hg;
int L_hg;
int L_agr;
}

parameters {
real<lower=0,upper=1> theta_S;
real<lower=0,upper=1> theta_L1;
real<lower=0,upper=1> theta_L2;
}

model {
target += uniform_lpdf(theta_S | 0, 1);
target += uniform_lpdf(theta_L1 | 0, 1);
target += uniform_lpdf(theta_L2 | 0, 1);
target += binomial_lpmf(N_hg | N, theta_S);
target += binomial_lpmf(L_hg | N_hg, theta_L1);
target += binomial_lpmf(L_agr | N - N_hg, theta_L2);
}
"


###############################

# Independent model for the parameter inference

M_G1_ind =
  "
data {
int N;
int N_hg;
int L_hg;
int L_agr;
}

parameters {
real<lower=0,upper=1> theta_S;
real<lower=0,upper=1> theta_L1;
}

model {
theta_S ~ uniform(0, 1);
theta_L1 ~ uniform(0, 1);

N_hg ~ binomial(N, theta_S);
L_hg ~ binomial(N_hg, theta_L1);
L_agr ~ binomial(N - N_hg, theta_L1);
}
"





# Log scale independent model for the marginal Ln estimation

M_G1_ind_Log =
  "
data {
int N;
int N_hg;
int L_hg;
int L_agr;
}

parameters {
real<lower=0,upper=1> theta_S;
real<lower=0,upper=1> theta_L1;
}

model {
target += uniform_lpdf(theta_S | 0, 1);
target += uniform_lpdf(theta_L1 | 0, 1);

target += binomial_lpmf(N_hg | N, theta_S);
target += binomial_lpmf(L_hg | N_hg, theta_L1);
target += binomial_lpmf(L_agr | N - N_hg, theta_L1);
}
"