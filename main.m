clear 
addpath '/Users/benoitmuller/Documents/GitHub/opti-manifolds-rotations-estimation/data generation'
d = 3;
m = 10;
ma = 1;
kappa1 = 5;
kappa2 = 0;
q = 0.8;
problem = build_problem(d, m, ma, kappa1, kappa2, q);
checkgradient(problem);