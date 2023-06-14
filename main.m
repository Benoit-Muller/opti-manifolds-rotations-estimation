%% Question 10,11,12
clear 
addpath '/Users/benoitmuller/Documents/GitHub/opti-manifolds-rotations-estimation/data generation'
d = 3;
m = 10;
ma = 1;
kappa1 = 5;
kappa2 = 0;
q = 0.8;
problem = build_problem(d, m, ma, kappa1, kappa2, q);
%checkmanifold(problem.M);

%% Question 16
%checkgradient(problem);
%saveas(gcf,'graphics/q16_checkgradient.pdf')
checkhessian(problem);