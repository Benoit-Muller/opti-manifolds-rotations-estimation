%% Question 10,11,12

clear
rng(05011998)
fprintf("\n––– Question 10,11,12 –––\n")
addpath '/Users/benoitmuller/Documents/GitHub/opti-manifolds-rotations-estimation/data generation'
d = 3;
m = 10;
ma = 1;
kappa1 = 5;
kappa2 = 0;
q = 0.8;
G = "E-R";

prob = build_problem(d, m, ma, kappa1, kappa2, q, G);
checkmanifold(prob.M);

%% Question 16

fprintf("\n––– Question 16 –––\n")
checkgradient(prob);
saveas(gcf,'graphics/q16_checkgradient.pdf')
checkhessian(prob);
saveas(gcf,'graphics/q16_checkhessian.pdf')

%% Question 22

fprintf("\n––– Question 22 –––\n")

X0 = initialization(prob);
cost0 = prob.cost(X0);
repeat=100;
cost_random = zeros(repeat,1);
for i=1:repeat
    cost_random(i) = prob.cost(prob.M.rand());
end
[sigma, mu] = std(cost_random);
fprintf("Initial cost : %f\n", cost0)
fprintf("Random cost  : %f +- %f\n", mu, sigma)

%% Question 23

fprintf("\n––– Question 23 –––\n")

X0 = initialization(prob);
mse0 = prob.MSE(X0);
repeat=100;
mse_random = zeros(repeat,1);
for i=1:repeat
    mse_random(i) = prob.MSE(prob.M.rand());
end
[sigma, mu] = std(mse_random);
fprintf("Initial MSE         : %f\n", mse0)
fprintf("Random MSE          : %f +- %f\n", mu, sigma)
fprintf("Expected random MSE : 10.5797 (= 2/3 * pi^2 + 4)  \n")

%% Question 24

fprintf("\n––– Question 24 –––\n")
% a) Generate a synthetic problem
clear
q=1; d=3; m=5; ma=1; kappa1=10; kappa2=10; G="complete";
prob = build_problem(d, m, ma, kappa1, kappa2, q, G);
prob.R(:,:,1) = eye(d);

% b) Choose a random initial starting point
X0 = prob.M.rand();

% c) Run RGD (with backtracking) and RTR
option = prob.option;
option.debug = 0;
option.verbosity =1;

[~, ~, info_rgd, option_rgd] = steepestdescent(prob, X0, option);
[~, ~, info_rtr, option_rtr] = trustregions(prob, X0, option);

% d) Plot the norm of the gradient
figure()
semilogy([info_rgd.iter], [info_rgd.gradnorm],'.-');
p = polyfit([info_rgd.iter], log([info_rgd.gradnorm]),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.gradnorm],'.-');
title("Question 24.d)", "Convergence of the gradient norm for random initialization")
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR")
xlabel("Iteration")
ylabel("Gradient norm")
hold off
saveas(gcf,'graphics/q24d.pdf')

% e) Plot the function value
figure()
plot([info_rgd.iter], -[info_rgd.cost],'.-');
hold on;
plot([info_rtr.iter], -[info_rtr.cost],'.-');
title("Question 24.e)", "Convergence of the loglikelyhood for random initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("Loglikelyhood")
hold off
saveas(gcf,'graphics/q24e.pdf')

% Plot the MSE
figure()
plot([info_rgd.iter], [info_rgd.mse],'.-');
hold on;
plot([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 24.mse", "Convergence of the mse for random initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("MSE")
hold off
saveas(gcf,'graphics/q24mse.pdf')

%% f) Use the spectral initialization
X0 = prob.init();

% g) Run RGD and RTR
option = prob.option;
option.debug = 0;
option.verbosity =1;

[~, ~, info, option_rgd] = steepestdescent(prob, X0, option);
[~, ~, info_rtr, option_rtr] = trustregions(prob, X0, option);

% h) Plot the norm of the gradient and function value

% Gradient
figure()
semilogy([info_rgd.iter], [info_rgd.gradnorm],'.-');
p = polyfit([info_rgd.iter], log([info_rgd.gradnorm]),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.gradnorm],'.-');
title("Question 24.h)", "Convergence of the gradient norm for guessed initialization")
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR")
xlabel("Iteration")
ylabel("Gradient norm")
hold off
saveas(gcf,'graphics/q24h_grad.pdf')

% Function
figure()
plot([info_rgd.iter], -[info_rgd.cost],'.-');
hold on;
plot([info_rtr.iter], -[info_rtr.cost],'.-');
title("Question 24.h)", "Convergence of the loglikelyhood for guessed initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("Loglikelyhood")
hold off
saveas(gcf,'graphics/q24h_fun.pdf')

% MSE
figure()
plot([info_rgd.iter], [info_rgd.mse],'.-');
hold on;
plot([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 24.mse", "Convergence of the mse for guessed initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("MSE")
hold off
saveas(gcf,'graphics/q24h_mse.pdf')


%% Question 25

fprintf("\n––– Question 25 –––\n")
% a) Generate a synthetic problem
clear
q=0.7; d=3; m=40; ma=5; kappa1=5; kappa2=0; G="complete";
prob = build_problem(d, m, ma, kappa1, kappa2, q, G);


% b) Choose a random initial starting point
X0 = prob.M.rand();

% c) Run RGD (with backtracking) and RTR
option = prob.option;
option.debug = 0;
option.verbosity =2;

[~, ~, info_rgd, option_rgd] = steepestdescent(prob, X0, option);
[~, ~, info_rtr, option_rtr] = trustregions(prob, X0, option);

% d) Plot the norm of the gradient
figure()
semilogy([info_rgd.iter], [info_rgd.gradnorm],'.-');
p = polyfit([info_rgd.iter], log([info_rgd.gradnorm]),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.gradnorm],'.-');
title("Question 25.d)", "Convergence of the gradient norm for random initialization")
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR")
xlabel("Iteration")
ylabel("Gradient norm")
hold off
saveas(gcf,'graphics/q25d.pdf')

% e) Plot the function value
figure()
plot([info_rgd.iter], -[info_rgd.cost],'.-');
hold on;
plot([info_rtr.iter], -[info_rtr.cost],'.-');
title("Question 25.e)", "Convergence of the loglikelyhood for random initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("Loglikelyhood")
hold off
saveas(gcf,'graphics/q25e.pdf')

% Plot the MSE
figure()
plot([info_rgd.iter], [info_rgd.mse],'.-');
hold on;
plot([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 25.mse", "Convergence of the mse for random initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("MSE")
hold off
saveas(gcf,'graphics/q25mse.pdf')

%% f) Use the spectral initialization
X0 = prob.init();

% g) Run RGD and RTR
option = prob.option;
option.debug = 0;
option.verbosity =1;

[~, ~, info, option_rgd] = steepestdescent(prob, X0, option);
[~, ~, info_rtr, option_rtr] = trustregions(prob, X0, option);

% h) Plot the norm of the gradient and function value

% Gradient
figure()
semilogy([info_rgd.iter], [info_rgd.gradnorm],'.-');
p = polyfit([info_rgd.iter], log([info_rgd.gradnorm]),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.gradnorm],'.-');
title("Question 25.h)", "Convergence of the gradient norm for guessed initialization")
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR")
xlabel("Iteration")
ylabel("Gradient norm")
hold off
saveas(gcf,'graphics/q25h_grad.pdf')

% Function
figure()
plot([info_rgd.iter], -[info_rgd.cost],'.-');
hold on;
plot([info_rtr.iter], -[info_rtr.cost],'.-');
title("Question 25.h)", "Convergence of the loglikelyhood for guessed initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("Loglikelyhood")
hold off
saveas(gcf,'graphics/q25h_fun.pdf')

% MSE
figure()
plot([info_rgd.iter], [info_rgd.mse],'.-');
hold on;
plot([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 25.mse", "Convergence of the mse for guessed initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("MSE")
hold off
saveas(gcf,'graphics/q25h_mse.pdf')

%% Question 26

fprintf("\n––– Question 26 –––\n")
% a) Generate a synthetic problem
clear
q=0.7; d=3; m=40; ma=5; kappa1=5; kappa2=0; G="E-R";
prob = build_problem(d, m, ma, kappa1, kappa2, q, G);

% f) Use the spectral initialization
X0 = prob.init();

% g) Run RGD and RTR
option = prob.option;
option.debug = 0;
option.verbosity =2;

[~, ~, info_rtr, option_rtr] = trustregions(prob, X0, option);

% h) Plot the norm of the gradient and function value

% Gradient
figure()
semilogy([info_rtr.iter], [info_rtr.gradnorm],'.-');
title("Question 26.h)", "RTR: Convergence of the gradient norm for guessed initialization")
xlabel("Iteration")
ylabel("Gradient norm")
saveas(gcf,'graphics/q26h_grad.pdf')

% Function
figure()
plot([info_rtr.iter], -[info_rtr.cost],'.-');
title("Question 26.h)", "RTR: Convergence of the loglikelyhood for guessed initialization")
xlabel("Iteration")
ylabel("Loglikelyhood")
saveas(gcf,'graphics/q26h_fun.pdf')

% MSE
figure()
plot([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 26.mse", "RTR: Convergence of the mse for guessed initialization")
xlabel("Iteration")
ylabel("MSE")
hold off
saveas(gcf,'graphics/q26h_mse.pdf')


