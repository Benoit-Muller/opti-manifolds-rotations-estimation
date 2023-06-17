addpath '/Users/benoitmuller/Documents/GitHub/opti-manifolds-rotations-estimation/data generation'
clear
rng(05011998)

%% Question 10,11,12

fprintf("\n––– Question 10,11,12 –––\n")
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
fprintf("\nCheck Gradient:\n")
checkgradient(prob);
print("graphics/q16_checkgradient", '-depsc')
fprintf("\nCheck Hessian:\n")
figure()
checkhessian(prob);
print('graphics/q16_checkhessian', '-depsc')

%% Question 22

fprintf("\n––– Question 22 –––\n")

X0 = prob.init();
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

disp("– Random initialization")
disp("RGD:")
[~, ~, info_rgd, ~] = steepestdescent(prob, X0, option);
disp("RTR:")
[~, ~, info_rtr, ~] = trustregions(prob, X0, option);

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
print('graphics/q24d', '-depsc')

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
print('graphics/q24e', '-depsc')

% e) Plot the function value + min reached

step = min([min([info_rgd.cost]), min([info_rtr.cost])]);
figure()
semilogy([info_rgd.iter], [info_rgd.cost]-step,'.-');
p = polyfit([info_rgd.iter], log([info_rgd.cost]-step+ eps()),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.cost]-step,'.-');
semilogy([info_rgd.iter],0*[info_rgd.iter]+ cost(prob,prob.R)-step)
title("Question 24.e)", sprintf("Convergence of the loglikelyhood (plus %f) for random initialization",-step))
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR","Cost of true R")
xlabel("Iteration")
ylabel("Loglikelyhood minus a constant")
hold off
print('graphics/q24e_step', '-depsc')

% Plot the MSE
figure()
semilogy([info_rgd.iter], [info_rgd.mse],'.-');
hold on;
semilogy([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 24.mse", "Convergence of the mse for random initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("MSE")
hold off
print('graphics/q24mse', '-depsc')

%% f) Use the spectral initialization
X0 = prob.init();

% g) Run RGD and RTR
option = prob.option;
option.debug = 0;
option.verbosity =1;

disp("– Spectral initialization")
disp("RGD:")
[~, ~, info_rgd, ~] = steepestdescent(prob, X0, option);
disp("RTR:")
[~, ~, info_rtr, ~] = trustregions(prob, X0, option);

% h) Plot the norm of the gradient and function value

% Gradient
figure()
semilogy([info_rgd.iter], [info_rgd.gradnorm],'.-');
p = polyfit([info_rgd.iter], log([info_rgd.gradnorm]),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.gradnorm],'.-');
title("Question 24.h)", "Convergence of the gradient norm for spectral initialization")
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR")
xlabel("Iteration")
ylabel("Gradient norm")
hold off
print('graphics/q24h_grad', '-depsc')

% Function
step = min([min([info_rgd.cost]), min([info_rtr.cost])]);
figure()
semilogy([info_rgd.iter], [info_rgd.cost]-step,'.-');
p = polyfit([info_rgd.iter], log([info_rgd.cost]-step+ eps()),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.cost]-step,'.-');
semilogy([info_rgd.iter],0*[info_rgd.iter]+ cost(prob,prob.R)-step)
title("Question 24.h)", sprintf("Convergence of the loglikelyhood (plus %f) for spectral initialization",-step))
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR","Cost of true R")
xlabel("Iteration")
ylabel("Loglikelyhood minus a constant")
hold off
print('graphics/q24h_fun', '-depsc')

% MSE
figure()
semilogy([info_rgd.iter], [info_rgd.mse],'.-');
hold on;
semilogy([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 24.h)", "Convergence of the MSE for spectral initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("MSE")
hold off
print('graphics/q24h_mse', '-depsc')


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
disp("– Random initialization")
disp("RGD:")
[~, ~, info_rgd, ~] = steepestdescent(prob, X0, option);
disp("RTR:")
[~, ~, info_rtr, ~] = trustregions(prob, X0, option);

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
print('graphics/q25d', '-depsc')

% e) Plot the function value
step = min([min([info_rgd.cost]), min([info_rtr.cost])]);
figure()
semilogy([info_rgd.iter], [info_rgd.cost]-step,'.-');
p = polyfit([info_rgd.iter], log([info_rgd.cost]-step+ eps()),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.cost]-step,'.-');
semilogy([info_rgd.iter],0*[info_rgd.iter]+ cost(prob,prob.R)-step)
title("Question 25.h)", sprintf("Convergence of the loglikelyhood (plus %f) for random initialization",-step))
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR","Cost of true R")
xlabel("Iteration")
ylabel("Loglikelyhood minus a constant")
hold off
print('graphics/q25e', '-depsc')

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
print('graphics/q25mse', '-depsc')

%% f) Use the spectral initialization
X0 = prob.init();

% g) Run RGD and RTR
option = prob.option;
option.debug = 0;
option.verbosity =2;
disp("– Spectral initialization")
disp("RGD:")
[~, ~, info_rgd, ~] = steepestdescent(prob, X0, option);
disp("RTR:")
[~, ~, info_rtr, ~] = trustregions(prob, X0, option);

% h) Plot the norm of the gradient and function value

% Gradient
figure()
semilogy([info_rgd.iter], [info_rgd.gradnorm],'.-');
p = polyfit([info_rgd.iter], log([info_rgd.gradnorm]),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.gradnorm],'.-');
title("Question 25.h)", "Convergence of the gradient norm for spectral initialization")
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR")
xlabel("Iteration")
ylabel("Gradient norm")
hold off
print('graphics/q25h_grad', '-depsc')

% Function
step = min([min([info_rgd.cost]), min([info_rtr.cost])]);
figure()
semilogy([info_rgd.iter], [info_rgd.cost]-step,'.-');
p = polyfit([info_rgd.iter], log([info_rgd.cost]-step+ eps()),1);
line = exp(polyval(p,[info_rgd.iter]));
hold on;
semilogy([info_rgd.iter], line);
semilogy([info_rtr.iter], [info_rtr.cost]-step,'.-');
semilogy([info_rgd.iter],0*[info_rgd.iter]+ cost(prob,prob.R)-step)
title("Question 26.h)", sprintf("Convergence of the loglikelyhood (plus %f) for spectral initialization",-step))
legend("RGD",sprintf("linear interpolation of slope %f",p(1)),"RTR","Cost of true R")
xlabel("Iteration")
ylabel("Loglikelyhood minus a constant")
hold off
print('graphics/q25h_fun', '-depsc')

% MSE
figure()
plot([info_rgd.iter], [info_rgd.mse],'.-');
hold on;
plot([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 25.mse", "Convergence of the mse for spectral initialization")
legend("RGD","RTR")
xlabel("Iteration")
ylabel("MSE")
hold off
print('graphics/q25h_mse', '-depsc')

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

[~, ~, info_rtr, ~] = trustregions(prob, X0, option);

% h) Plot the norm of the gradient and function value

% Gradient
figure()
semilogy([info_rtr.iter], [info_rtr.gradnorm],'.-');
title("Question 26.h)", "RTR: Convergence of the gradient norm for spectral initialization")
xlabel("Iteration")
ylabel("Gradient norm")
print('graphics/q26h_grad', '-depsc')

% Function
figure()
plot([info_rtr.iter], -[info_rtr.cost],'.-');
hold on;
plot([info_rtr.iter],0*[info_rtr.iter]+ cost(prob,prob.R))
legend("RTR","cost of true R")
title("Question 26.h)", "RTR: Convergence of the loglikelyhood for spectral initialization")
xlabel("Iteration")
ylabel("Loglikelyhood")
print('graphics/q26h_fun', '-depsc')

% MSE
figure()
plot([info_rtr.iter], [info_rtr.mse],'.-');
title("Question 26.mse", "RTR: Convergence of the mse for spectral initialization")
xlabel("Iteration")
ylabel("MSE")
hold off
print('graphics/q26h_mse', '-depsc')

%% Question 27

fprintf("\n––– Question 27 –––\n")
% Generate a synthetic problem
clear
d=3; m=20; ma=5; kappa1=10; kappa2=0; G="E-R";
qq = linspace(0.5,1,5);
repeat = 2;
for j=1:repeat
    fprintf("%i/%i\n",j,repeat)
    for i=1:length(qq)
        prob = build_problem(d, m, ma, kappa1, kappa2, qq(i), G);
        option = prob.option;
        option.verbosity =0;
        [X, ~, info, ~] = trustregions(prob, prob.M.rand(), option);
        if info(end).gradnorm > option.tolgradnorm
            warn("RTR did not reach Gradient norm tolerance")
        end
        rand(i,j) = prob.cost(X);
        [X, ~, info, ~] = trustregions(prob, prob.init(), option);
        if info(end).gradnorm > option.tolgradnorm
            warn("RTR did not reach Gradient norm tolerance")
        end
        init(i,j) = prob.cost(X);
    end
end
diff = rand-init;
%%
alpha=0.05;
[sigma.rand,mu.rand] = std(rand,0,2);
[sigma.init,mu.init] = std(init,0,2);
[sigma.diff,mu.diff] = std(diff,0,2);

figure()

p = plot(qq,mu.rand,'.-'); hold on;
col= p.Color;
interval1 = mu.rand + sigma.rand * norminv(1-alpha/2)/sqrt(repeat) ;
interval2 = mu.rand - sigma.rand * norminv(1-alpha/2)/sqrt(repeat) ;
qq1 = [qq, fliplr(qq)];
inBetween = [interval1', fliplr(interval2')];
fill(qq1, inBetween, col,'FaceAlpha',0.3,"EdgeColor","none");

p=plot(qq,mu.init,'.-');
col= p.Color;
interval1 = mu.init + sigma.init * norminv(1-alpha/2)/sqrt(repeat) ;
interval2 = mu.init - sigma.init * norminv(1-alpha/2)/sqrt(repeat)  ;
inBetween = [interval1', fliplr(interval2')];
fill(qq1, inBetween, col,'FaceAlpha',0.1,"EdgeColor","none");


p=plot(qq,mu.diff,'.-');
col= p.Color;
interval1 = mu.diff + sigma.diff * norminv(1-alpha/2)/sqrt(repeat) ;
interval2 = mu.diff - sigma.diff * norminv(1-alpha/2)/sqrt(repeat)  ;
inBetween = [interval1', fliplr(interval2')];
fill(qq1, inBetween, col,'FaceAlpha',0.1,"EdgeColor","none");

legend("random","","spectral","","difference","")
title("Question 27.a : Cost of RTR according to noise",sprintf("Confidence interval %g",alpha))
xlabel("q")
ylabel("cost")
print('graphics/q27a', '-depsc')


