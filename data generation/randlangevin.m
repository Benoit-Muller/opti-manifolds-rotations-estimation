function Z = randlangevin(n, kappa)
% function Z = randlangevin(n, kappa)
% 
% n is an integer >= 2; kappa is a nonnegative vector of length N.
% Z is an n-by-n-by-N matrix such that each n-by-n slice Z(:,:,i) is a
% random matrix on SO(n) distributed according to the isotropic Langevin
% distribution with concentration kappa(i)*eye(n) around the mean eye(n):
%
% pdf:  (1/c) * exp( kappa(i) * trace( Z(:,:,i) ) )
%       with c as obtained from langevinnormalization(n, kappa(i))
%
% To generate matrices around a mean R in SO(n), simply multiply each slice
% by R (on the left or on the right).
%
% This is a simple acceptance/rejection algorithm with poor performance,
% based on Chikuse's notes in the book Statistics on Special Manifolds
% (where the description is given for O(n) and not for SO(n)). Perhaps a
% better algorithm would be based on the Metropolis-Hastings algorithm,
% like they mention in Chiuso's paper 'Wide sense estimation on the special
% orthogonal group'.

    kappa = kappa(:);
    N = length(kappa);

    Z = zeros(n, n, N);
    
    exact = isinf(kappa);
    Z(:, :, exact) = repmat(eye(n), [1, 1, nnz(exact)]);
    todo = find(~exact);
    
    ntodo = length(todo);
    while ntodo > 0
        Z(:, :, todo) = randrot(n, ntodo);
        intodo = rand(ntodo, 1) >= ...
                        exp(kappa(todo).*(multitrace(Z(:, :, todo))-n));
        todo = todo(intodo);
        ntodo = length(todo);
    end

end