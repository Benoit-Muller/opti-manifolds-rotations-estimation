% This code is dedicated to producing a random problem instance. 
% Manopt is needed to run this code.

clear all;
close all;
clc;

% INPUTS:
%
%  n = 2, 3 or 4: work with rotations in SO(n), i.e., n-by-n orthog. matrices
%     * Note: we denoted this by d in the project pdf.
%  N : number of rotations to synchronize
%     * Note: we denoted this by m in the project pdf.
%  M : number of measurements of rotation ratios (aka number of edges in
%      measurement graph)
%  I, J : length-M vectors of indices between 1 and N (encodes the edges of
%         the measurment graph)
%  H : n-by-n-by-M matrix; Each matrix H(:, :, k) in SO(n) is a measurement
%      of the ratio Ra*Rb', where Ra is the I(k)'s rotation matrix to
%      synchronize and Rb is the J(k)'s rotation matrix to synchronize
%  kappa1,2 : length-M vectors of confidences in the M measurements
%  p : probability that a measurement is not an outlier
%     * Note: we denoted this by q in the project pdf.
%  A : vector of indices of the anchors of the problem, i.e.,
%      indices of the known rotation matrices; if omitted, replaced by [1].
%     * Note: we denoted this by J in the project pdf.
%  Ra : nxnx|A| matrix with anchored rotations; if omitted,
%       replaced by the identity matrix eye(n).
% 
% OUTPUTS:
%
%  problem : a structure containing all the given information plus some
%            precomputed data.

% Below is an example with, 100 rotations in SO(3), 
% one anchor, 90% edge density, kappa1=5, kappa2 = 0, p=q=0.8.

fprintf('Generating a problem instance... ');

% Synchronize N rotations in SO(n). Rtrue contains the true rotations we
% are looking for. Of course, these are not known for a real problem
% instance.
n = 3;
N = 100;
Rtrue = randrot(n, N); % ground truth

% Rotations with indices in A are anchored. Ra contains the rotations of
% these anchors. That is: Ra(:, :, k) contains the nxn rotation matrix
% associated to the node A(k), which is anchored. If there are no anchors,
% simply let A = [1] and Ra = eye(n), that is: artificially fix any of the
% rotations to an arbitrary value.
m = 1;
A = 1:m;
Ra = Rtrue(:, :, A);

% For this test problem instance, we generate a random Erdos-Renyi graph.
% For each edge in the graph, we will have a measurement of the relative
% rotation between the adjacent nodes. The data is presented this way:
% I, J are two column vectors of length M, where M is the number of edges.
% There is an edge between nodes I(k) and J(k) for k = 1 : M.
% The graph is symmetric, so that we only define each edge once. That is:
% if the matrix [I J] contains the row [a b], then it does not contain the
% row [b a]. ERp is the edge density in the Erdos-Renyi graph.
% For a complete graph, you may use: [I J] = find(triu(ones(N), 1));
ERp = 0.9;
[I, J] = erdosrenyi(N, ERp);
M = length(I);

% Pick noise parameters and generate noise (random rotation matrices) according
% to these parameters. The measurements are stored in H, a 3D matrix such
% that each slice H(:, :, k) is an nxn rotation matrix which corresponds to
% a measurement of the relative rotation Ri Rj^T,
% with Ri = Rtrue(:, :, I(k)) and Rj = Rtrue(:, :, J(k)).
% The measurement H(:, :, k) is distributed around the real relative
% rotation with parameters kappa1(k), kappa2(k) and p(k), where kappa1,
% kappa2 and p are vectors of length M.
kappa1 = 5.0*ones(M, 1);
kappa2 = 0.0*ones(M, 1);
p = 0.8*ones(M, 1);
Z = randlangevinmixture(n, kappa1, kappa2, p);
H = multiprod(Rtrue(:, :, I), multiprod(Z, multitransp(Rtrue(:, :, J))));

% Put all the data together in a structure which describes the
% synchronization problem. Here, kappa1, kappa2 and p need not be the same
% as those used in generating the measurements.
problem = buildproblem(n, N, M, I, J, H, kappa1, kappa2, p, A, Ra, Rtrue);

fprintf('done.\n');



