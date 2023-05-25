function problem = build_problem(d, m, ma, kappa1, kappa2, q)
    problem = build_data(d, m, ma, kappa1, kappa2, q);
    problem.M = manyrotationsfactory(d,m,problem.A,problem.Ra);
    problem.cost = @(X) cost(problem,M.add_anchors(X));
    problem.grad = @(X) grad(problem,M.add_anchors(X));
    problem.costgrad = @(X) costgrad(problem,M.add_anchors(X));
end


function data = build_data(d, m, ma, kappa1, kappa2, q)
    % build the data using a random Erdos-Renyi graph
    % with edge probability ERp = 0.9;
    % output a structure data with fields:
    % 
    % 
    %       d               : int
    %       m               : int
    %       M               : int, |E| in pdf
    %       I,J             : size (M), E={{I(k),J(k)}}_k
    %       H               : size (d,d,M)
    %       kappa1,kappa2   : floats
    %       q               : float in [0,1]
    %       R               : size(d,d,M) true rotations
    %       birth           : (year,month,day,hour,second) moment of birth of the data 
    %       A               : size (M) indices of the anchors
    %       Ra              : size (d,d,|A|) the anchors Ra= R(A)
    %       c1,c2           : float, normalization constants 
    %       deg             : size (M), degrees of the graph (1:M,I,J)
    %       disconnected    : boolean, true if the graph is connected
    %       maskI,maskJ     : sparse (m,m) adjency matrix of the graph1
    
    % rotations:
    Rtrue = randrot(d, m);
    A = 1:ma;
    Ra = Rtrue(:, :, A);

    % graph:
    ERp = 0.9;
    [I, J] = erdosrenyi(m, ERp);
    M = length(I);

    % data:
    kappa1 = kappa1*ones(M, 1);
    kappa2 = kappa2*ones(M, 1);
    q = q*ones(M, 1);
    Z = randlangevinmixture(d, kappa1, kappa2, q);
    H = multiprod(Rtrue(:, :, I), multiprod(Z, multitransp(Rtrue(:, :, J))));
    
    % build the data:
    prob = buildproblem(d, m, M, I, J, H, kappa1, kappa2, q, A, Ra, Rtrue);
    
    % rename fields according to pdf:
    data.d = d;
    data.m = m;
    data.M = prob.M;
    data.I = I;
    data.J = J;
    data.H = prob.H;
    data.kappa1 = kappa1;
    data.kappa2 = kappa2;
    data.q = prob.p;
    data.R = Rtrue;
    data.birth = prob.dob;
    data.A = A;
    data.c1 = prob.c1;
    data.c2 = prob.c2;
    data.deg = prob.d;
    data.disconnected = prob.disconnected;
    data.maskI = prob.maskI;
    data.maskJ = prob.maskJ;
end










