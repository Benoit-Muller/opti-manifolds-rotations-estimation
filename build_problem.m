function problem = build_problem(d, m, ma, kappa1, kappa2, q, G)
    if ma>m
        error("ma > m impossible.")
    end
    if q<0 || 1<q
        error("q must be in [0,1].")
    end
    problem = build_data(d, m, ma, kappa1, kappa2, q,G);
    M = manyrotationsfactory(d,m,problem.A,problem.Ra);
    problem.M = M;

    problem.cost = @(X) cost(problem,M.add_anchors(X));
    problem.grad = @(X) grad(problem,M.add_anchors(X));
    problem.costgrad = @(X) costgrad(problem,M.add_anchors(X));
    problem.hess = @(X,S) hessian(problem,M.add_anchors(X),M.add_zeros(S));

    problem.MSE = @(X) MSE(problem,X);

    problem.init = @() initialization(problem);
    problem.option.Delta_bar = 2*sqrt(d*(m-ma));
    problem.option.Delta0 = problem.option.Delta_bar / 8;
    problem.option.statsfun = statsfunhelper('mse', @(X) problem.MSE(X));
    problem.option.tolgradnorm = 1e-6 / problem.cardE;
    problem.options.minstepsize = 0;
end


function data = build_data(d, m, ma, kappa1, kappa2, q, G)
    % build the data using a random Erdos-Renyi graph
    % with edge probability ERp = 0.9;
    % anchors are in the first indices
    % output a structure data with fields:
    % 
    % 
    %       d               : int, dimension of rotations
    %       m               : int, number of rotations
    %       cardE           : int, |E| in pdf
    %       I,J             : size (cardE), E={{I(k),J(k)}}_k
    %       H               : size (d,d,cardE)
    %       k1,k2   : floats
    %       q               : float in [0,1]
    %       R               : size(d,d,cardE) true rotations
    %       birth           : (year,month,day,hour,second) moment of birth of the data 
    %       A               : size (length(A)) indices of the anchors
    %       Ra              : size (d,d,|A|) the anchors Ra= R(A)
    %       c1,c2           : float, normalization constants 
    %       deg             : size (cardE), degrees of the graph (1:cardE,I,J)
    %       disconnected    : boolean, true if the graph is connected
    %       maskI,maskJ     : sparse (m,m) adjency matrix of the graph1
    
    % rotations:
    Rtrue = randrot(d, m);
    A = 1:ma;
    Ra = Rtrue(:, :, A);

    % graph:
    switch G
        case "E-R"
            ERp = 0.9;
            [I, J] = erdosrenyi(m, ERp);
            cardE = length(I);
        case "complete"
            I=[];
            J=[];
            for i=1:(m-1)
                for j=(i+1):m
                    I = [I;i];
                    J = [J;j];
                end
            end
            cardE = length(I);
    end

    % data:
    kappa1 = kappa1*ones(cardE, 1);
    kappa2 = kappa2*ones(cardE, 1);
    q = q*ones(cardE, 1);
    Z = randlangevinmixture(d, kappa1, kappa2, q);
    H = multiprod(Rtrue(:, :, I), multiprod(Z, multitransp(Rtrue(:, :, J))));
    
    % build the data:
    prob = buildproblem(d, m, cardE, I, J, H, kappa1, kappa2, q, A, Ra, Rtrue);
    
    % rename fields according to pdf:
    data.d = d;
    data.m = m;
    data.ma = ma;
    data.cardE = prob.M;
    data.I = I;
    data.J = J;
    data.H = prob.H;
    data.k1 = kappa1(1);
    data.k2 = kappa2(1);
    data.q = prob.p(1);
    data.R = Rtrue;
    data.birth = prob.dob;
    data.A = A;
    data.Ra = Ra;
    data.c1 = prob.c1(1);
    data.c2 = prob.c2(1);
    data.deg = prob.d;
    data.disconnected = prob.disconnected;
    data.maskI = prob.maskI;
    data.maskJ = prob.maskJ;
end










