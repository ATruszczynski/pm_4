function lab8()
    Q = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 0];
    inv(Q)
    c = [0; 0; 0; 0];
    [X, Y] = dane5();
    
%     Y = -Y
    
    size(Y,2)
    asd = find(Y < 0)
    
    X1 = X.';
    X1(:,4) = ones(size(X1,1), 1);
    Y1 = Y.';
    
    
    A = X1 .* Y1;
    A = -A;
    b = -ones(size(A,1),1);
    
    x0 = zeros(4,1);

    e = 1e-10;
    ea = 1e-8;
    
    options = optimoptions('quadprog','Display','iter', 'ConstraintTolerance', e, 'OptimalityTolerance', e);
    
    [a,fval,exitflag,output,lambda] = quadprog(Q,c,A,b,[],[],[],[], [], options);
    
    a
    lambda.ineqlin
    
    m = size(Y,2)
    
    Q = (X.' * X ) .* (Y.' * Y) 
    C = diag(Y)
    Q = C.' * (X.' * X) * C
    det(Q)
    
    c = ones(size(X, 2), 1)
    A = Y
    b = 0
    lb = zeros(size(X, 2), 1)
    
    [x,fval,exitflag,output,llambda] = quadprog(Q, -c, [], [], A, b, lb, [], [], options)
    
    norm(lambda.ineqlin - x)
    llambda.upper
    llambda.ineqlin
    llambda.eqlin
    
%     x0 = -ones(size(A,2),1)
    x0 = ones(size(A, 2), 1)
    [lllambdy, f_opt, exitflagg, it, LL] = IPM(Q, -c, A, b, lb, x0, e, 100, 100000)
    
    lambda.ineqlin
    norm(lllambdy - x)
    lllambdy
    norm(lambda.ineqlin - lllambdy)
    a;
    A = X .* Y;
    aa = A * lllambdy;
    
    inds = find(lllambdy > e);
    ind = inds(1);
    
    
    
    A = X1 .* Y1;
    d = A(ind, :);
    
    L = d(4);
    R = 1 - d(1:3) * aa(1:3);
    
    b = linsolve(L, R);
    
    a
    aa = [aa; b]
    
    norm(aa - a)
end