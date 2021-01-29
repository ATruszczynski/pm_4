function lab8()
    H = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 0];
    cc = [0; 0; 0; 0];
    [X, Y] = dane1();
    
    X1 = X.';
    X1(:,4) = ones(size(X1,1), 1);
    Y1 = Y.';
    
    A = X1 .* Y1;
    A = -A;
    b = -ones(size(A,1),1);

    e = 1e-10;
    
    options = optimoptions('quadprog','Display','iter', 'ConstraintTolerance', e, 'OptimalityTolerance', e);
    
    [a,fval0,exitflag,output,lambda] = quadprog(H,cc,A,b,[],[],[],[], [], options);
    
    C = diag(Y);
    Q = C.' * (X.' * X) * C;
    
    c = ones(size(X, 2), 1);
    A = Y;
    b = 0;
    lb = zeros(size(X, 2), 1);
    
    [x,fval,exitflag,output,llambda] = quadprog(Q, -c, [], [], A, b, lb, [], [], options);
    
    [lllambdy, f_opt, exitflagg, it, LL] = IPM(Q, -c, A, b, e, 200, 100000)
    
    if exitflagg == 1
        quadLambda = [llambda.eqlin; llambda.lower]
        norm(LL - quadLambda)
        norm(fval - f_opt)
        norm(lllambdy - x)
        norm(1/2 * a'*H*a + cc' * a + (1/2 * x' * Q * x - c' * x))
        
        A = X .* Y;
        aa = A * lllambdy;
        
        inds = find(lllambdy > e);
        ind = inds(1);
        
        A = X1 .* Y1;
        d = A(ind, :);
    
        L = d(4);
        R = 1 - d(1:3) * aa(1:3);
        b = linsolve(L, R);
        aa = [aa; b];
        
        norm(aa - a)
        aa - b
        
        wykres(X, Y, aa)
        set(gcf,'color','w');
        
        
    end
    it
    exitflag
    exitflagg
    ori = 0
end