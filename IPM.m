function [lambdy, f_opt, exitflag, it, LL]=IPM(Q, c, A, b, lb, x0, e, MAX_IT, M)
    [m , n] = size(A);
    
    r = 10;
    dr = 0.2;
    beta = 0.999;
    
    eye_n = ones(n, 1);
    eye_m = ones(m, 1);
    x = x0;
    X = diag(x);
    y = getLambdas(Q, c, A, b, x, e, r, eye_n);
    y = ones(m, 1);
%     y = zeros(m, 1);
    z = r * (X \ eye_n);
    y = ones(m, 1);
    z = ones(n, 1);
    it = 0;
    
    x = -x;
%     y = -y
    z = -z;
    
    while 1
        X = diag(x);
        X_inv = inv(X);
        Z = diag(z);
        Y = diag(y);
%         Y_inv = inv(Y);
        
        % krok 1
        
        
        
        exitflag = checkEndConds(Q, c, A, b, x, y, z, e, M);
        
        if exitflag == 1
            disp('Ok')
        elseif exitflag == -1
            disp('Primal unbound')
        elseif exitflag == -2
            disp('Dual unbound')
        end
        
        if it > MAX_IT
            exitflag = -3
            disp('Timed out')
        end
        
        if exitflag ~= 0
            break
        end
        
        % krok 2
%         r = r * dr;
        r = dr * (z.' * x) / (n + m);

        % krok 3

%         L = zeros(2*n + m);
%         L(1:m, 1:n) = A;
%         L(m+1:m+n, 1:n) = -Q;
%         L(n+m+1:2*n+m, 1:n) = Z;
%         L(m+1:m+n, n+1:n+m) = A.';
%         L(m+1:m+n, n+m+1:2*n+m) = eye(n);
%         L(n+m+1:2*n+m, n+m+1:2*n+m) = X;
% 
%         R = zeros(2*n+m, 1);
%         R(1:m, 1) = b - A*x;
%         R(m+1:m+n, 1) = c - A.' *  y - z + Q*x;
%         R(m+n+1:m+2*n, 1) = r*eye_n - X * Z * eye_n;
% 
%         delta = linsolve(L, R);
% 
%         delta_x = delta(1:n, 1);
%         delta_y = delta(n+1:n+m, 1);
%         delta_z = delta(n+m+1:2*n+m, 1);
        
        %%%%%%%%%%%%%%%%%%%
        L = zeros(n+m, n+m);
        L(1:n, 1:n) = -(X \ Z + Q);
        L(1:n, n+1:m+n) = A.';
        L(n+1:n+m, 1:n) = A;
        
        R = zeros(n + m, 1);
        R(1:n, 1) = c - A.' * y - r * (X \ eye_n) + Q * x;
        R(n+1:n+m, 1) = b - A * x + r * (Y \ eye_m);
        
        delta = linsolve(L, R);
%         ddelta = inv(L) * R
        
        delta_x = delta(1:n, 1);
        delta_y = delta(n+1:n+m, 1);
        
        delta_z = X \ (r * eye_n - X * Z * eye_n - Z * delta_x);
        delta = [delta; delta_z];
        
        ori = 1

        % krok 4

        alpha = 1;

        xyz = [x; y; z];
        alphas = [alpha];

        for i = 1:size(xyz, 1)
            if delta(i, 1) < 0
                alphas = [alphas, -beta * xyz(i, 1) / delta(i, 1)];
            end
        end

        step = min(alphas);

        % krok 5

        x = x + step * delta_x;
        y = y + step * delta_y;
        z = z + step * delta_z;
        
        f_opt = 1/2 * x.' * Q * x + c.' * x;
        
        ori = 1
    
        it = it + 1
    end

    lambdy = x
    f_opt = 0
    it = 0
    LL = 0
end

function llambda = getLambdas(Q, c, A, b, x, e, r, eye_n)
   X = diag(x);
   d = A * x - b;
   aa = x;
   llambda = [];
   ind = [];

    for i = 1:size(A,1)
        if d(i, :) > -e * 1e4
            ind = [ind, i];
        end
    end
    
    L_A = A(ind,:).';
    L_b = Q * x + c - r * (X \ eye_n);
    tmp_llambda = linsolve(L_A, L_b);

    llambda = zeros(size(A,1), 1);
    llambda(ind, :) = tmp_llambda;
end

function flag = checkEndConds(Q, c, A, b, x, y, z, e, M)
    ro = b - A * x;
    sigma = c - A.' * y - z + Q * x;
    gamma = z.' * x;
    
    max_x = max(x);
    max_y = max(y);
    
    flag = 0;
    
    if max_x > M
        flag = -1;
    end
    
    if max_y > M
        flag = -2;
    end
    
    ro_norm = sum(abs(ro));
    sigma_norm = sum(abs(sigma));
    
    if ro_norm < e && sigma_norm < e && gamma < e
        flag = 1;
    end
    
end

function sp = startingPoint()
    
end