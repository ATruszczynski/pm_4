function [lambdy, f_opt, exitflag, it, LL]=IPM(Q, c, A, b, e, MAX_IT, M)
    [m , n] = size(A);
    
    r = 1;
    dr = 0.2;
    dr_scale = 0.2;
    beta = 0.999;
    
    x0 = startingPoint(Q, c, A, b, r);
    it = 0;
    
    eye_n = ones(n, 1);
    x = x0;
    X = diag(x);
%     y = getLambdas(Q, c, A, b, x, e, r, eye_n);
    y = ones(m, 1);
    y = zeros(m, 1);
    z = r * (X \ eye_n);
    
    while 1
        X = diag(x);
        Z = diag(z);
        
        warning('off','last') % matlab narzeka na odwracanie prawie zerowych macierzy
        
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
            exitflag = -3;
            disp('Timed out')
        end
        
        if exitflag ~= 0
            break
        end
        
        % krok 2
        
        r = dr * (z.' * x) / (n + m);
        dr = dr_scale * dr;

        % krok 3
        
        L = zeros(n+m, n+m);
        L(1:n, 1:n) = -(X \ Z + Q);
        L(1:n, n+1:m+n) = A.';
        L(n+1:n+m, 1:n) = A;
        
        R = zeros(n + m, 1);
        R(1:n, 1) = c - A.' * y - r * (X \ eye_n) + Q * x;
        R(n+1:n+m, 1) = b - A * x;        
        delta = linsolve(L, R);
        
        delta_x = delta(1:n, 1);
        delta_y = delta(n+1:n+m, 1);
        
        delta_z = X \ (r * eye_n - X * Z * eye_n - Z * delta_x);

        % krok 4

        alpha = 1;
        alphas = [alpha];
        
        delta = [delta_x; delta_z];
        xyz = [x; z];

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
        
        f_opt = 1/2 * x.' * Q * x + c.' * x - r * sum(log(x));
    
        it = it + 1;
    end

    lambdy = x;
    LL = [-y; z];
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

function sp = startingPoint(Q, c, A, b, r)
    n = size(A,2);
    
    sp = ones(n, 1);
    
    pos = find(A > 0);
    neg = find(A < 0);
    
    sp(neg) = (length(pos))/(length(neg)) * sp(neg);
end