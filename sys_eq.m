% Systemetic RK Runge-Kutta Method
% Systemetic two-step Adam's Bashforth
% Systemetic two-step Adam's Moulton
% Systemetic three-step Adam's Bashforth
% Systemetic three-step Adam's Moulton
% Systemetic four-step Adam's Bashforth
% Systemetic four-step Adam's Moulton



epsilon = 4;
%f = @(t,v) [-v(1)*(lambda - v(2) + v(1)), - v(2)* (v(2) - v(1))];
%f = @(t,v) [1+v(1)/t]
%f = @(t,v) [v(2), epsilon * (1-v(1)^2)* v(2) - v(1)];
%f = @(t,v) [v(2),v(3),-1/2*v(1)*v(3)];
%f = @(t,v) [v(2), -v(1) - 2*v(2)];
%f = @(t,v) [v(2),v(3),-4*v(3)-5*v(2)];
%f = @(t,v) [v(2), v(3), (sin(v(1)^2)+cos(t)^3-v(3))/(1+t^2)];

% Optimal RK2
% A = [0,0;
%     2/3,0];
% B = [1/4,3/4];
% C = [0, 2/3]

% A = [0,0,0,0;
%     1/3,0,0,0;
%     -1/3,1,0,0;
%     1,-1,1,0];
% B = [1/8,3/8,3/8,1/8];
% C = [0,1/3,2/3,1];

g = 9.8;
L = 1;

a = 4;
b = 5;

for k = 1:100
    lambda = (a+b)/2;
    f = @(t,v) [v(2), (-lambda+sqrt(2+cos(t)))*v(1)];

    N = 30001;

    % RK4
    A = [0,0,0,0;
    1/2,0,0,0;
    0,1/2,0,0;
    0,0,1,0];
    B = [1/6,1/3,1/3,1/6];
    C = [0,1/2,1/2,1];

    ini_t = 0;
    alpha = [1,0];
    h = 1/N;

    [w,t] = RK(f,A,B,C,N,h,alpha,ini_t);
    disp(w(N+1,1));
    if(w(N+1,1) > 0)
        a = (a+b)/2;
    else
        b = (a+b)/2;
    end
end
disp(lambda);

function [w,t] = RK(f,A,B,C,N,h,alpha,ini_t)
    s = length(B);
    n = length(alpha);

    for i = 1:N+1
        t(i) = ini_t + (i - 1)*h;
    end

    for p = 1:n
        w(1,p) = alpha(p);
    end

    for i = 1:N
        for p = 1:n
            tmp(p) = 0;
        end
        for j = 1:s
            for p = 1:n
                sum(p) = 0;
            end
            for l = 1:j-1
                for p = 1:n
                    sum(p) = sum(p) + A(j,l) .* k(l,p);
                end
            end
            tmp2 = f(t(i) + C(j) * h, w(i,:) + sum);
            for p = 1:n
                k(j,p) = h * tmp2(p);
                tmp(p) = tmp(p) + k(j,p) .* B(j);
            end
        end
        for p = 1:n
            w(i+1,p) = w(i,p) + tmp(p);
        end
    end
end

function [w,t] = AB2(f,N,h,alpha,ini_t)
    n = length(alpha);
    for i = 1:N+1
        t(i) = ini_t + (i - 1) * h;
    end

    % The initial w is determined by RK
    A = [0,0;
        2/3,0];
    B = [1/4,3/4];
    C = [0, 2/3];
    tmp_w = RK(f,A,B,C,2,h,alpha,ini_t);

    for p = 1:n
        w(1,p) = tmp_w(1,p);
        w(2,p) = tmp_w(2,p);
    end

    for i = 2:N
        f_i = f(t(i),w(i,:));
        f_i1 = f(t(i-1),w(i-1,:));
        for p = 1:n
            w(i+1,p) = w(i,p) + h * (3/2 * f_i(p) - 1/2 * f_i1(p));
        end
    end
end

function [w,t] = AB3(f,N,h,alpha,ini_t)
    n = length(alpha);
    for i = 1:N+1
        t(i) = ini_t + (i - 1) * h;
    end

    for p = 1:n
        w(1,p) = alpha(p);
    end

    for i = 1:N
        f_i = f(t(i),w(i,:));
        f_i1 = f(t(i-1),w(i-1,:));
        f_i2 = f(t(i-2),w(i-2,:));
        for p = 1:n
            w(i+1,p) = w(i,p) + h * (23/12 * f_i(p) - 4/3 * f_i1(p) + 5/12 * f_i2(p));
        end
    end
end

function [w,t] = AB4(f,N,h,alpha,ini_t)
    n = length(alpha);
    for i = 1:N+1
        t(i) = ini_t + (i - 1) * h;
    end

    for p = 1:n
        w(1,p) = alpha(p);
    end

    for i = 1:N
        f_i = f(t(i),w(i,:));
        f_i1 = f(t(i-1),w(i-1,:));
        f_i2 = f(t(i-2),w(i-2,:));
        f_i3 = f(t(i-3),w(i-3,:));
        for p = 1:n
            w(i+1,p) = w(i,p) + h * (55/24 * f_i(p) - 59/24 * f_i1(p) + 37/24 * f_i2(p) - 9/24 * f_i3(p));
        end
    end
end