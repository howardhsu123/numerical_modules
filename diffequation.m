%
%
% E Euler Method
% 
% RK Runge-Kutta Method

% f = @(t,x) (-t) .* tan(x)/(1+t.^2) ;
% x = @(t) asin(sqrt(1./(2+2.*t.^2)));
% N = 4;
% alpha = pi/4;
% a = 0;
% b = 1;
% 
% [w,t] = E(f,N,alpha,a,b);
% 
% disp(t);
% disp(w);
% disp(abs(w - x(t)));
% plot(t,w);
% hold on
% plot(t,x(t));

% f = @(t,x) (4.*x+t.^5.*exp(t))./t;
% alpha = 0;
% a = 1;
% b = 6;
% N = 5./(1/64);
% A = [0,0,0,0;
%     1/3, 0,0,0;
%     -1/3,1,0,0;
%     1,-1,1,0];
% B = [1/8,3/8,3/8,1/8];
% C = [0,1/3,2/3,1];
% 
% x = @(t) t.^4.*(exp(t) - exp(1));
% [w,t] = RK(f,A,B,C,N,alpha,a,b);
% 
% err = abs(x(t)-w);
% disp(err);

f = @(t,x) (3.*t - x./t);
alpha = 1;
a = 1;
b = 6;
h = [1/10, 1/100, 1/1000, 1/10000, 1/100000, 1/1000000];
N = (b-a)./h;
A = [0,0;
    1/2,0];
B = [0,1];
C = [0,1/2];
x = @(t) t.^2 ;

for k = 1:6
    [w,t] = RK(f,A,B,C,N(k),alpha,a,b);

    err = abs(x(t)-w);

    %disp(err);

    for i = 1:N+1
        max(k) = -1;
        if(err(i) > max(k))
            max(k) = err(i);
            index = i;
        end
    end

    disp(max(k));

end

disp(polyfit(log(h),log(max),1));



function [w,t] = E(f,N,alpha,a,b)
    h = (b-a)./N;
    w(1) = alpha;
    for i = 1:N+1
        t(i) = a + (i - 1).*h;
    end
    for i = 1:N
        w(i+1) = w(i)+h.*f(t(i),w(i));
    end
end

function [w,t] = RK(f,A,B,C,N,alpha,a,b)
    s = length(B);
    h = (b-a)./N;
    for i = 1:N+1
        t(i) = a + (i - 1).*h;
    end
    w(1) = alpha;
    for i = 1:N
        tmp = 0;
        for j = 1:s
            sum = 0;
            for l = 1:j-1
                sum = sum + A(j,l) .* k(l);
            end
            k(j) = h .* f(t(i) + C(j).*h, w(i) + sum);
            tmp = tmp + k(j) * B(j);
        end
        w(i+1) = w(i) + tmp;
    end
end