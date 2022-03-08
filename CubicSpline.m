format("long");
%xi = [0.00,0.25,0.50,0.75,1.00];
%fi = [0.00,0.176777,0.5,0.530330,0.000];

%xi = [300, 400, 500, 600, 700, 800, 900, 1000, 1100];
%fi = [0.024, 0.035, 0.046, 0.058, 0.067, 0.083, 0.097, 0.111, 0.125];

%xi = [-1.0,-0.5,0.0,0.5,1.0];
%fi = [0,0.82436,1.00,0.9098,0.73576];
%dfa = 2.71828;
%dfb = -0.36788;

%xi = [0,10,20,30,40,50,60,70,80,90,100];
%fi = [1.788,1.307,1.003,0.799,0.657,0.548,0.467,0.405,0.355,0.316,0.283]
%fi = [0.611,1.227,2.337,4.242,7.375,12.34,19.92,31.16,47.35,70.11,101.3]
%fi = [1402,1447,1482,1509,1529,1542,1551,1553,1554,1550,1543]
%xi = [0,500,1000,1500,2000,2500,3000]
%fi = [1.2255,1.1677,1.1120,1.0583,1.0067,0.9570,0.9092] 
%[a1,b1,c1,d1] = nature(xi,fi);
%[a2,b2,c2,d2] = NotAKnot(xi,fi);
%disp(Ans(xi,a1,b1,c1,d1,800));
%disp(Ans(xi,a2,b2,c2,d2,800));
%disp(Ans(xi,a1,b1,c1,d1,1600));
%disp(Ans(xi,a2,b2,c2,d2,1600));
%disp(Ans(xi,a1,b1,c1,d1,2350));
%disp(Ans(xi,a2,b2,c2,d2,2350));
%disp(Ans(xi,a1,b1,c1,d1,2790));
%disp(Ans(xi,a2,b2,c2,d2,2790));


xi = [1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50];
fi = [72.1, 65.0, 48.7, 43.6, 38.9 35.4, 28.7];

[a,b,c,d] = nature(xi,fi);

disp(Int(xi,a,b,c,d));


% max = 0;
% for x = 0:0.001:19
%     key = abs(Ans(xi,a,b,c,d,x)-6+2^cos(0.4*x));
%     if(key > max)
%         max = key;
%     end
% end
% disp(max);

%x = [0:0.1:100];
%f = (x+1) .* exp(-x);
%hold on;
%plot(x,f);

% for i = 1:10
%     x = [xi(i):0.01:xi(i+1)];
%     y = a(i) + b(i) * (x-xi(i)) + c(i) * (x-xi(i)).^2 + d(i) * (x-xi(i)).^3;
%     hold on;
%     plot(x, y);
%     hold on;
%     plot(xi(i),fi(i),'r*');
% end


function [a,b,c,d] = Clamped(xi, fi, dfa, dfb)
    n = length(xi);
    a = fi;
    hi = diff(xi);
    di = 2*(hi(1:n-2)+hi(2:n-1));
    di = [2*hi(1), di, 2*hi(n-1)];
    A = diag(di) + diag(hi,-1) + diag(hi,1);
    %disp(A);
    da = diff(a);
    f(1) = da(1)*3/hi(1) - 3*dfa;
    f(2:n-1) = 3*diff(da./hi);
    f(n) = -da(n-1)*3/hi(n-1)+3*dfb;
    f = f';
    %disp(f);
    c = A\f;
    for i = 1:n-1
        b(i) = (a(i+1) - a(i))/hi(i) - (2*c(i) + c(i+1)) / 3*hi(i);
        d(i) = (c(i+1) - c(i))/3/hi(i);
    end
    a = a(1:n-1);
    c = c(1:n-1);
    a = a';
    b = b';
    d = d';
end

function [a,b,c,d] = NotAKnot(xi, fi)
    n = length(xi);
    a = fi;
    hi = diff(xi);
    di = 2*(hi(1:n-2)+hi(2:n-1));
    di = [hi(2), di, hi(n-2)];
    A = diag(di) + diag(hi,-1) + diag(hi,1);
    A(1,2) = -(hi(1) + hi(2));
    A(1,3) = hi(1);
    A(n,n-2) = hi(n-1);
    A(n,n-1) = -(hi(n-1) + hi(n-2));
    %disp(A);
    da = diff(a);
    f(1) = 0;
    f(2:n-1) = 3*diff(da./hi);
    f(n) = 0;
    f = f';
    %disp(f);
    c = A\f;
    for i = 1:n-1
        b(i) = (a(i+1) - a(i))/hi(i) - (2*c(i) + c(i+1)) / 3*hi(i);
        d(i) = (c(i+1) - c(i))/3/hi(i);
    end
    a = a(1:n-1);
    c = c(1:n-1);
    a = a';
    b = b';
    d = d';
end

function [a,b,c,d] = nature(xi, fi)
    n = length(xi);
    a = fi;
    hi = diff(xi);
    di = 2*(hi(1:n-2)+hi(2:n-1));
    di = [1, di, 1];
    A = diag(di) + diag(hi,-1) + diag(hi,1);
    A(1,2) = 0;
    A(n,n-1) = 0;
    %disp(A);
    da = diff(a);
    f(1) = 0;
    f(2:n-1) = 3*diff(da./hi);
    f(n) = 0;
    f = f';
    %disp(f);
    c = A\f;
    for i = 1:n-1
        b(i) = (a(i+1) - a(i))/hi(i) - (2*c(i) + c(i+1)) / 3*hi(i);
        d(i) = (c(i+1) - c(i))/3/hi(i);
    end
    a = a(1:n-1);
    c = c(1:n-1);
    a = a';
    b = b';
    d = d';
end

function t = Ans(xi,a,b,c,d,x)
    n = length(xi);
    for i = 1:n-1
        if (x >= xi(i) && x <= xi(i+1))
            t = a(i) + b(i) * (x-xi(i)) + c(i) * (x-xi(i)).^2 + d(i) * (x-xi(i)).^3;
        end
    end
end

function t = Int(xi,a,b,c,d)
    n = length(xi);
    t = 0;
    for i = 1:n-1
        t = t + a(i)*(xi(i+1) - xi(i)) + 1/2*b(i) * (xi(i+1)-xi(i))^2 + 1/3*c(i) * (xi(i+1)-xi(i)).^3 + 1/4* d(i) * (xi(i+1)-xi(i)).^4;
    end
end