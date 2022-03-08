% spec
%
% v NCQO Newton's Cotes Quadrature Open
% NCQO(f,a,b,n)
%
% v NCQC Newton's Cotes Quadrature Closed
% NCQC(f,a,b,n)
%
% v SQ Simpson's Quadrature
% SQ(f,a,b)
%
% v CNCQO Composite Newton's Cotes Quadrature Open
% CNCQO(f,a,b,n,m) slice to n intervals
% 
% v CNCQC Composite Newton's Cotes Quadrature Closed
% CNCQC(f,a,b,n,m)
%
% v CSQ Composite Simpson's Quadrature
% CSQ(f,a,b,n,m)
%
% v G2Q Gaussian two-points Quadrature
% G2Q(f,a,b)
%
% v G3Q Gaussian three-points Quadrature
% G3Q(f,I)
%
% v CG2Q Composite Gaussian two-points Quadrature
% CG2Q(f,I,n)
%
% v CG3Q Composite Gaussian three-points Quadrature
% CG3Q(f,I,n)
format("long");

f = @(x) sqrt(x^2+1)*sin(x^3-x+2);
a = 10;
b = 20;
exact = CG3Q(f,a,b,10000000);
disp(exact);
m = [1:20];
for n = 1:20
    E(n) = abs(CG2Q(f,a,b,n) - exact);
end

%disp(polyfit(log(m),log(E),1));

function y = NCQO(f,a,b,n)
    if n == 0
        y = (b-a)*f((a+b)/2);
    elseif n == 1
        y = (b-a)/2*(f((2*a+b)/3)+f((a+2*b)/3));
    elseif n == 2
        y = (b-a)/3*(2*f((3*a+b)/4)-f((2*a+2*b)/4)+2*f((a+3*b)/4));
    elseif n == 3
        y = (b-a)/24*(11*f((4*a+b)/5)+f((3*a+2*b)/5)+f((2*a+3*b)/5)+11*f((a+4*b)/5));
    end
end

function y = NCQC(f,a,b,n)
    if n == 1
        y = (b-a)/2*(f(a)+f(b));
    elseif n == 2
        y = (b-a)/6*(f(a)+f((a+b)/2)+f(b));
    elseif n == 3
        y = (b-a)/8*(f(a)+3*f((2*a+b)/3)+3*f((a+2*b)/3)+f(b));
    elseif n == 4
        y = (b-a)/90*(7*f(a)+32*f((3*a+b)/4)+12*f((a+b)/2)+32*f((a+3*b)/4)+7*f(b));
    end
end

function y = SQ(f,a,b,n)
    if n == 1/3
        y = (b-a)/6*(f(a)+4*f((a+b)/2)+f(b));
    elseif n == 3/8
        y = (b-a)/8*(f(a)+3*f((2*a+b)/3)+3*f((a+2*b)/3)+f(b));
    end
end

function y = CNCQO(f,a,b,n,m)
    h = (b-a)/n;
    sum = 0;
    for i = 0:n-1
        sum = sum + NCQO(f,a+i*h,a+(i+1)*h,m);
    end
    y = sum;
end

function y = CNCQC(f,a,b,n,m)
    h = (b-a)/n;
    sum = 0;
    for i = 0:n-1
        sum = sum + NCQC(f,a+i*h,a+(i+1)*h,m);
    end
    y = sum;
end

function y = CSQ(f,a,b,n,m)
    h = (b-a)/n;
    sum = 0;
    for i = 0:n-1
        sum = sum + SQ(f,a+i*h,a+(i+1)*h,m);
    end
    y = sum;
end

function y = G2Q(f,a,b)
    y = (b-a)/2 * (f((b+a)/2 - sqrt(1/3)*(b-a)/2) + f((b+a)/2 + sqrt(1/3)*(b-a)/2));
end

function y = G3Q(f,a,b)
    y = (b-a)/2 * (5/9*f((b-a)/2*(-sqrt(3/5))+(a+b)/2) + 8/9*f((a+b)/2) + 5/9*f((b-a)/2*(sqrt(3/5))+(a+b)/2));
end

function y = CG2Q(f,a,b,n)
    h = (b-a)/n;
    sum = 0;
    for i = 0:n-1
        sum = sum + G2Q(f,a+i*h,a+(i+1)*h);
    end
    y = sum;
end

function y = CG3Q(f,a,b,n)
    h = (b-a)/n;
    sum = 0;
    for i = 0:n-1
        sum = sum + G3Q(f,a+i*h,a+(i+1)*h);
    end
    y = sum;
end

