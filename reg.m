
%y = [0.45,0.7,0.95,1.2,1.85,1.45,1.55,1.7,1.85,2.15,2.2,2.4,2.65];
%x = [12,25,33,43,50,52,59,67,73,82,87,90,95];

% Delta = [0,4,6,8,10,12,14,16,18,20,22,24];
% P = [0,1.25,1.85,2.4,3.05,3.64,4.25,4.85,5.45,6.05,6.7,7.25];
% A = (0.507 / 2) ^ 2 * pi;
% L = 2;
% epsilon = P/A;
% sigma = Delta / L;

t = [0,1,2,3,4,5,6,7,8,9,10];
n = [5.03,4.71,4.4,3.97,3.88,3.62,3.3,3.15,3.08,2.92,2.7];

[a,b] = reg(t,1./n);
disp(a);
disp(b);
plot(epsilon,sigma);

function [a,b] = reg(x,y)
    n = length(y);
    sumxy = 0;
    sumx = 0;
    sumy = 0;
    sumx2 = 0;
    for i = 1:n
        sumxy = sumxy + x(i)*y(i);
        sumx = sumx + x(i);
        sumy = sumy + y(i);
        sumx2 = sumx2 + x(i)^2;
    end
    b = (n*sumxy - sumx * sumy)/(n*sumx2 - sumx^2);
    a = 1/n * (sumy - sumx * b);
end