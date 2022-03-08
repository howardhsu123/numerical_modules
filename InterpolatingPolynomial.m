
%y = LagrangeP([300,400,500,600,700,800,900,1000,1100],[0.024,0.035,0.046,0.058,0.067,0.083,0.097,0.111,0.125],750)
format('long');

n = 10;
for j = 1:n
    xi(j) = cos((2*j-1)*pi/20);
    yi(j) = 1/(1+xi(j)+4*xi(j)*xi(j));
end


% xi2 = [xi(1),xi(2),xi(3),xi(5),xi(7)];
% yi2 = [yi(1),yi(2),yi(3),yi(5),yi(7)];

% max = 0;
% for i = 1:101
%     P = Neville(xi,yi,-1+(i-1)/200);
%     P = P(n,n);
%     f = 1/(1+(-1+(i-1)/200)+4*(-1+(i-1)/200)*(-1+(i-1)/200));
%     if abs(P-f) > max
%         max = abs(P-f);
%     end
% end

Q1 = Neville(xi,yi,1/3);
Q2 = DD(xi,yi);
y = value(Q2,1/3,xi);
disp(Q1(n,n));
disp(y);
% Q3 = DD(xi2,yi2);
% disp(Q1(n,n));
% disp(Q1(4+6-1,4));
% disp(Q2(n,n));
% disp(Q3(5,5));
% disp(max);

for n = 1:1000
    for j = 1:n
        xi(j) = cos((2*j-1)*pi/20);
        yi(j) = 1/(1+xi(j)+4*xi(j)*xi(j));
    end

    max = 0;
    for i = 1:101
        P = Neville(xi,yi,-1+(i-1)/200);
        P = P(n,n);
        f = 1/(1+(-1+(i-1)/200)+4*(-1+(i-1)/200)*(-1+(i-1)/200));
        if abs(P-f) > max
            max = abs(P-f);
        end
    end

    if max > 10e20
        break;
    end
end
%disp(max);

function sum = LagrangeP(xi, fi, x)
    n = length(xi);
    sum = 0;
    for i = 1:n
        toadd = fi(i);
        for j = 1:n
            if (i ~= j)
                toadd = toadd * (x - xi(j)) / (xi(i) - xi(j));
            end
        end
        sum = sum + toadd;
    end
end

function Q = Neville (xi,fi,x)
    m = length(xi);
    Q = zeros(m);
    for j = 1:m
        for i = j:m
            if j == 1
                Q(i,j) = fi(i);
            else
                Q(i,j) = ((x - xi(i-j+1))*Q(i,j-1) -(x - xi(i)) * Q(i-1,j-1))/(xi(i) - xi(i-j+1)) ;
            end           
        end
    end
end

function Q = DD(xi,fi)
    m = length(xi);
    Q = zeros(m);
    for j = 1:m
        for i = j:m
            if j == 1
                Q(i,j) = fi(i);
            else
                Q(i,j) = (Q(i,j-1) - Q(i-1,j-1))/(xi(i) - xi(i-j+1));
            end           
        end
    end
end

function y = value(Q,x,xi)
    m = length(Q);
    sum = 0;
    for i = 1:m
        totime = Q(i,i);
        for j = 1:i-1
            totime = totime .* (x - xi(j));
        end
        sum = sum + totime;
    end
    y = sum;
end