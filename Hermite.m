format("long");
%xi = [0,0.75,1.5,2.0856,2.676193,3.219694,3.748513,4.279179,4.821254,5.000];
%fi = [0,445.903683,842.695315,1095.211197,1295.955674,1437.602773,1542.363644,1621.280769,1680.890649,1696.803710];
%fi1 = [600,573.579644,477.074216,384.947629,296.576145,226.796410,171.475176,127.808738,93.728061,84.473801];

xi = [1,2,3];
fi = xi .* exp(-xi);
fi1 = [0,-exp(-2),-2*exp(-3)];

%Q = H_Interpolate(xi,fi,fi1);
x = 1.5;
%x = [1:0.01:3];
y = value(Q,x,xi);
y2 = H_cubic(xi,x,fi);
%plot(x,y2);
%hold on;
%plot(x,y);

%disp(Q);
disp(y);
disp(y2);

function Q = H_Interpolate(xi,fi,fi1)
    m = length(xi);
    m = 2*m;
    Q = zeros(m);
    for j = 1:m
        for i = j:m
            if j == 1
                Q(i,j) = fi(floor((i-1)/2)+1);
            else
                if (j == 2 && mod(i,2) == 0)
                    Q(i,j) = fi1(i/2);
                else
                    Q(i,j) = (Q(i,j-1) - Q(i-1,j-1))/(xi(floor((i - 1)/2)+1) - xi(floor((i - j)/2)+1));
                end
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
            totime = totime .* (x - xi(floor((j-1)/2)+1));
        end
        sum = sum + totime;
    end
    y = sum;
end

function y = H_cubic(xi, x, fi)
    n = length(xi);
    y = 0;
    for j = 1:n-1
        for i = x
            if (i > xi(j) && i < xi(j+1))
                y = y + (1 - 2.*(x - xi(j))/(xi(j) - xi(j+1))) .* (((x - xi(j+1))/(xi(j) - xi(j+1))).^2) .* fi(j);
                y = y + (1 - 2.*(x - xi(j+1))/(xi(j+1) - xi(j))) .* (((x - xi(j+1))/(xi(j+1) - xi(j))).^2) .* fi(j+1);
                y = y + (x - xi(j)) .* (((x - xi(j+1))/(xi(j) - xi(j+1))).^2) .* fi(j);
                y = y + (x - xi(j+1)) .* (((x - xi(j))/(xi(j+1) - xi(j))).^2) .* fi(j+1);
            end
        end
    end
end