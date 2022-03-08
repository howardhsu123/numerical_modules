
function LagrangeP(xi,a,x)
    n = length(a);
    y = a(n) * (x - xi(n-1)) + a(n-1);
    for i = n-2:-1:1
        y = y.*(x - xi(i)) + a(i);
    end
end