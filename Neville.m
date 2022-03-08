xi = [-1,0,1,2];
yi = [5,1,1,11];

Q = Nevilles(xi,yi,1.5);
disp(Q)

function Q = Nevilles(xi, fi, x)
Q(:,1) = fi;
n = length(xi);
a = x - xi;
for j = 2:n
    for i = j:n
        Q(i,j) = a(i-j+1) * Q(i,j-1) + a(i) * Q(i-1,j-1) /(xi(i) - xi(i-j+1));
    end
end
end
