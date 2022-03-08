A = [-3 2 -1 ;
    6 8 1 ;
    4 2 7 ];

[R, C] = size(A);
L = zeros(R);
U = zeros(R);
for i = 1:R
    L(i, 1) = A(i, 1);
    U(i, i) = 1;
end
for j = 2:R
    U(1, j) = A(1, j) / L(1, 1);
end
        
for i = 2:R
    for j = 2:i
        L(i, j) = A(i, j) - L(i, 1:j - 1) * U(1:j - 1, j);
    end
            
    for j = i + 1:R
        U(i, j) = (A(i, j) - L(i, 1:i - 1) * U(1:i - 1, j)) / L(i, i);
    end
end

disp(L)
disp(U)