%--INPUTS------------------------------------------------------------------
%
% U          -     upper triangular matrix
% b          -     column vector
%
%--OUTPUT------------------------------------------------------------------
%
% b          -     solution x of Ux=b, stored in b
%
% Define U and b outside of function.
%format("long");
U=[1.0000   0   0 ;
   -0.6667    1.0000   0;
   0.3333     -0.0833    1.0000];

b=[-3.0000   1.2500   -5.6667 ];
% Call function
b = myBackwardSubstitution(U,b);
function b=myBackwardSubstitution(U,b)
    d=size(U,1);
for i=1:d
    summ = 0;
    for j=1:i-1
        summ = summ + U(i,j)*b(j);
    end
    b(i) = (b(i)-summ)/U(i,i);
end
disp(b);
end
