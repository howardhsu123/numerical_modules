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
U=[1 2 3;
   0 1 2; 
   0 0 1]
b=[6 3 1];
% Call function
b = myBackwardSubstitution(U,b);
function b=myBackwardSubstitution(U,b)
    d=size(U,1);
for i=1:d-1
    summ = 0;
    for j=d-i+1:d
        summ = summ + U(d-i,j)*b(j);
    end
    b(d-i) = (b(d-i)-summ)/U(d-i,d-i);
end
disp(b);
end
