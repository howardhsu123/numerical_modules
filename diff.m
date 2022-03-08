
% h = 0.1;
% x0 = 2;
% y = [log(x0),log(x0+h),log(x0+2*h)];
% 
% y2 = sec_diff(y,h);
% disp(-1/(x0^2));
% disp(y2);

x = pi/2;
h = 0.001;
f = [sin(x - h),sin(x)];
f1 = first_diff(f,h);
disp(f1);

function f2 = sec_diff(f,h)
    f2 = (f(1) - 2*f(2) + f(3))/(h^2);
end

function f1 = first_diff(f,h)
    f1 = (f(2) - f(1))/h;
end