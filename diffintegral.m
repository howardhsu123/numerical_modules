% f = @(x) sin(cos(x));
% df = @(a,h)(f(a+h) - f(a-h))/2/h;
% 
% dfext = @(a) cos(cos(a)) .*(-sin(a));
% x = [0:0.1:3];
% err = abs(df(x,0.01) - dfext(x));
% plot(x, df(x,0.001));
% hold on;
% plot(x, dfext(x), '-');   

f = @(x) (sqrt(x.^2 + 5) + exp(1./x)).^sin(sign(2-x).*(abs(2-x)).^(1/5));
df = @(x,h) (-11.*f(x) + 18.*f(x+h) - 9.*f(x+2.*h)+ 2.*f(x+3.*h))./(6.*h);

h = [0.1,0.09,0.08,0.07,0.06,0.05];
d = df(5,h);
exact = df(5,0.00000000000001);
err = abs(d - exact);

disp(err);

% static = (f(10) - f(8))/(10-8);
% a = 8;
% b = 10;
% for i = 1:2000000
%     c = (a + b)/2;
%     tmp = df(c,0.00000000000001);
%     if(tmp > static)
%         b = c;
%     else
%         a = c;
%     end
% end
% 
% disp(c)