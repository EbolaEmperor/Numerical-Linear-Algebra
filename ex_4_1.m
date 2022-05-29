n = 100;
h = 1.0/n;
epsilon = 0.0001;
a = 0.5;
A = zeros(n+1,n+1);
b = a*h*h*ones(n+1,1);
b(1) = 0;
b(n+1) = 1;
s = zeros(n+1,1);
for i = 0:n
    s(i+1) = (1-a)/(1-exp(-1/epsilon))*(1-exp(-i/(n*epsilon)))+a*i/n;
end

A(1,1) = 1;
for i = 2:n
    A(i,i-1) = epsilon;
    A(i,i) = -(2*epsilon+h);
    A(i,i+1) = epsilon+h;
end
A(n+1,n+1) = 1;

eps = 1e-7;
[x1,step1] = jacobi(A,b,eps);
disp([step1 vecnorm(x1-s)]);
[x2,step2] = gauss_seidel(A,b,eps);
disp([step2 vecnorm(x2-s)]);
[x3,step3] = sor(A,b,1.01,eps);
disp([step3 vecnorm(x3-s)]);

% m = 110;
% X = zeros(m,1);
% Y = zeros(m,1);
% for i = 1:m-1
%     X(i) = 1+0.0001*i;
%     [~, Y(i)] = sor(A,b,X(i),eps);
% end
% X(m) = X(m-1)+0.00008;
% [~, Y(m)] = sor(A,b,X(m),eps);
% plot(X,Y)
