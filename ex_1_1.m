n = 55;

A = zeros(n, n);
b = zeros(n, 1);
A(1,1) = 6; A(1,2) = 1;
b(1) = 7;
for i = 2:n-1
    A(i,i-1) = 8;
    A(i,i) = 6;
    A(i,i+1) = 1;
    b(i) = 15;
end
A(n,n-1) = 8; A(n,n) = 6;
b(n) = 14;

root = ones(n,1);

x = solveEquationWithLU(A,b);
disp(x.');
disp(vecnorm(x-root));

y = solveEquationWithPLU(A,b);
disp(y.');
disp(vecnorm(y-root));