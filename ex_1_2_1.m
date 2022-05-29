A = zeros(100,100);
A(1,1) = 10;
for i = 2:100
    A(i,i) = 10;
    A(i-1,i) = 1;
    A(i,i-1) = 1;
end
b = unifrnd(0,20,100,1);
x1 = solveEquationWithCholesky(A,b);
x2 = solveEquationWithCholeskyImproved(A,b);
x3 = solveEquationWithLU(A,b);
x4 = solveEquationWithPLU(A,b);
disp([vecnorm(A*x1-b),vecnorm(A*x2-b),vecnorm(A*x3-b),vecnorm(A*x4-b)]);
