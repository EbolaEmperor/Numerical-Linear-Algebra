A = zeros(100,100);
A(1,1) = 10;
for i = 2:100
    A(i,i) = 10;
    A(i-1,i) = 1;
    A(i,i-1) = 1;
end
x = ones(100,1);
b = A*x;
x1 = solveEquationWithCholesky(A,b);
x2 = solveEquationWithCholeskyImproved(A,b);
x3 = solveEquationWithLU(A,b);
x4 = solveEquationWithPLU(A,b);
x5 = solveLSwithQR(A,b);
disp([vecnorm(x1-x),vecnorm(x2-x),vecnorm(x3-x),vecnorm(x4-x),vecnorm(x5-x)]);
