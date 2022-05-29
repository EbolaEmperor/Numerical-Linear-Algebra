n = 40;
A = hilbert(n);
b = zeros(n,1);
s = ones(n,1);
for i = 1:n
    b(i) = sum(A(i,1:n));
end
x1 = solveEquationWithCholesky(A,b);
x2 = solveEquationWithCholeskyImproved(A,b);
x3 = solveEquationWithLU(A,b);
x4 = solveEquationWithPLU(A,b);
disp([vecnorm(x1-s),vecnorm(x2-s),vecnorm(x3-s),vecnorm(x4-s)]);