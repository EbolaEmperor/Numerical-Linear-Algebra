B = xlsread('ch3data.xlsx');
A = [ones(28,1) B(1:28,1:11)];
b = B(1:28, 12);
x = solveLSwithQR(A,b);
disp(x.');
disp(vecnorm(A*x-b));