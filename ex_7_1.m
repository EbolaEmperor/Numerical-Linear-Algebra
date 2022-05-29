n = 100;
A = zeros(n,n);
for i = 1:n-1
    A(i,i) = 4;
    A(i+1,i) = 1;
    A(i,i+1) = 1;
end
A(n,n) = 4;
[D,Q] = symmetricEigen(A)
max(max(abs(A*Q-Q*D)))