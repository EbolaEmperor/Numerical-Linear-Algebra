[A,D,Q,step] = solve(5)

index = zeros(99,1);
err = zeros(99,1);
steps = zeros(99,1);
for n = 2:100
    index(n-1) = n;
    [A,D,Q,step] = solve(n);
    err(n-1) = norm(A*Q-Q*D,inf);
    steps(n-1) = step;
end
subplot(1,2,1)
plot(index, err);
xlabel('(1) residual of n','fontname', 'Times new roman','interpreter', 'latex');
box off
subplot(1,2,2)
plot(index, steps);
xlabel('(2) iterations of n','fontname', 'Times new roman','interpreter', 'latex');
box off

function [A,D,Q,step] = solve(n)
    A = zeros(n);
    for i = 1:n-1
        A(i,i) = 4;
        A(i,i+1) = 1;
        A(i+1,i) = 1;
    end
    A(n,n) = 4;
    [D,Q,step] = symEigJacobi(A);
end