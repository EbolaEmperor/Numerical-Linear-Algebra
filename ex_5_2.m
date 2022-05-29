N = 100;
index = zeros(N,1);
err = zeros(N,1);
step = zeros(N,1);
for n = 1:N
    A = hilb(n);
    b = zeros(n,1);
    for i = 1:n
        b(i) = sum(A(i,1:n))/3;
    end
    index(i) = i;
    [x,step(i)] = CG(A,b,1e-14);
    err(i) = vecnorm(x-ones(n,1)/3);
end
subplot(1,2,1)
plot(index, err);
xlabel('(1) residual of n','fontname', 'Times new roman','interpreter', 'latex');
box off
subplot(1,2,2)
plot(index, step);
xlabel('(2) iterations of n','fontname', 'Times new roman','interpreter', 'latex');
box off
