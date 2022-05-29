n = 5;

s = zeros(n,1);
for k = 0:n-1
    tmp = n;
    for i = 1:k
        tmp = tmp * (n+i) / i * (n-i) /i;
    end
    if mod(n+k,2)==0
        tmp = -tmp;
    end
    s(k+1) = tmp;
end

b = ones(n,1);
A = hilb(n);
x = zeros(n,1);
r = A*x-b;
p = -r;
step = 1;
Y(1) = vecnorm(x-s);
while vecnorm(r) >= 1e-6
    alpha = - (r.'*p)/(p.'*A*p);
    x = x + alpha*p;
    step = step + 1;
    Y(step) = vecnorm(x-s);
    tmp = r.'*r;
    r = r + alpha*A*p;
    beta = (r.'*r)/tmp;
    p = -r + beta*p;
end
disp(x);

X = linspace(1,step,step);
figure
plot(X,Y(1,1:step));