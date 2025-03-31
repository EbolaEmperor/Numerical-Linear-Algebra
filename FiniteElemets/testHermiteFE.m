clc
clear
clf

uFun = @(x) sin(pi*x).^2;
duFun = @(x) pi * sin(2*pi*x);
d2uFun = @(x) 2*pi*pi * cos(2*pi*x);
fFun = @(x) -8*pi^4 * cos(2*pi*x);

n = 16;
[u, d, K] = fem4th_Hermite(n, fFun);
x = (0: 1/n : 1)';
err1 = [L2Error(u, d, uFun), H1Error(u, d, duFun), H2Error(u, d, d2uFun)]

n = n * 2;
[u, d] = fem4th_Hermite(n, fFun);
x = (0: 1/n : 1)';

err2 = zeros(1, 3);
set(gcf,'Units','centimeters','Position',[6 6 50 15]);
subplot(1, 3, 1);
title('$u(x)$', 'Interpreter', 'latex');
hold on
plot(x, u, 'o');
err2(1) = L2Error(u, d, uFun);

subplot(1, 3, 2);
title("$u'(x)$", 'Interpreter', 'latex');
hold on
plot(x, d, 'o');
err2(2) = H1Error(u, d, duFun);

subplot(1, 3, 3);
title("$u''(x)$", 'Interpreter', 'latex');
hold on
err2(3) = H2Error(u, d, d2uFun);
err2

order = log2(err1 ./ err2)

function e = L2Error(u, d, uFun)
    n = length(u) - 1;
    p1 = @(x) (1 - x).^2 .* (1 + 2*x);
    p2 = @(x) x.^2 .* (3 - 2*x);
    p3 = @(x,a,b) (x - a).^2 .* (x - b) ./ ((b - a).^2);
    e = 0;
    for i = 1 : n
        a = (i - 1) / n;
        b = i / n;
        uf = @(x) u(i) * p1((x-a)./(b-a)) ...
                + u(i+1) * p2((x-a)./(b-a)) ...
                + d(i) * p3(x, b, a) ...
                + d(i+1) * p3(x, a, b);
        plot(a:(b-a)/20:b, uf(a:(b-a)/20:b), 'LineWidth', 1.5);
        errf = @(x) (uf(x) - uFun(x)).^2;
        e = e + integral(errf, a, b);
    end
    e = sqrt(e);
end

function e = H1Error(u, d, uFun)
    n = length(u) - 1;
    dp1 = @(x) 6 * x .* (x - 1);
    dp2 = @(x) 6 * x .* (1 - x);
    dp3 = @(x,a,b) (x - a) .* (3*x - a - 2*b) ./ ((b - a).^2);
    e = 0;
    for i = 1 : n
        a = (i - 1) / n;
        b = i / n;
        uf = @(x) u(i) * dp1((x-a)./(b-a)) ./ (b-a) ...
                + u(i+1) * dp2((x-a)./(b-a)) ./ (b-a) ...
                + d(i) * dp3(x, b, a) ...
                + d(i+1) * dp3(x, a, b);
        plot(a:(b-a)/20:b, uf(a:(b-a)/20:b), 'LineWidth', 1.5);
        errf = @(x) (uf(x) - uFun(x)).^2;
        e = e + integral(errf, a, b);
    end
    e = sqrt(e);
end

function e = H2Error(u, d, uFun)
    n = length(u) - 1;
    dp1 = @(x) 6 * (2*x - 1);
    dp2 = @(x) 6 * (1 - 2*x);
    dp3 = @(x,a,b) (6*x - 4*a - 2*b) ./ ((b - a).^2);
    e = 0;
    for i = 1 : n
        a = (i - 1) / n;
        b = i / n;
        uf = @(x) u(i) * dp1((x-a)./(b-a)) ./ ((b-a).^2) ...
                + u(i+1) * dp2((x-a)./(b-a)) ./ ((b-a).^2) ...
                + d(i) * dp3(x, b, a) ...
                + d(i+1) * dp3(x, a, b);
        plot(a:(b-a)/20:b, uf(a:(b-a)/20:b), 'LineWidth', 1.5);
        errf = @(x) (uf(x) - uFun(x)).^2;
        e = e + integral(errf, a, b);
    end
    e = sqrt(e);
end