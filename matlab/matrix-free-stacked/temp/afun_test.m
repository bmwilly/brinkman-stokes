u = linspace(1, nu + np, nu + np)';

w1 = afun_single(u, p);
w2 = afun_dynamic(u, p);
norm(w1 - w2) / norm(w1)