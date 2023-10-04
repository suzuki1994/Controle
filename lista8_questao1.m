T = 0.2;
Gs = tf([1], [1, 2]);
Gz = c2d(Gs, T);
As = tf([1],[1,5]);
Bs = 5;
ABs = parallel(As, Bs);
GHs = Gs * ABs;
GHz = c2d(GHs, T);
Fs = tf([5],[1, 4]);
Fz = c2d(Fs, T);
%Equacoes recursivas 
% G
% c[k] = 0.1648*e[k-1] + 0.6703*c[k-1]
% GH
% v[k] = 0.837*e[k-1] - 0.2952*e[k-2] + 1.038*v[k-1] - 0.2466*v[k-2]
% F
% x[k] =  0.6883*r[k-1] + 0.4493*x[k-1]
% e = x - v
kmax=50;
t = T * (0:(kmax-1));
r = ones(1, kmax);

c = zeros(1, kmax);
v = zeros(1, kmax);
x = zeros(1, kmax);
e = zeros(1, kmax);
for k = 3:kmax
    v(k) = 0.837 * e(k-1) - 0.2952 * e(k-2) + 1.038 * v(k-1) - 0.2466 * v(k-2);
    x(k) = 0.6883 * r(k-1) + 0.4493 * x(k-1);
    e(k) = x(k) - v(k);
    c(k) = 0.1648 * e(k-1) + 0.6703 * c(k-1);
end
figure(1)
stairs(t, c);
FTMF = minreal((Gz/(1+GHz))*Fz);
[xf, yf] = step(FTMF, t);
figure(2)
plot (yf,xf)
