T = 0.15;
Gs = tf(1, [1, 1]);
Gz = c2d(Gs, T);
Hs = tf([1],[1, 0]);
GHs = Gs*Hs;
GHz = c2d(GHs, T);
Cs = tf([5.175, -4.455],[1, -0.5388], T);
%C:
%e[k] = r[k] - v[k]
%u[k] = 5.175*e[k] - 4.455*e[k-1] + 0.5388*u[k-1]
%G:
%c[k] = 0.1393*u[k-1] + 0.8607c[k-1]
%GH:
%v[k] = 0.01071*u[k1] + 1.861*v[k-1] + 0.01019*u[k-2] - 0.8607*v[k-2]
Tf=5;
kmax = floor(5/T)+1;
t = T * (0:(kmax-1));
r = T * (0:(kmax-1));
u = zeros(1, kmax);
c = zeros(1, kmax);
v = zeros(1, kmax);
e = r - v;
for k = 3:kmax
    v(k) = 0.01071 * u(k-1) + 1.861 * v(k-1) + 0.01019 * u(k-2) - 0.8607 * v(k-2);
    e(k) = r(k) - v(k);
    u(k) = 5.175 * e(k) - 4.455 * e(k-1) + 0.5388 * u(k-1);
    c(k) = 0.1393 * u(k-1) + 0.8607 * c(k-1);
end
plot(t, c);
grid on;
