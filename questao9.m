T= 0.5;
Gs=tf([1,1],[1,4]);
Gz=c2d(Gs, T)

GHs= tf([2,2],[1,4,0]);
GHz = c2d(GHs, T)

Fs = tf([5],[1, 0]);
Fz = c2d(Fs, T)

% recursiva
% G
% c[k] = e[k] - 0.7838*e[k-1] + 0.1353*c[k-1]
% 
% GH
% v[k] = 0.5742*e[k-1] - 0.3581*e[k-2] + 1.135*v[k-1] - 0.1353*v[k-2]
% 
% F
% x[k] =  2.5*r[k-1] + x[k-1]
% 
% e = x - v

kmax=25;
t = T * (0:(kmax-1));
r = ones(1, kmax);

c = zeros(1, kmax);
v = zeros(1, kmax);
x = zeros(1, kmax);
e = zeros(1, kmax);


for k = 3:kmax
    v(k) = 0.5742 * e(k-1) - 0.3581 * e(k-2) + 1.135 * v(k-1) - 0.1353 * v(k-2);
    x(k) = 2.5 * r(k-1) + x(k-1);
    e(k) = x(k) - v(k);
    c(k) = e(k) - 0.7838 * e(k-1) + 0.1353 * c(k-1);
end

% Para mostrar em etapas (steps)
figure(1)
stairs(t, c,'*r');
hold on 
% Para mostrar em tempo contínuo
plot(t, c);
hold off
