%Matheus Sandim Gonçalves
clear;
clc;
T=0.5;

Gs=tf([1 1],[1 4]);
Hs=tf(2,[1 0]);
Fs=tf(5,[1 0]);

Gz=c2d(Gs,T);
GHz=c2d(Gs*Hs,T);
Fz=c2d(Fs,T);

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

v(1)=0;
x(1)=0;
e(1)=x(1)-v(1);
c(1)=e(1);

v(2) = 0.5742 * e(1) + 1.135 * v(1);
x(2) = 2.5 * r(1) + x(1);
e(2) = x(2) - v(2);
c(2) = e(2) - 0.7838 * e(1) + 0.1353 * c(1);

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

