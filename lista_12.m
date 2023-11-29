close all; clear all; clc;
T = 0.1;

G1s = tf(1,[1 0]);
G2s = tf(1,[1 3]);
Hs = tf(1,[1 6]);

G1z = c2d(G1s,T);
G2z = c2d(G2s,T);
Hz = c2d(Hs,T);

Gz = minreal(G1z+G2z);
FTMA = Gz*Hz;
FTMF = Gz/(1+(Gz*Hz));

Tf = 35;
k = 0:Tf/T;
rk = ones(length(k));

ak(1) = 0;
bk(1) = 0;
ck(1) = ak(1) + bk(1);
vk(1) = 0;
ek(1) = rk(1) - vk(1);

for j=2:length(k)
    ak(j) = 0.1*ek(j-1) + 1*ak(j-1);
    bk(j) = 0.08639*ek(j-1) + 0.7408*bk(j-1);
    ck(j) = ak(j) + bk(j);
    vk(j) = 0.0752*ck(j-1) + 0.5488*vk(j-1);
    ek(j) = rk(j) - vk(j);
end

figure(1); hold on;
step(FTMF,'b');
stairs(k*T,ck,'r');
hold off; legend('FT de malha fechada','Equação Recursiva C(kT)');

figure(2); hold on;
stairs(k*T,ak);
stairs(k*T,bk);
stairs(k*T,ck);
stairs(k*T,vk);
stairs(k*T,ek);
hold off;
legend('ak','bk','ck','vk','ek');

[n,d] = tfdata(minreal(tf([1 -1],1,T)*FTMA/T),'v');
Kp = polyval(n,1)/polyval(d,1);
ess = 1/Kp;
