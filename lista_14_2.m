close all; clear all; clc;
T = 0.2;

Gs = tf(1,[1 5]);
Gz = c2d(Gs,T);

Cz = tf(0.1813, [1 -0.8187],T);
FTMA = Cz*Gz;
FTMF = feedback(Cz*Gz,1);
step(FTMF);

[n,d] = tfdata(Gz,'v');

Tf = 5;
k = 0:Tf/T;
rk = ones(length(k));

uk(1) = 0;
ck(1) = 0;
ek(1) = rk(1) - ck(1);

for j=2:length(k)
    uk(j) = 0.1813*ek(j-1) + 0.8187*uk(j-1);
    ck(j) = n(2)*uk(j-1) - d(2)*ck(j-1);
    ek(j) = rk(j) - ck(j);
end

figure(2); hold on;
stairs(k*T,uk);
stairs(k*T,ck);
stairs(k*T,ek);
hold off; legend('uk','ck','ek');

% Considerando degrau unitário
[n2,d2] = tfdata(minreal(FTMA),'v'); 
Kp = polyval(n2,1)/polyval(d2,1);
ess = 1/(1+Kp);

fprintf('Kp: %.2f\n', Kp);
fprintf('ess: %.4f\n', ess);
