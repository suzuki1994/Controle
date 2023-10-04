close all; clear; clc;
T = 0.15;

Cz = tf([5.175 -4.455],[1 -0.5388],T);
Gs = tf(1,[1 1]);
Hs = tf(1,[1 0]);

Gz = c2d(Gs, T);
GHz = c2d(Gs*Hs, T);

FTMA = Cz*GHz;
FTMA2=minreal(Cz*GHz);
FTMF = minreal(Cz*Gz/(1+(Cz*GHz)));
FTMF2=feedback(FTMA2,1);
[Cz_n , Cz_d ] = tfdata(Cz , 'v');
[Gz_n , Gz_d ] = tfdata(Gz , 'v');
[GHz_n, GHz_d] = tfdata(GHz, 'v');

Tf = 5;
k = 0:Tf/T;
rk = k*T;

ck(1) = 0;
vk(1) = 0;
ek(1) = rk(1) - vk(1);
uk(1) = Cz_n(1)*ek(1);

ck(2) = Gz_n(2)*uk(1) - Gz_d(2)*ck(1);
vk(2) = GHz_n(2)*uk(1) - GHz_d(2)*vk(1);
ek(2) = rk(2) - vk(2);
uk(2) = Cz_n(1)*ek(2) + Cz_n(2)*ek(1) - Cz_d(2)*uk(1);


for j=3:length(k)
    ck(j) = Gz_n(2)*uk(j-1) - Gz_d(2)*ck(j-1);
    vk(j) = GHz_n(2)*uk(j-1) + GHz_n(3)*uk(j-2) - GHz_d(2)*vk(j-1) - GHz_d(3)*vk(j-2);
    ek(j) = rk(j) - vk(j);
    uk(j) = Cz_n(1)*ek(j) + Cz_n(2)*ek(j-1) - Cz_d(2)*uk(j-1);
end

figure(2);
hold on;
stairs(k*T,ck); % Bloco G
stairs(k*T,vk); % Bloco GH
stairs(k*T,ek); % Erro
stairs(k*T,uk); % Bloco C
plot(k*T,rk);
legend('Bloco G - c(kt)','Bloco GH - v(kt)','Erro - e(kt)','Bloco C - u(kt)');
hold off;

figure(3); hold on;
stairs(k*T,ck);
[y,t] = lsim(FTMF,rk);
plot(t,y); title('Resposta à rampa');
hold off; legend('Saída G(kt)','Saída FTMF');

% Não era necessário calcular
aux = minreal(tf([1 -1],[1 0],T)*FTMA/T);
[n_erro,d_erro] = tfdata(aux,'v');
Kv = polyval(n_erro,1)/polyval(d_erro,1);
ess = 1/Kv;
fprintf('Erro de regime permanente à rampa: %f\n',ess);