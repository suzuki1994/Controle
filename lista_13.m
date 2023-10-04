close all; clear all; clc;

T = 0.15;
zeta = 0.7;
wn = 2.5;

Gs_desejada = tf(wn^2,[1 2*zeta*wn wn^2]);

Gs = tf(1,[1 1]);
Gz = c2d(Gs,T);
Hs = tf(1,[1 0]);
Hz = c2d(Hs,T);
GHz = c2d(Gs*Hs,T);

s1 = -zeta*wn + 1i*wn*sqrt(1-(zeta^2));
z1 = exp(s1*T);

polos = pole(Gz);
G2 = minreal(tf([1 -polos(1)],1,T)*GHz);

[n_g2,d_g2] = tfdata(G2,'v');
fi2 = angle(polyval(n_g2,z1)/polyval(d_g2,z1));
fi1 = -pi - fi2;
beta = (imag(z1) - real(z1)*tan(-fi1)) /tan(-fi1);
G1 = tf(1,[1 beta],T);

% Condição de módulo
[n_g1g2, d_g1g2] = tfdata(G1*G2,'v');
kc = 1/abs(polyval(n_g1g2,z1)/polyval(d_g1g2,z1));

Cz = kc*tf([1 -polos(1)],1,T)*G1;

FTMA = minreal(Cz*GHz);
FTMF = minreal(Cz*Gz/(1+(Cz*GHz)));

% Condição de módulo e ângulo
[n_FTMA, d_FTMA] = tfdata(FTMA, 'v');
mod = abs(polyval(n_FTMA,z1)/polyval(d_FTMA,z1));
ang = angle(polyval(n_FTMA,z1)/polyval(d_FTMA,z1));
fprintf('Condição de módulo: %.2f\n', mod);
fprintf('Condição de ângulo: %.2f\n\n', ang*180/pi);

% Verificação dos polos de malha fechada
polos_FTMF = pole(FTMF);
z_mod = exp(-T*zeta*wn);
z_ang = T*wn*sqrt(1-(zeta^2));
[x,y] = pol2cart(z_ang, z_mod);
fprintf('Polos desejados: %f %fi\n',x,y);
fprintf('Polos ftma:      %f %fi\n\n', real(polos_FTMF(1)), imag(polos_FTMF(1)));

[n_Gz , d_Gz ] = tfdata(Gz ,'v');
[n_GHz, d_GHz] = tfdata(GHz,'v');
[n_Cz , d_Cz ] = tfdata(Cz ,'v');

Tf = 5;
k = 0:Tf/T;
rk = k*T;

ck(1) = 0; 
vk(1) = 0; 
ek(1) = rk(1) - vk(1);
uk(1) = n_Cz(1)*ek(1);

ck(2) = n_Gz(2)*uk(1) - d_Gz(2)*ck(1); 
vk(2) = n_GHz(2)*uk(1) - d_GHz(2)*vk(1); 
ek(2) = rk(2) - vk(2);
uk(2) = n_Cz(1)*ek(2) + n_Cz(2)*ek(1) - d_Cz(2)*uk(1);
    
for j=3:length(k)
    ck(j) = n_Gz(2)*uk(j-1) - d_Gz(2)*ck(j-1); % Bloco G
    vk(j) = n_GHz(2)*uk(j-1) + n_GHz(3)*uk(j-2) - d_GHz(2)*vk(j-1) - d_GHz(3)*vk(j-2); % Bloco GH
    ek(j) = rk(j) - vk(j); % Erro
    uk(j) = n_Cz(1)*ek(j) + n_Cz(2)*ek(j-1) - d_Cz(2)*uk(j-1); % Bloco C
end

figure(1); hold on;
[y,t] = lsim(FTMF,rk); plot(t,y);
stairs(k*T,ck);
hold off;

figure(2); hold on;
stairs(k*T,ck);
stairs(k*T,vk);
stairs(k*T,ek);
stairs(k*T,uk);
plot(k*T,rk);
hold off; legend('ck','vk','ek','uk');

aux = minreal(tf([1 -1],[1 0],T)*FTMA/T);
[n_erro,d_erro] = tfdata(aux,'v');
Kv = polyval(n_erro,1)/polyval(d_erro,1);
ess = 1/Kv;
fprintf('Erro de regime permanente à rampa: %f\n',ess);





