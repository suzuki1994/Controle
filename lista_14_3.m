close all; clear all; clc;

% Valores da questão
T = 0.15;
zeta = 0.6;
wn = 5;
Gs = tf(1,[1 5]);
Hs = tf(1,[1 0]);

% Conversão para digital
Gz = c2d(Gs,T);
GHz = c2d(Gs*Hs,T);

% Planta desejada e planta original não controlada
Gs_desejada = tf(wn^2,[1 2*zeta*wn wn^2]); 
FTMF_original = minreal(Gz/(1+GHz));

figure(1); hold on;
step(FTMF_original);
step(Gs_desejada);
hold off;

% Período bom de amostragem
wd = wn*sqrt(1-(zeta^2));
ws = 10*wd;
Ts = 2*pi/ws;

% Parâmetros desejados
Mp = exp(-pi*zeta/sqrt(1-(zeta^2)));
ts5 = 3/(zeta*wn);
fprintf('Settling Time ideal: %f\n',ts5);
fprintf('Overshoot ideal: %f\n\n', Mp*100);

% Polo do controlador
s1 = -zeta*wn + 1i*wn*sqrt(1-(zeta^2));
z1 =  exp(s1*T);

% Construção da parcial G2 da planta controlada
polos = pole(GHz);
G2 = minreal(tf([1 -polos(2)],1,T)*GHz);

% Condição de ângulo
[n2,d2] = tfdata(G2,'v');
fi2 = angle(polyval(n2,z1)/polyval(d2,z1));
fi1 = (-pi - fi2);
beta = (imag(z1) - real(z1)*tan(-fi1)) /tan(-fi1);
G1 = tf(1, [1 beta], T);

% Condição de módulo
[n3,d3] = tfdata(minreal(G1*G2),'v');
K = 1/abs(polyval(n3,z1)/polyval(d3,z1));

% Verificações das especificações
Cz = minreal(tf([1 -polos(2)],1,T)*G1*K);
FTMA = G1*G2*K;
FTMF = minreal(Cz*GHz/(1+(Cz*GHz)));
%FTMF = feedback(Cz*GHz,1);
S_info = stepinfo(FTMF,'SettlingTimeThreshold',0.05);

% Verificação dos parâmetros da planta controlada
fprintf('Settling time: %.5f\n', S_info.SettlingTime);
fprintf('Overshoot: %.2f\n', S_info.Overshoot);
fprintf('Peak: %.4f\n\n', S_info.Peak); 

% Comparação da Planta desejada com a planta controlada
figure(2); hold on;
step(FTMF); step(Gs_desejada);
hold off;

% Verificação das condições de módulo e ângulo
[n4,d4] = tfdata(minreal(FTMA),'v');
c_modulo = abs(polyval(n4,z1)/polyval(d4,z1));
c_angulo = angle(polyval(n4,z1)/polyval(d4,z1))*180/pi;
fprintf('Condição de módulo: %.2f\n',c_modulo);
fprintf('Condição de ângulo: %.2f\n\n',c_angulo);

polos2 = pole(FTMF);
fprintf('Polo de malha fechada do sistema: %f %fi\n',real(polos2(1)), imag(polos2(1)));

mod = exp(-T*zeta*wn);
ang = T*wn*sqrt(1-(zeta^2));
[x,y] = pol2cart(ang,mod);
fprintf('Polo de malha fechada desejado: %f %fi\n',x,y);
