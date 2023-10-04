%projeto 1- controlador digital
%Matheus Sandim Gonçalves
clear;
clc;
format long
tp=24.8e-3;%instante de pico
delta1=112e-3;%regime permanente
delta2=504e-3;%sobressinal
Mp= delta1/delta2; %sobressinal
zeta = sqrt((log(Mp)^2) / ((pi^2) + log(Mp)^2));%fator de amortecimento
Wn = (pi) / (tp * sqrt(1 - zeta^2));% frequência natural
num = Wn^2;%numerador da função de tranferência da planta analógica
den = [1, 2*zeta*Wn, Wn^2];%denominador da função de tranfêrencia da planta analógica
ts5 = 3/(zeta * Wn); %tempo de acomodação
Gs=tf(num,den);%função de tranferência da planta analógica
t = linspace(0, 0.15, 1000); % Vetor de tempo
[yf, xf] = step(Gs, t);%gráfico da função de tranferência
figure(1)
plot(xf, yf);
title('Resposta do sistema a um degrau unitário');
grid on;
Mp2 = (max(yf) - yf(end)) / (yf(end) - yf(1)) * 100;%sobressinal do gráfico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%especificações do projeto controlador digital
Mpc=(2*length('Matheus'))/100; %sobressinal para o meu nome
Ts5c=length('MatheusSandimGoncalves')*1e-3;%tempo de acomodação para o meu nome completo
zeta2 = sqrt((log(Mpc)^2) / ((pi^2) + log(Mpc)^2));%fator de amortecimento para controlador
Wnc = 3/(zeta2*Ts5c);%frequência natural do controlador
wd=Wnc*sqrt(1-zeta2^2);%frequência amortecida
wa=8.8*wd;%frequência de amostragem
Ta=(2*pi)/wa;% período de anostragem
%cálculo dos pólos
Gz=c2d(Gs,Ta);
z_modulo=exp(-Ta*zeta2*Wnc);% modulo de z
z_ang=Ta*Wnc*sqrt(1-zeta2^2);% angulo de z
z_real=z_modulo*cos(z_ang);
z_img=z_modulo*sin(z_ang);
z1=complex(z_real,z_img);% Z complexo
Cz = tf(Gz.num(1), 1, Ta);% numerador de Gz
Ci = tf(1, [1 -1], Ta);% integrador z-1
G2z=Cz*Ci;%junção dos dois anteriores
angulo_Gc=angle(evalfr(G2z,z1));%obter o angulo de Gc em função de z1
disp(['Ângulo de H(z) em graus: ', num2str(rad2deg(angulo_Gc))]);% angulo do controlador
G1z= -angulo_Gc-pi;%Angulo de G1z
beta=(z_img-(z_real)*tan(-G1z))/tan(-G1z);% calculo do beta
Cp=tf(1,[1 beta],Ta);
FTMA_temp=Cp*G2z;
[Gz_n,Gz_d]=tfdata(Gz, 'v');
Kc=1/abs(evalfr(FTMA_temp,z1));%calculo de Kc
C_z=tf(Kc*Gz_d,[1 beta-1 -beta],Ta);%função de transferencia C (controlador)
FTMA=minreal(C_z*Gz); % função de tranferencia de malha aberta
FTMF=feedback(FTMA,1);% função de transferencia de malha fechada
polos=pole(FTMF);% polos da malha fechada
figure(2)
pzmap(FTMF)% lugar das raizes da função de tranferencia de malha fechada
title('Pólos e zeros de FTMF');
Kmax = ceil(0.1/Ta);
k = (0:Ta:0.1+1);
[yf2, t]=step(FTMF,k);

figure(3)
plot(k, yf2, 'r*');% sistema controlado
hold on;
plot(k, yf2);
title('Resposta do sistema para degrau unitario');
grid on;
hold on;
plot(xf, yf, 'b-');% grafico do sistema sem controle
xlim([0, 0.1]);
legend({'Controlada com *', 'controlada','original'});
hold off
figure (4)
step(FTMF)
figure (5)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%recursiva
kmax=ceil(0.2/Ta);
%tk = Ta * (0:(kmax-1));

r = ones(1, kmax);
r =[r 1.5*r];
%yk = ones(1, kmax);
%e = ones(1, kmax);
kmax=(2*kmax)-1;
tk = Ta * (0:(kmax-1));
yk(1)=0;
e(1) = r(1) - yk(1);
u(1) = 2.576 * e(1);
yk(2) = 0.09157 * u(1) + 1.5 * yk(1);
e(2) = r(2) - yk(2);
u(2) = 2.576 * e(2) - 3.864 * e(1) + 1.203 * u(1);


for k = 3:kmax
    yk(k) = 0.09157 * u(k-1) + 1.5 * yk(k-1) + 0.08015 * u(k-2) - 0.672 * yk(k-2);% bloco G
    e(k) = r(k) - yk(k);
    u(k) = 2.576 * e(k) - 3.864 * e(k-1) + 1.731 * e(k-2) + 1.203 * u(k-1) -0.2027 * u(k-2);%controlador
end

plot(tk, yk);%grafico recursivo
hold on
plot(tk, u);
title('Ação do controlador para uma entrada tipo degrau 1.5V ');
hold off