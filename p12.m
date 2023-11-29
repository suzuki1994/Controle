% Exemplo de cálculo alternativo da condição de angulo
% para determinar o polo do controlador
clear
clc

% periodo de amostragem
T=0.15;

% função de transferência da planta
Gp=tf(1,[1 5]);
Gz=c2d(Gp,T);

% polo desejado
zeta=0.6;
wn=4;
s1=-zeta*wn+j*wn*sqrt(1-zeta^2);
z1=exp(s1*T);

% Controlador
% Gd(z) = Kc*(z+alfa)/(z+beta)

% Função de transferência de laço aberto
% FTMA=Gdz*Gz

% condição de angulo --> angulo(FTMA)=+-180 graus quando z =z1 
% o ganho do controlador não interfere no angulo

% dividindo a FTMA de forma diferente
% FTMA = G1z*G2z
% onde G1z=1/(z+beta)
% e G2z=(z+alfa)*Gpz

% considalfaerando que o zero do controlador (alfa) cancela um polo da planta
polos=pole(Gz)
alfa = -polos(1)

G2z=minreal(tf([1 alfa],1,T)*Gz)
[n2,d2]=tfdata(G2z,'v')

% fi2 é o ângulo de G2z quando z=z1
fi2=angle(polyval(n2,z1)/polyval(d2,z1))
% fi1 é o ângulo de G1z quando z=z1
fi1=-pi-fi2;
beta=(imag(z1)-real(z1)*tan(-fi1))/tan(-fi1)

% calculo k
