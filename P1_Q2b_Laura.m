% análise do erro em regime permanante
% entrada tipo rampa
clc
clear
format long

T=0.2;
Gs=tf(1,[1 2 0])

Gz=c2d(Gs,T)

Gc = tf([1 -0.6703], [1 2*0.25 0.25*0.25],T)
FTMA=4.74*Gz*Gc

% importante, erro em regime permanente é calculado a partir 
% da função de transferência de laço aberto
aux=minreal(tf([1 -1],[1 0],T)*FTMA);  % onde tf([1 -1],[1 0],T) = (z-1)/z = 1-z^-1
[n,d]=tfdata(aux,'v')
kv=(polyval(n,1)/polyval(d,1))/T   % substituindo 1 em aux
ess=1/kv

