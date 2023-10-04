clc
clear
 
% Período de amostragem
T = 0.5;

% Função de Transferência Contínua
Gs=tf([1 1],[1 4])
Hs=tf(2,[1 0])
GHs=Gs*Hs

% Função de Transferência Discreta
Gz=c2d(Gs,T)
GHz=c2d(GHs,T)

k = 0:25;       % kmax = 25

% Entrada: Rampa Unitária
r = k*T;

% Condição Inicial
c(1) = 0;   % para k = 0
e(1) = 0;
v(1) = 0;
r(1) = 0;

 for j=1:length(k)-2
     v(j+2)=0.5742*e(j+1)-0.3581*e(j)+1.135*v(j+1)-0.1353*v(j);
     e(j)=r(j)-v(j);
     c(j+1)=e(j+1)-0.7838*e(j)+0.1353*c(j);
 end 
 plot(k*T,c,'ok')
 hold on


% for j=3:length(k)
%      v(j)=0.5742*e(j-1)-0.3581*e(j-2)+1.135*v(j-1)-0.1353*v(j-2);
%      e(j)=r(j)-v(j);
%      c(j)=e(j) - 0.7838*e(j-1) + 0.1353*c(j-1);
% end
%  
% plot(k*T,c,'ok')

t = 0:T:25*T;

FTMF = minreal(Gz/(1+GHz))
x = lsim(FTMF, t);
plot(k*T, x)
hold off
