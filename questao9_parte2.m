T = 0.2;
A = tf(1, [1, 2, 0]);
Az = c2d(A, T)
pole(Az)
G2z = tf([0.01758, 0.01539], [1,-1], T)
T=0.2;
F = tf(1, [1, 2, 0]);
Fz = c2d(F, T)
zeta = 0.6;
wn = 4;
z_ang = T * wn * sqrt(1 - (zeta^2));
fprintf('Ângulo de Z: %f\n', z_ang);

z_mod = exp(-T * zeta * wn);
fprintf('Módulo de Z: %f\n', z_mod);

% Transformando agora para a forma retangular
z_rect = z_mod * exp(1i * z_ang);
fprintf('Z na forma retangular: %f + %fi\n', real(z_rect), imag(z_rect));

angle_g1 = angle(freqresp(G2z, z_rect)) - pi;
fprintf('Ângulo em radianos: %f\n', angle_g1);
angle_degrees = rad2deg(angle_g1);
fprintf('Ângulo em graus: %f\n', angle_degrees);

tan_angle_g1 = tan(angle_g1);
beta = (imag(z_rect) - (real(z_rect) * tan_angle_g1)) / tan_angle_g1;
fprintf('Valor de beta calculado: %f\n', beta);

G1z=tf([1],[1, beta],T);
FTMA= G1z*G2z;
Kc= 1/abs(freqresp(FTMA, z_rect));
rlocus(FTMA, Kc);
hold on;

% Overlay a marker at the point z_rect
plot(z_rect, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Customize plot labels and legend
xlabel('Real Part');
ylabel('Imaginary Part');
title('Root Locus Plot');
legend('Root Locus', 'Z-rect Point');

% Adjust the plot limits for better visualization
xlim([-2, 2]);
ylim([-2, 2]);

% Display the grid
grid on;

% Customize line colors and styles for the root locus
line_handles = findobj(gcf, 'Type', 'line');
set(line_handles(1:end-1), 'Color', 'b', 'LineStyle', '-');

% Customize marker properties for the z_rect point
marker_handles = findobj(gcf, 'Type', 'scatter');
set(marker_handles, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'SizeData', 100);

hold off;
