T = 0.5;
z_rect = 0.4 + 0.3*1j;
abs_z_rect = abs(z_rect); % Magnitude of z_rect
phase_z_rect = angle(z_rect); % Phase angle of z_rect
G1z = tf([0.2], [1,-0.6],T);
Hz = tf([1,0], [1,-1],T);
Gc = tf([1,-0.6], [1,-0.25],T);
FTMF = minreal((Gc*G1z)/(1+(Gc*G1z*Hz)));
FTMA = minreal(Gc*G1z*Hz);
Kc = 1 / abs(freqresp(FTMA, z_rect));
zeta = (1 + (log(abs_z_rect) / phase_z_rect)^(-2))^(-1/2);
rlocus(FTMA)
hold on;

% Plot the desired point 'z_rect' with a red 'x' marker
plot(z_rect, 'xr', 'MarkerSize', 10);

% Customize the appearance of the plot
title('Root Locus Plot');
xlabel('Real');
ylabel('Imaginary');
grid on;

% Add a legend for the red 'x' marker
legend('Root Locus', 'Desired Point');

% Optionally, adjust the plot limits for better visualization
% xlim([-2, 2]);
% ylim([-2, 2]);

% Annotate the Kc and zeta values on the plot
text(real(z_rect) + 0.05, imag(z_rect), ['Kc = ', num2str(Kc)], 'Color', 'r');
text(real(z_rect) + 0.05, imag(z_rect) - 0.1, ['zeta = ', num2str(zeta)], 'Color', 'r');

hold off;