%% Setup
clear; home;

load('TissDataAPOPAI2021.mat');

%% Example plot of subject 1, thigh, experiment 1
nm=sub1.plotTissNames{4};
data=sub1.(nm)(1);

figure(1000); clf;
errorbar(data.CW.lambda, data.CW.II(:, 1), data.CW.II_err(:, 1)); hold on;
errorbar(data.CW.lambda, data.CW.II(:, 2), data.CW.II_err(:, 2));
errorbar(data.CW.lambda, data.CW.II(:, 3), data.CW.II_err(:, 3));
errorbar(data.CW.lambda, data.CW.II(:, 4), data.CW.II_err(:, 4)); hold off;
set(gca, "YScale", "log");
xlim([data.CW.lambda(1), data.CW.lambda(end)]);
legend('1A', '1B', '2B', '2A');
xlabel('$\lambda$ (nm)', 'Interpreter','latex');
ylabel('$I$ (arb./mm$^2$)', 'Interpreter','latex');

title('Subject 1, Thigh, Experiment 1: CW Data');



figure(1010); clf;
subplot(2, 1, 1);
errorbar(data.FD.lambda, data.FD.II(:, 1), data.FD.II_err(:, 1), 'o'); 
hold on;
errorbar(data.FD.lambda, data.FD.II(:, 2), data.FD.II_err(:, 2), 'o');
errorbar(data.FD.lambda, data.FD.II(:, 3), data.FD.II_err(:, 3), 'o');
errorbar(data.FD.lambda, data.FD.II(:, 4), data.FD.II_err(:, 4), 'o'); 
hold off;
xlim([data.FD.lambda(end), data.FD.lambda(1)]+[-50, 50]);
legend('1A', '1B', '2B', '2A');
xlabel('$\lambda$ (nm)', 'Interpreter','latex');
ylabel('$I$ (arb./mm$^2$)', 'Interpreter','latex');

title('Subject 1, Thigh, Experiment 1: FD Data');


subplot(2, 1, 2);
errorbar(data.FD.lambda, data.FD.PP(:, 1), data.FD.PP_err(:, 1), 'o'); 
hold on;
errorbar(data.FD.lambda, data.FD.PP(:, 2), data.FD.PP_err(:, 2), 'o');
errorbar(data.FD.lambda, data.FD.PP(:, 3), data.FD.PP_err(:, 3), 'o');
errorbar(data.FD.lambda, data.FD.PP(:, 4), data.FD.PP_err(:, 4), 'o'); 
hold off;
xlim([data.FD.lambda(end), data.FD.lambda(1)]+[-50, 50]);
xlabel('$\lambda$ (nm)', 'Interpreter','latex');
ylabel('$\phi$ (rad)', 'Interpreter','latex');