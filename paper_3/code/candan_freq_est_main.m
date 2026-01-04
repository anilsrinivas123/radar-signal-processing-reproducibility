% mc_candan_vs_crlb.m
% Monte Carlo: Candan 3-DFT-sample estimator vs CRLB
% SNR: -10 dB .. +10 dB, 1000 trials per SNR, N=256
clear; close all; clc;

% Parameters

f0 = 1234.567;        % true tone (Hz)
fs = 2.1*f0;            % sampling freq (Hz)
N = 1024;              % # samples
omega0 = 2*pi*f0/fs;  % radians/sample
snr_db_vec = -10:1:10;
n_trials = 10000;


% Precompute time indices
n = (0:N-1)';

% Storage
mse_freq = zeros(size(snr_db_vec));
mse_freq_trials = zeros(length(snr_db_vec), n_trials);

fprintf('Starting Monte Carlo: N=%d, trials=%d per SNR\n', N, n_trials);
tic;
for ii = 1:length(snr_db_vec)
    snr_db = snr_db_vec(ii);
    snr_lin = 10^(snr_db/10);
    % Signal amplitude A such that signal power per sample = 1
    A = 1; 
    noise_var = A^2 / snr_lin;      % E[|noise|^2] per complex sample
    for tr = 1:n_trials
        phase = 2*pi*rand; % random initial phase
        s = A * exp(1j*(omega0*n + phase));
        % complex AWGN: real and imag with variance = noise_var/2
        noise = sqrt(noise_var/2) * (randn(N,1) + 1j*randn(N,1));
        x = s + noise;
        fhat = candan_three_point_freq(x, fs);  % returns Hz
        mse_freq_trials(ii, tr) = (fhat - f0)^2;
    end
    mse_freq(ii) = mean(mse_freq_trials(ii, :));
    if mod(ii,5)==0
        fprintf('  Done SNR = %+3d dB (%d/%d)\n', snr_db, ii, length(snr_db_vec));
    end
end
toc;

% CRLB (angular frequency variance) and convert to Hz^2
snr_lin_vec = 10.^(snr_db_vec/10);
crlb_omega = 12 ./ (snr_lin_vec .* N .* (N^2 - 1));      % Var(omega) (rad^2/sample^2)
crlb_freq = crlb_omega .* (fs/(2*pi)).^2;              % Var(freq) (Hz^2)

% Plot MSE vs CRLB
figure('Color','w','Units','normalized','Position',[0.05 0.1 0.6 0.6]);
semilogy(snr_db_vec, mse_freq, 'o-', 'LineWidth', 1.5, 'MarkerSize',6);
hold on;
semilogy(snr_db_vec, crlb_freq, 'x--', 'LineWidth', 1.5, 'MarkerSize',6);
grid on;
xlabel('SNR (dB)');
ylabel('MSE (Hz^2) (log scale)');
title(sprintf('Candan 3-point estimator vs CRLB (N=%d, trials=%d)', N, n_trials));
legend('Candan (MC MSE)','CRLB (Hz^2)','Location','southwest');
set(gca,'FontSize',11);

% Plot ratio in dB
ratio_db = 10*log10(mse_freq ./ (crlb_freq + eps));
figure('Color','w','Units','normalized','Position',[0.7 0.1 0.25 0.5]);
plot(snr_db_vec, ratio_db, 'o-', 'LineWidth',1.5, 'MarkerSize',6);
grid on;
xlabel('SNR (dB)');
ylabel('10 log10( MSE / CRLB ) (dB)');
title('Performance gap to CRLB');
set(gca,'FontSize',11);

% Print a few summary points
fprintf('\nSummary examples:\n');
example_idx = [1, 6, 11, length(snr_db_vec)]; % -10, -5, 0, +10 dB indices
for k = example_idx
    fprintf(' SNR %+3d dB: MSE_freq = %.3e Hz^2, CRLB_freq = %.3e Hz^2, ratio = %.2f\n', ...
        snr_db_vec(k), mse_freq(k), crlb_freq(k), mse_freq(k)/crlb_freq(k));
end