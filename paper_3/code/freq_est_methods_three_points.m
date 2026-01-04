clc; clear; close all;

%% Signal parameters
fs = 200e9;          % sampling freq
f0 = 71e9;           % true tone freq (Hz)
N = 256;
t = (0:N-1)'/fs;

%% SNR sweep settings
SNR_vec = 200:2:240;    % SNR range (dB)
MC = 2000;              % Monte-Carlo runs

%% Arrays to store RMSE
rmse_parabolic = zeros(size(SNR_vec));
rmse_jacobsen  = zeros(size(SNR_vec));
rmse_macleod   = zeros(size(SNR_vec));
rmse_candan    = zeros(size(SNR_vec));
rmse_m3        = zeros(size(SNR_vec));
rmse_m4        = zeros(size(SNR_vec));

%% Loop over SNR values
for ss = 1:length(SNR_vec)
    SNR_dB = SNR_vec(ss);
    
    err_p = zeros(MC,1);
    err_j = zeros(MC,1);
    err_m = zeros(MC,1);
    err_c = zeros(MC,1);
    err_3 = zeros(MC,1);
    err_4 = zeros(MC,1);

    for mc = 1:MC
        
        %% Generate clean tone
        x_clean = 0.9*sin(2*pi*f0*t + 0.3);

        %% Add noise
        Px = mean(abs(x_clean).^2);             
        Pn = Px / (10^(SNR_dB/10));             
        noise = sqrt(Pn) * randn(size(x_clean));
        x = x_clean + noise; 

        if isreal(x), x = hilbert(x); end

        %% FFT
        X = fft(x,N);
        [~, krel] = max(abs(X));
        k = krel - 1;

        wrap = @(k) mod(k, N);

        R_m1 = X(wrap(k-1)+1);
        R_0  = X(wrap(k)+1);
        R_p1 = X(wrap(k+1)+1);

        %% ===== Method 1: Parabolic interpolation ======
        delta1 = (abs(R_p1)-abs(R_m1)) / (4*abs(R_0) - 2*abs(R_m1) - 2*abs(R_p1));
        f1 = (k + delta1)*fs/N;

        %% ===== Method 2: Jacobsen ======
        delta2 = real((R_m1 - R_p1)/(2*R_0 - R_m1 - R_p1));
        f2 = (k + delta2)*fs/N;

        %% ===== Method 3: MacLeod ======
        d = real(R_m1*conj(R_0)-R_p1*conj(R_0)) / ...
            real(2*(abs(R_0))^2 + R_m1*conj(R_0) + R_p1*conj(R_0));
        delta5 = (sqrt(1+8*d^2)-1)/(4*d);
        f5 = (k + delta5)*fs/N;

        %% ===== Method 4: Candan ======
        f7 = candan_three_point_freq(x, fs);

        %% ===== Method 5: Method-3 (Hamming) ======
        w = hamming(N);
        xw = x .* w;
        Xw = fft(xw,N);
        [~, krel] = max(abs(Xw));
        k = krel - 1;
        R_m1 = Xw(wrap(k-1)+1);
        R_0  = Xw(wrap(k)+1);
        R_p1 = Xw(wrap(k+1)+1);
        P = 1.22;   % Hamming constant
        delta3 = P*(abs(R_p1)-abs(R_m1))/(abs(R_0)+abs(R_p1)+abs(R_m1));
        f3 = (k + delta3)*fs/N;

        %% ===== Method 6: Method-4 (Hamming) ======
        Q = 0.60;   % Hamming constant
        delta4 = real(Q*(R_p1-R_m1)/(2*R_0+R_m1+R_p1));
        f4 = (k + delta4)*fs/N;

        %% Store errors
        err_p(mc) = (f1 - f0)^2;
        err_j(mc) = (f2 - f0)^2;
        err_m(mc) = (f5 - f0)^2;
        err_c(mc) = (f7 - f0)^2;
        err_3(mc) = (f3 - f0)^2;
        err_4(mc) = (f4 - f0)^2;
    end

    %% Compute RMSE for this SNR
    rmse_parabolic(ss) = sqrt(mean(err_p));
    rmse_jacobsen(ss)  = sqrt(mean(err_j));
    rmse_macleod(ss)   = sqrt(mean(err_m));
    rmse_candan(ss)    = sqrt(mean(err_c));
    rmse_m3(ss)        = sqrt(mean(err_3));
    rmse_m4(ss)        = sqrt(mean(err_4));
end

%% ===== Plot RMSE curves ======
figure; hold on; grid on;

% Define markers and colors
markers = {'o','s','d','^','v','>'};
colors  = lines(6);   % generates 6 visually distinct colors

% Plot with unique color + marker
plot(SNR_vec, rmse_parabolic, '-o', 'Color', colors(1,:), 'MarkerSize', 6, 'LineWidth', 1.8);hold on;
plot(SNR_vec, rmse_jacobsen,  '-s', 'Color', colors(2,:), 'MarkerSize', 6, 'LineWidth', 1.8);hold on;
plot(SNR_vec, rmse_macleod,   '-d', 'Color', colors(3,:), 'MarkerSize', 6, 'LineWidth', 1.8);hold on;
plot(SNR_vec, rmse_candan,    '-^', 'Color', colors(4,:), 'MarkerSize', 6, 'LineWidth', 1.8);hold on;
plot(SNR_vec, rmse_m3,        '-v', 'Color', colors(5,:), 'MarkerSize', 6, 'LineWidth', 1.8);hold on;
plot(SNR_vec, rmse_m4,        '->', 'Color', colors(6,:), 'MarkerSize', 6, 'LineWidth', 1.8);hold on;

xlabel('SNR (dB)');
ylabel('RMSE of Frequency Estimation (Hz)');
title('Frequency Estimation RMSE vs SNR');

legend('Parabolic', 'Jacobsen', 'MacLeod', 'Candan', ...
       'Method-3 (Ham)', 'Method-4 (Ham)', 'Location','northEast');

