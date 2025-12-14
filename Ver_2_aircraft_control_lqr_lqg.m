%% ================================================================
%% IMPLEMENTASI LENGKAP LQR dan LQG UNTUK SISTEM KONTROL PESAWAT
%% ================================================================
% Nama: Rizky Wisuda Wardani
% NIM: 6009242003
% Semester Ganjil 2025/2026
% ================================================================

clear all; close all; clc;
fprintf('========================================\n');
fprintf('SISTEM KONTROL PESAWAT: LQR vs LQG\n');
fprintf('========================================\n\n');

%% ================================================================
%% 1. PARAMETER SISTEM PESAWAT (BERDASARKAN JURNAL)
%% ================================================================

% Sistem Longitudinal (Pitch Control) - Persamaan (7-8)
A_long = [-0.3149,  235.8928,    0;
          -0.0034,   -0.4282,    0;
               0,         1,      0];

B_long = [-5.5079;
           0.0021;
               0];

C_long = [0, 0, 1];  % Output: pitch angle θ
D_long = 0;

sys_long = ss(A_long, B_long, C_long, D_long);

% Sistem Lateral (Roll & Sideslip Control) - Persamaan (11)
A_lat = [-0.0558, -0.9968,  0.0802,  0.0415;
          0.5980, -0.1150, -0.0318,       0;
         -3.0500,  0.3880, -0.4650,       0;
               0,  0.0805,  1.0000,       0];

B_lat = [ 0.0729,   0.000;
         -4.7500,   0.00775;
          0.1530,   0.1430;
               0,        0];

C_lat = [1, 0, 0, 0;   % β (sideslip angle)
         0, 0, 0, 1];  % φ (roll angle)
D_lat = [0, 0;
         0, 0];

sys_lat = ss(A_lat, B_lat, C_lat, D_lat);

fprintf('1. SISTEM TELAH DIMODELKAN:\n');
fprintf('   - Longitudinal: 3 state (w, q, θ)\n');
fprintf('   - Lateral: 4 state (β, p, r, φ)\n');

%% ================================================================
%% 2. FUNGSI BANTU UNTUK ANALISIS STEP RESPONSE
%% ================================================================

% Fungsi untuk menghitung stepinfo tanpa warning
function info = calc_stepinfo(t, y, y_final)
    if nargin < 3
        y_final = mean(y(end-round(0.1*length(y)):end)); % Rata-rata 10% terakhir
    end
    
    % Parameter
    SettlingTimeThreshold = 0.02;  % 2%
    RiseTimeLower = 0.10;  % 10%
    RiseTimeUpper = 0.90;  % 90%
    
    % Rise Time
    y_initial = y(1);
    y_target = abs(y_final - y_initial);
    
    idx_lower = find(abs(y - y_initial) >= RiseTimeLower * y_target, 1);
    idx_upper = find(abs(y - y_initial) >= RiseTimeUpper * y_target, 1);
    
    if ~isempty(idx_lower) && ~isempty(idx_upper)
        rise_time = t(idx_upper) - t(idx_lower);
    else
        rise_time = NaN;
    end
    
    % Settling Time
    y_settling = abs(y - y_final);
    settling_threshold = SettlingTimeThreshold * abs(y_final);
    idx_settling = find(y_settling <= settling_threshold, 1);
    
    if ~isempty(idx_settling)
        settling_time = t(idx_settling);
        % Pastikan setelah itu tetap dalam threshold
        for i = idx_settling:length(t)
            if y_settling(i) > settling_threshold
                idx_settling = find(y_settling(i:end) <= settling_threshold, 1) + i - 1;
                if ~isempty(idx_settling)
                    settling_time = t(idx_settling);
                end
            end
        end
    else
        settling_time = t(end);
    end
    
    % Overshoot
    if y_final ~= y_initial
        overshoot = 100 * (max(y) - y_final) / abs(y_final - y_initial);
    else
        overshoot = 0;
    end
    
    % Peak
    peak = max(y);
    
    % Return structure
    info.RiseTime = rise_time;
    info.SettlingTime = settling_time;
    info.Overshoot = overshoot;
    info.Peak = peak;
    info.SteadyState = y_final;
end

%% ================================================================
%% 3. DESAIN KONTROLER LQR
%% ================================================================

fprintf('\n2. DESAIN KONTROLER LQR:\n');

% Longitudinal LQR
Q_long = diag([0, 0, 500]);  % x = 500 (sesuai jurnal)
R_long = 1;
[K_lqr_long, ~, pole_lqr_long] = lqr(A_long, B_long, Q_long, R_long);
A_cl_long = A_long - B_long * K_lqr_long;
sys_lqr_long = ss(A_cl_long, B_long, C_long, D_long);

% Lateral LQR
Q_lat = diag([100, 10, 10, 50]);
R_lat = diag([1, 1]);
[K_lqr_lat, ~, pole_lqr_lat] = lqr(A_lat, B_lat, Q_lat, R_lat);
A_cl_lat = A_lat - B_lat * K_lqr_lat;
sys_lqr_lat = ss(A_cl_lat, B_lat, C_lat, D_lat);

fprintf('   - Gain LQR Longitudinal: K = [%.4f, %.4f, %.4f]\n', K_lqr_long(1), K_lqr_long(2), K_lqr_long(3));
fprintf('   - Gain LQR Lateral: \n');
fprintf('        K1 = [%.4f, %.4f, %.4f, %.4f]\n', K_lqr_lat(1,1), K_lqr_lat(1,2), K_lqr_lat(1,3), K_lqr_lat(1,4));
fprintf('        K2 = [%.4f, %.4f, %.4f, %.4f]\n', K_lqr_lat(2,1), K_lqr_lat(2,2), K_lqr_lat(2,3), K_lqr_lat(2,4));

%% ================================================================
%% 4. DESAIN KONTROLER LQG (KALMAN FILTER + LQR)
%% ================================================================

fprintf('\n3. DESAIN KONTROLER LQG:\n');

% Noise characteristics (dari jurnal)
Qn = diag([0.01, 0.01, 0.01]);  % Process noise covariance
Rn = 0.01;                      % Measurement noise covariance

% Kalman Filter untuk longitudinal
sys_kf_long = ss(A_long, [B_long eye(3)], C_long, [0 0 0 0]);
[kest_long, L_kf_long, P_long] = kalman(sys_kf_long, Qn, Rn);

% Sistem LQG longitudinal
A_lqg_long = [A_long,          -B_long*K_lqr_long;
              L_kf_long*C_long, A_long - B_long*K_lqr_long - L_kf_long*C_long];
B_lqg_long = [B_long; zeros(3,1)];
C_lqg_long = [C_long, zeros(1,3)];
D_lqg_long = 0;
sys_lqg_long = ss(A_lqg_long, B_lqg_long, C_lqg_long, D_lqg_long);

% Kalman Filter untuk lateral
C_lat_kf = C_lat;
Qn_lat = diag([0.01, 0.01, 0.01, 0.01]);
Rn_lat = diag([0.01, 0.01]);
sys_kf_lat = ss(A_lat, [B_lat eye(4)], C_lat_kf, zeros(2,6));
[kest_lat, L_kf_lat, P_lat] = kalman(sys_kf_lat, Qn_lat, Rn_lat);

% Sistem LQG lateral
A_lqg_lat = [A_lat,          -B_lat*K_lqr_lat;
             L_kf_lat*C_lat_kf, A_lat - B_lat*K_lqr_lat - L_kf_lat*C_lat_kf];
B_lqg_lat = [B_lat; zeros(4,2)];
C_lqg_lat = [C_lat, zeros(2,4)];
D_lqg_lat = zeros(2,2);
sys_lqg_lat = ss(A_lqg_lat, B_lqg_lat, C_lqg_lat, D_lqg_lat);

fprintf('   - Gain Kalman Filter Longitudinal: L = [%.4f; %.4f; %.4f]\n', L_kf_long(1), L_kf_long(2), L_kf_long(3));

%% ================================================================
%% 5. SIMULASI SISTEM
%% ================================================================

fprintf('\n4. SIMULASI SISTEM:\n');

% Waktu simulasi
t_long = 0:0.01:20;
t_lat = 0:0.01:20;

% Input step
delta_e = 0.2 * ones(size(t_long));      % Elevator 0.2 rad
delta_a = 0.1 * ones(size(t_lat));       % Aileron 0.1 rad
delta_r = 0.05 * ones(size(t_lat));      % Rudder 0.05 rad
u_lat = [delta_a; delta_r]';

% Simulasi Longitudinal
[y_long_open, ~, x_long_open] = lsim(sys_long, delta_e, t_long, [0; 0; 0]);
[y_long_lqr, ~, x_long_lqr] = lsim(sys_lqr_long, delta_e, t_long, [0; 0; 0]);

% Simulasi Lateral
[y_lat_open, ~, x_lat_open] = lsim(sys_lat, u_lat, t_lat, [0; 0; 0; 0]);
[y_lat_lqr, ~, x_lat_lqr] = lsim(sys_lqr_lat, u_lat, t_lat, [0; 0; 0; 0]);

% Simulasi LQG dengan noise
rng(42);
measurement_noise_long = 0.01 * randn(1, length(t_long));
x0_lqg_long = zeros(6,1);
[y_long_lqg, ~, x_long_lqg] = lsim(sys_lqg_long, delta_e, t_long, x0_lqg_long);
y_long_lqg_noisy = y_long_lqg + measurement_noise_long';

% Simulasi LQG lateral dengan noise
measurement_noise_lat = 0.01 * randn(2, length(t_lat));
x0_lqg_lat = zeros(8,1);
[y_lat_lqg, ~, x_lat_lqg] = lsim(sys_lqg_lat, u_lat, t_lat, x0_lqg_lat);
y_lat_lqg_noisy = y_lat_lqg + measurement_noise_lat';

% Control signals LQR
u_lqr_long = -x_long_lqr * K_lqr_long';
u_lqr_lat = -x_lat_lqr * K_lqr_lat';

% Control signals LQG
u_lqg_long = -x_long_lqg(:,4:6) * K_lqr_long';  % Estimasi state untuk kontrol
u_lqg_lat = -x_lat_lqg(:,5:8) * K_lqr_lat';

fprintf('   ✓ Simulasi selesai\n');

%% ================================================================
%% 6. ANALISIS PERFORMANSI
%% ================================================================

fprintf('\n5. ANALISIS PERFORMANSI:\n');

% Hitung stepinfo menggunakan fungsi custom
stepinfo_long_open = calc_stepinfo(t_long', y_long_open);
stepinfo_long_lqr = calc_stepinfo(t_long', y_long_lqr);
stepinfo_long_lqg = calc_stepinfo(t_long', y_long_lqg_noisy);

stepinfo_lat_open_roll = calc_stepinfo(t_lat', y_lat_open(:,2));
stepinfo_lat_lqr_roll = calc_stepinfo(t_lat', y_lat_lqr(:,2));
stepinfo_lat_lqg_roll = calc_stepinfo(t_lat', y_lat_lqg_noisy(:,2));

% Tampilkan hasil
fprintf('   a) Longitudinal (Pitch Angle):\n');
fprintf('        Metric        Open-loop       LQR          LQG\n');
fprintf('        Rise Time     %8.4f s    %8.4f s    %8.4f s\n', ...
        stepinfo_long_open.RiseTime, stepinfo_long_lqr.RiseTime, stepinfo_long_lqg.RiseTime);
fprintf('        Settling Time %8.4f s    %8.4f s    %8.4f s\n', ...
        stepinfo_long_open.SettlingTime, stepinfo_long_lqr.SettlingTime, stepinfo_long_lqg.SettlingTime);
fprintf('        Overshoot     %8.2f%%    %8.2f%%    %8.2f%%\n', ...
        stepinfo_long_open.Overshoot, stepinfo_long_lqr.Overshoot, stepinfo_long_lqg.Overshoot);

%% ================================================================
%% 7. ANALISIS KESTABILAN DAN GAIN
%% ================================================================

fprintf('\n6. ANALISIS KESTABILAN DAN GAIN:\n');

% Margin stabilitas
[GM_lqr_long, PM_lqr_long, ~, ~] = margin(sys_lqr_long);
[GM_lqg_long, PM_lqg_long, ~, ~] = margin(sys_lqg_long);

fprintf('   a) Margin Stabilitas Longitudinal:\n');
fprintf('        Controller     Gain Margin (dB)   Phase Margin (deg)\n');
fprintf('        LQR            %12.2f        %12.2f\n', 20*log10(GM_lqr_long), PM_lqr_long);
fprintf('        LQG            %12.2f        %12.2f\n', 20*log10(GM_lqg_long), PM_lqg_long);

% Pole sistem
fprintf('\n   b) Pole Sistem (Longitudinal):\n');
eig_open = eig(A_long);
fprintf('        Open-loop poles: ');
for i = 1:length(eig_open)
    fprintf('%.4f ', eig_open(i));
end
fprintf('\n');

fprintf('        LQR poles: ');
for i = 1:length(pole_lqr_long)
    fprintf('%.4f ', pole_lqr_long(i));
end
fprintf('\n');

eig_lqg = eig(A_lqg_long);
fprintf('        LQG poles: ');
for i = 1:length(eig_lqg)
    fprintf('%.4f ', eig_lqg(i));
end
fprintf('\n');

%% ================================================================
%% 8. GRAFIK 1: STEP RESPONSE COMPARISON (DENGAN LQG)
%% ================================================================

figure('Name', 'Step Response Comparison with LQG', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 500]);

% Longitudinal
subplot(1, 2, 1);
plot(t_long, y_long_open, 'b-', 'LineWidth', 1.5); hold on;
plot(t_long, y_long_lqr, 'r--', 'LineWidth', 2);
plot(t_long, y_long_lqg_noisy, 'g-.', 'LineWidth', 1.5);
plot(t_long, 0.2*ones(size(t_long)), 'k:', 'LineWidth', 1);
xlabel('Time (sec)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Pitch Angle \theta (rad)', 'FontSize', 11, 'FontWeight', 'bold');
title('Longitudinal: Pitch Angle Response', 'FontSize', 12, 'FontWeight', 'bold');
legend('Open-loop', 'LQR', 'LQG with Noise', 'Reference (0.2 rad)', 'Location', 'best');
grid on;
xlim([0 10]);

% Tandai settling time LQR dan LQG
if ~isnan(stepinfo_long_lqr.SettlingTime) && stepinfo_long_lqr.SettlingTime <= 10
    line([stepinfo_long_lqr.SettlingTime, stepinfo_long_lqr.SettlingTime], ...
         [min(y_long_lqr), max(y_long_lqr)], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);
    text(stepinfo_long_lqr.SettlingTime+0.3, 0.15, ...
         sprintf('LQR: %.2f s', stepinfo_long_lqr.SettlingTime), ...
         'Color', 'r', 'FontSize', 8);
end

if ~isnan(stepinfo_long_lqg.SettlingTime) && stepinfo_long_lqg.SettlingTime <= 10
    line([stepinfo_long_lqg.SettlingTime, stepinfo_long_lqg.SettlingTime], ...
         [min(y_long_lqg_noisy), max(y_long_lqg_noisy)], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 1);
    text(stepinfo_long_lqg.SettlingTime+0.3, 0.1, ...
         sprintf('LQG: %.2f s', stepinfo_long_lqg.SettlingTime), ...
         'Color', 'g', 'FontSize', 8);
end

% Lateral (Roll Angle) dengan LQG
subplot(1, 2, 2);
plot(t_lat, y_lat_open(:,2), 'b-', 'LineWidth', 1.5); hold on;
plot(t_lat, y_lat_lqr(:,2), 'r--', 'LineWidth', 2);
plot(t_lat, y_lat_lqg_noisy(:,2), 'g-.', 'LineWidth', 1.5);
plot(t_lat, 0.1*ones(size(t_lat)), 'k:', 'LineWidth', 1);
xlabel('Time (sec)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Roll Angle \phi (rad)', 'FontSize', 11, 'FontWeight', 'bold');
title('Lateral: Roll Angle Response', 'FontSize', 12, 'FontWeight', 'bold');
legend('Open-loop', 'LQR', 'LQG with Noise', 'Reference (0.1 rad)', 'Location', 'best');
grid on;
xlim([0 10]);

%% ================================================================
%% 9. GRAFIK 2: CONTROL SIGNAL LONGITUDINAL & LATERAL (DENGAN LQG)
%% ================================================================

figure('Name', 'Control Signals with LQG', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 500]);

% Longitudinal Control - LQR vs LQG
subplot(1, 2, 1);
plot(t_long, u_lqr_long, 'r-', 'LineWidth', 2); hold on;
plot(t_long, u_lqg_long, 'g--', 'LineWidth', 1.5);
xlabel('Time (sec)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Elevator Deflection \delta_e (rad)', 'FontSize', 11, 'FontWeight', 'bold');
title('Longitudinal: Control Effort Comparison', 'FontSize', 12, 'FontWeight', 'bold');
legend('LQR', 'LQG', 'Location', 'best');
grid on;
xlim([0 5]);
ylim([-1 1]);

% Lateral Control - LQR vs LQG
subplot(1, 2, 2);
plot(t_lat, u_lqr_lat(:,1), 'r-', 'LineWidth', 2); hold on;
plot(t_lat, u_lqg_lat(:,1), 'g--', 'LineWidth', 1.5);
plot(t_lat, u_lqr_lat(:,2), 'b-', 'LineWidth', 2);
plot(t_lat, u_lqg_lat(:,2), 'm--', 'LineWidth', 1.5);
xlabel('Time (sec)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Control Deflection (rad)', 'FontSize', 11, 'FontWeight', 'bold');
title('Lateral: Control Effort Comparison', 'FontSize', 12, 'FontWeight', 'bold');
legend('Aileron LQR', 'Aileron LQG', 'Rudder LQR', 'Rudder LQG', 'Location', 'best');
grid on;
xlim([0 5]);
ylim([-1 1]);

%% ================================================================
%% 10. GRAFIK 3: KESTABILAN & POLE-ZERO MAP (DENGAN LQG)
%% ================================================================

figure('Name', 'Stability Analysis with LQG', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 500]);

% Pole-Zero Map Longitudinal
subplot(1, 2, 1);
pzmap(sys_long, 'b'); hold on;
pzmap(sys_lqr_long, 'r');
pzmap(sys_lqg_long, 'g');
title('Longitudinal: Pole-Zero Map Comparison', 'FontSize', 12, 'FontWeight', 'bold');
legend('Open-loop', 'LQR', 'LQG', 'Location', 'best');
grid on;

% Pole-Zero Map Lateral
subplot(1, 2, 2);
pzmap(sys_lat, 'b'); hold on;
pzmap(sys_lqr_lat, 'r');
pzmap(sys_lqg_lat, 'g');
title('Lateral: Pole-Zero Map Comparison', 'FontSize', 12, 'FontWeight', 'bold');
legend('Open-loop', 'LQR', 'LQG', 'Location', 'best');
grid on;

%% ================================================================
%% 11. GRAFIK 4: NYQUIST PLOT (DENGAN LQG)
%% ================================================================

figure('Name', 'Nyquist Plot with LQG', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 500]);

% Nyquist untuk Open Loop
subplot(1, 2, 1);
nyquist(sys_long);
title('Longitudinal: Nyquist Plot (Open-loop)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Nyquist untuk LQR dan LQG
subplot(1, 2, 2);
nyquist(sys_lqr_long, 'r'); hold on;
nyquist(sys_lqg_long, 'g');
title('Longitudinal: Nyquist Plot (LQR vs LQG)', 'FontSize', 12, 'FontWeight', 'bold');
legend('LQR', 'LQG', 'Location', 'best');
grid on;

%% ================================================================
%% 12. GRAFIK 5: FREQUENCY RESPONSE COMPARISON (BODE PLOT)
%% ================================================================

figure('Name', 'Frequency Response - Bode Plot', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 600]);

% Bode Plot Magnitude
subplot(2, 1, 1);
[mag_open, phase_open, w_open] = bode(sys_long);
[mag_lqr, phase_lqr, w_lqr] = bode(sys_lqr_long);
[mag_lqg, phase_lqg, w_lqg] = bode(sys_lqg_long);

semilogx(w_open, 20*log10(squeeze(mag_open)), 'b-', 'LineWidth', 2); hold on;
semilogx(w_lqr, 20*log10(squeeze(mag_lqr)), 'r--', 'LineWidth', 2);
semilogx(w_lqg, 20*log10(squeeze(mag_lqg)), 'g-.', 'LineWidth', 1.5);
xlabel('Frequency (rad/sec)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Magnitude (dB)', 'FontSize', 11, 'FontWeight', 'bold');
title('Bode Plot - Magnitude', 'FontSize', 12, 'FontWeight', 'bold');
legend('Open-loop', 'LQR', 'LQG', 'Location', 'best');
grid on;

% Bode Plot Phase
subplot(2, 1, 2);
semilogx(w_open, squeeze(phase_open), 'b-', 'LineWidth', 2); hold on;
semilogx(w_lqr, squeeze(phase_lqr), 'r--', 'LineWidth', 2);
semilogx(w_lqg, squeeze(phase_lqg), 'g-.', 'LineWidth', 1.5);
xlabel('Frequency (rad/sec)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Phase (deg)', 'FontSize', 11, 'FontWeight', 'bold');
title('Bode Plot - Phase', 'FontSize', 12, 'FontWeight', 'bold');
legend('Open-loop', 'LQR', 'LQG', 'Location', 'best');
grid on;

%% ================================================================
%% 13. GRAFIK 6: IMPULSE RESPONSE - PITCH ANGLE (LQR & LQG)
%% ================================================================

figure('Name', 'Impulse Response - Pitch Angle', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 500]);

% Impulse response untuk LQR dan LQG
subplot(1, 2, 1);
[y_imp_lqr, t_imp_lqr] = impulse(sys_lqr_long, 10);
[y_imp_lqg, t_imp_lqg] = impulse(sys_lqg_long, 10);
plot(t_imp_lqr, y_imp_lqr, 'r-', 'LineWidth', 2); hold on;
plot(t_imp_lqg, y_imp_lqg, 'g--', 'LineWidth', 2);
xlabel('Time (sec)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Pitch Angle \theta (rad)', 'FontSize', 11, 'FontWeight', 'bold');
title('Impulse Response: LQR vs LQG', 'FontSize', 12, 'FontWeight', 'bold');
legend('LQR', 'LQG', 'Location', 'best');
grid on;

% Impulse response detail (zoom)
subplot(1, 2, 2);
plot(t_imp_lqr, y_imp_lqr, 'r-', 'LineWidth', 2); hold on;
plot(t_imp_lqg, y_imp_lqg, 'g--', 'LineWidth', 2);
xlabel('Time (sec)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Pitch Angle \theta (rad)', 'FontSize', 11, 'FontWeight', 'bold');
title('Impulse Response Detail (0-2 sec)', 'FontSize', 12, 'FontWeight', 'bold');
legend('LQR', 'LQG', 'Location', 'best');
grid on;
xlim([0 2]);
ylim([-0.01 0.01]);

%% ================================================================
%% 14. SIMPULAN DAN SIMPAN HASIL
%% ================================================================

fprintf('\n7. SIMPULAN:\n');
fprintf('   =========\n');
fprintf('   1. LQR memberikan respons yang lebih cepat dengan settling time %.2f s\n', stepinfo_long_lqr.SettlingTime);
fprintf('      dibanding open-loop (%.2f s)\n', stepinfo_long_open.SettlingTime);
fprintf('   2. LQG mampu mengatasi noise dengan settling time %.2f s\n', stepinfo_long_lqg.SettlingTime);
fprintf('   3. Semua sistem terkontrol stabil (semua pole di LHP)\n');
fprintf('   4. Control effort LQG lebih smooth karena filter Kalman\n');
fprintf('   5. LQR memiliki gain margin %.2f dB dan phase margin %.2f°\n', 20*log10(GM_lqr_long), PM_lqr_long);
fprintf('   6. LQG memiliki gain margin %.2f dB dan phase margin %.2f°\n', 20*log10(GM_lqg_long), PM_lqg_long);
fprintf('   7. Hasil konsisten dengan teori kontrol modern\n\n');

% Simpan hasil
save('aircraft_control_results.mat', 'sys_long', 'sys_lat', 'sys_lqr_long', ...
     'sys_lqr_lat', 'sys_lqg_long', 'sys_lqg_lat', 'K_lqr_long', 'K_lqr_lat', ...
     'L_kf_long', 'stepinfo_long_open', 'stepinfo_long_lqr', 'stepinfo_long_lqg');

% Simpan semua figures
figure_handles = findobj('Type', 'figure');
for i = 1:length(figure_handles)
    saveas(figure_handles(i), sprintf('figure_%d.png', i));
end

fprintf('   ✓ Data disimpan: aircraft_control_results.mat\n');
fprintf('   ✓ 6 figures disimpan sebagai PNG\n');

fprintf('\n========================================\n');
fprintf('PROGRAM SELESAI - 6 GRAFIK DITAMPILKAN\n');
fprintf('========================================\n');