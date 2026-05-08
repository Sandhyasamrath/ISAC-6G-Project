clc; close all;
load('workspace_S1.mat');
load('workspace_S2.mat');

fprintf('============================================================\n');
fprintf('  S7 - CRB ANALYSIS + DETECTION PROBABILITY\n');
fprintf('============================================================\n\n');

fprintf('  [PART A] Cramer-Rao Bound\n\n');

pilot_freqs = ((pilot_idx - mean(pilot_idx)) * delta_f).';
B_eff       = sqrt(mean(pilot_freqs.^2));
sym_times   = (0:Nsym-1) * Tsym;
T_eff       = sqrt(mean(sym_times.^2));

fprintf('  B_eff (RMS BW)      : %.4f MHz\n', B_eff/1e6);
fprintf('  T_eff (RMS time)    : %.6f ms\n', T_eff*1e3);
fprintf('  Pilot subcarriers   : %d\n', Npilot_sc);

CRB_range_vec  = zeros(size(SNR_dB_vec));
CRB_vel_vec    = zeros(size(SNR_dB_vec));
RMSE_range_vec = zeros(size(SNR_dB_vec));
RMSE_vel_vec   = zeros(size(SNR_dB_vec));

for ii = 1:numel(SNR_dB_vec)
    snr_r = 10^(SNR_dB_vec(ii)/10) * (1-beta) * mean(target_rcs.^2) * BF_gain_mean_lin;

    crb_tau           = 1 / (8 * pi^2 * max(snr_r, 1e-12) * B_eff^2);
    CRB_range_vec(ii) = sqrt(crb_tau) * c / 2;

    crb_fd            = 1 / (8 * pi^2 * max(snr_r, 1e-12) * T_eff^2);
    CRB_vel_vec(ii)   = sqrt(crb_fd) * lambda / 2;

    eff_factor         = 1 + 3 / (snr_r + 0.1);
    RMSE_range_vec(ii) = CRB_range_vec(ii) * eff_factor;
    RMSE_vel_vec(ii)   = CRB_vel_vec(ii) * eff_factor;
end

snr_r_op  = SNR_lin * (1-beta) * mean(target_rcs.^2) * BF_gain_mean_lin;
crb_r_op  = interp1(SNR_dB_vec, CRB_range_vec, SNR_dB, 'linear', 'extrap');
crb_v_op  = interp1(SNR_dB_vec, CRB_vel_vec, SNR_dB, 'linear', 'extrap');
rmse_r_op = interp1(SNR_dB_vec, RMSE_range_vec, SNR_dB, 'linear', 'extrap');

fprintf('\n  At operating SNR = %d dB:\n', SNR_dB);
fprintf('  Radar SNR (scaled)  : %.2f dB\n', 10*log10(snr_r_op));
fprintf('  CRB Range           : %.4f m\n', crb_r_op);
fprintf('  CRB Velocity        : %.4f m/s\n', crb_v_op);
fprintf('  Practical RMSE_R    : %.4f m\n', rmse_r_op);

fprintf('\n  [PART B] Detection Probability - Reference Curves\n\n');
fprintf('  Pfa                 : %.0e\n', Pfa);
fprintf('  alpha_cfar          : %.2f\n', alpha_cfar);

Pd_swerling0 = zeros(size(SNR_dB_vec));
Pd_swerling1 = zeros(size(SNR_dB_vec));
Pd_swerling3 = zeros(size(SNR_dB_vec));

for ii = 1:numel(SNR_dB_vec)
    snr_r = 10^(SNR_dB_vec(ii)/10) * mean(target_rcs.^2);

    arg              = sqrt(2*snr_r) - sqrt(-2*log(Pfa + 1e-15));
    Pd_swerling0(ii) = 0.5 * erfc(-arg / sqrt(2));
    Pd_swerling1(ii) = Pfa^(1 / (1 + snr_r / alpha_cfar));
    Pd_swerling3(ii) = (1 + 1/(snr_r + 1e-6)) * Pfa^(1/(1 + snr_r));

    Pd_swerling0(ii) = min(max(Pd_swerling0(ii), 0), 1);
    Pd_swerling1(ii) = min(max(Pd_swerling1(ii), 0), 1);
    Pd_swerling3(ii) = min(max(Pd_swerling3(ii), 0), 1);
end

figure('Name', 'S7 - CRB + Detection Probability', 'Color', [0.12 0.12 0.15], ...
       'Position', [50 50 1100 500]);

subplot(1, 2, 1);
ax1 = gca; styleAx(ax1);
semilogy(SNR_dB_vec, RMSE_range_vec, 'o-', 'Color', [0.3 0.8 1.0], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.3 0.8 1.0], 'MarkerSize', 6, 'DisplayName', 'RMSE Range');
hold on;
semilogy(SNR_dB_vec, CRB_range_vec, '--', 'Color', [0.3 0.8 1.0], 'LineWidth', 1.3, 'DisplayName', 'CRB Range');
semilogy(SNR_dB_vec, RMSE_vel_vec, 's-', 'Color', [1.0 0.5 0.2], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [1.0 0.5 0.2], 'MarkerSize', 6, 'DisplayName', 'RMSE Velocity');
semilogy(SNR_dB_vec, CRB_vel_vec, '--', 'Color', [1.0 0.5 0.2], 'LineWidth', 1.3, 'DisplayName', 'CRB Velocity');
xline(SNR_dB, '--', 'Color', 'white');
yline(delta_R_eff, '-.', 'Color', [0.4 1.0 0.5], 'Label', sprintf('\\DeltaR=%.1fm', delta_R_eff));
title('CRB vs Practical RMSE', 'Color', [0.9 0.9 1.0], 'FontSize', 10, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'Color', 'white');
ylabel('Estimation Error', 'Color', 'white');
legend('Location', 'best', 'FontSize', 7);

subplot(1, 2, 2);
ax2 = gca; styleAx(ax2);
plot(SNR_dB_vec, Pd_swerling0, 'o-', 'Color', [0.3 0.8 1.0], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.3 0.8 1.0], 'MarkerSize', 6, 'DisplayName', 'Swerling 0');
hold on;
plot(SNR_dB_vec, Pd_swerling1, 's-', 'Color', [1.0 0.5 0.2], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [1.0 0.5 0.2], 'MarkerSize', 6, 'DisplayName', 'Swerling I');
plot(SNR_dB_vec, Pd_swerling3, '^-', 'Color', [0.5 1.0 0.4], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.5 1.0 0.4], 'MarkerSize', 6, 'DisplayName', 'Swerling III');
yline(0.9, '--', 'Color', 'white');
xline(SNR_dB, '--', 'Color', 'white');
ylim([0 1.05]); xlim([min(SNR_dB_vec) max(SNR_dB_vec)]);
title('Detection Probability vs SNR', 'Color', [0.9 0.9 1.0], 'FontSize', 10, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'Color', 'white');
ylabel('P_d', 'Color', 'white');
legend('Location', 'best', 'FontSize', 8);

sgtitle('S7 - CRB & Detection Probability', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');

save('workspace_S7.mat', 'B_eff', 'T_eff', 'pilot_freqs', 'sym_times', ...
    'CRB_range_vec', 'CRB_vel_vec', 'RMSE_range_vec', 'RMSE_vel_vec', ...
    'Pd_swerling0', 'Pd_swerling1', 'Pd_swerling3', 'crb_r_op', 'crb_v_op', 'snr_r_op');

fprintf('\n  [SAVED] workspace_S7.mat\n');
fprintf('  Next -> run S8_ParetoSurface.m\n');
fprintf('============================================================\n');

function styleAx(ax)
set(ax, 'Color', [0.14 0.14 0.18], 'XColor', [0.7 0.7 0.8], ...
    'YColor', [0.7 0.7 0.8], 'GridColor', [0.35 0.35 0.45], 'GridAlpha', 0.5);
grid(ax, 'on'); box(ax, 'on'); hold(ax, 'on');
end