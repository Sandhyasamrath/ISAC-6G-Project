clc; close all;

load('workspace_S1.mat');
load('workspace_S2.mat');
load('workspace_S3.mat');
load('workspace_S4.mat');
load('workspace_S5.mat');
load('workspace_S6.mat');
load('workspace_S7.mat');
load('workspace_S8.mat');

fprintf('============================================================\n');
fprintf('  S9 - FINAL COMBINED PLOTS\n');
fprintf('============================================================\n\n');

figure('Name', 'ISAC 6G - Final Results', 'NumberTitle', 'off', ...
       'Position', [30 30 1600 950], 'Color', [0.10 0.10 0.13]);

sgtitle({'Optimal Joint ISAC Waveform Design - 6G Networks', ...
    sprintf('\\alpha^* = %.2f | \\beta^* = %.2f | \\lambda^* = %.2f | SNR = %d dB | f_c = 28 GHz | %d\\times%d MIMO', ...
    opt_alpha, beta, opt_lambda, SNR_dB, Nt, Nr)}, 'Color', 'white', 'FontSize', 11, 'FontWeight', 'bold');

tcol = [0.90 0.90 1.00];
axBG = [0.14 0.14 0.18];
gcol = [0.35 0.35 0.45];

ax1 = subplot(2, 4, 1); 
styleAx(ax1, axBG, gcol);
plot(theta_scan, BP_comm_dB, '--', 'Color', [1.0 0.6 0.2], 'LineWidth', 2, 'DisplayName', 'Comm');
hold on;
plot(theta_scan, BP_isac_single_dB, '-', 'Color', [0.9 0.4 1.0], 'LineWidth', 2, 'DisplayName', 'ISAC Single');
for t = 1:Ntargets
    xline(target_angle(t), '--', 'Color', [0.4 1.0 0.5], 'LineWidth', 1);
end
xline(0, ':', 'Color', [1.0 0.6 0.2], 'LineWidth', 1.2);
ylim([-50 3]); xlim([-90 90]);
title('(1) ISAC Beampattern', 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold');
xlabel('Angle (deg)', 'Color', 'white');
ylabel('Gain (dB)', 'Color', 'white');
legend('Location', 'best', 'FontSize', 6);

ax2 = subplot(2, 4, 2); 
styleAx(ax2, axBG, gcol);
RD_norm = RD_power_dB - max(RD_power_dB(:));
imagesc(doppler_axis, range_axis, RD_norm, [-40 0]); 
axis xy;
colormap(ax2, hot);
colorbar;
hold on;
for t = 1:Ntargets
    plot(target_velocity(t), target_range(t), 'c^', 'MarkerSize', 9, ...
         'MarkerFaceColor', 'cyan', 'LineWidth', 1.5);
end
ylim([0 Rmax]);
title('(2) Range-Doppler Map', 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold');
xlabel('Velocity (m/s)', 'Color', 'white');
ylabel('Range (m)', 'Color', 'white');

ax3 = subplot(2, 4, 3); 
styleAx(ax3, axBG, gcol);
semilogy(SNR_dB_vec, RMSE_range_vec, 'o-', 'Color', [0.3 0.8 1.0], 'LineWidth', 1.5, ...
    'MarkerFaceColor', [0.3 0.8 1.0], 'MarkerSize', 5, 'DisplayName', 'RMSE Range');
hold on;
semilogy(SNR_dB_vec, CRB_range_vec, '--', 'Color', [0.3 0.8 1.0], 'LineWidth', 1.2, 'DisplayName', 'CRB Range');
semilogy(SNR_dB_vec, RMSE_vel_vec, 's-', 'Color', [1.0 0.5 0.2], 'LineWidth', 1.5, ...
    'MarkerFaceColor', [1.0 0.5 0.2], 'MarkerSize', 5, 'DisplayName', 'RMSE Velocity');
semilogy(SNR_dB_vec, CRB_vel_vec, '--', 'Color', [1.0 0.5 0.2], 'LineWidth', 1.2, 'DisplayName', 'CRB Velocity');
xline(SNR_dB, '--', 'Color', 'white');
title('(3) CRB vs RMSE', 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'Color', 'white');
ylabel('Error', 'Color', 'white');
legend('Location', 'best', 'FontSize', 6);

ax4 = subplot(2, 4, 4); 
styleAx(ax4, axBG, gcol);
plot(SNR_dB_vec, Pd_swerling0, 'o-', 'Color', [0.3 0.8 1.0], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.3 0.8 1.0], 'MarkerSize', 5, 'DisplayName', 'Swerling 0');
hold on;
plot(SNR_dB_vec, Pd_swerling1, 's-', 'Color', [1.0 0.5 0.2], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [1.0 0.5 0.2], 'MarkerSize', 5, 'DisplayName', 'Swerling I');
plot(SNR_dB_vec, Pd_swerling3, '^-', 'Color', [0.5 1.0 0.4], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.5 1.0 0.4], 'MarkerSize', 5, 'DisplayName', 'Swerling III');
yline(0.9, '--', 'Color', 'white');
xline(SNR_dB, '--', 'Color', 'white');
ylim([0 1.05]); 
xlim([min(SNR_dB_vec) max(SNR_dB_vec)]);
title('(4) Detection Probability', 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'Color', 'white');
ylabel('P_d', 'Color', 'white');
legend('Location', 'best', 'FontSize', 6);

ax5 = subplot(2, 4, 5); 
styleAx(ax5, axBG, gcol);
semilogy(SNR_dB_vec, BER_ref, '--', 'Color', [1.0 0.6 0.2], 'LineWidth', 1.5, 'DisplayName', 'AWGN Ref');
hold on;
semilogy(SNR_dB, BER_sim, 'p', 'Color', [0.2 1.0 0.4], 'MarkerSize', 14, ...
    'MarkerFaceColor', [0.2 1.0 0.4], 'DisplayName', sprintf('Op.Pt (%d dB)', SNR_dB));
yline(1e-3, '--', 'Color', 'white');
ylim([1e-6 1]); 
xlim([min(SNR_dB_vec) max(SNR_dB_vec)]);
title('(5) BER', 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'Color', 'white');
ylabel('BER', 'Color', 'white');
legend('Location', 'best', 'FontSize', 6);

ax6 = subplot(2, 4, 6);
[A_grid, L_grid] = meshgrid(alpha_sweep, lambda_sweep);
surf(A_grid, L_grid, U_surface, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
colormap(ax6, parula);
colorbar;
hold on;
plot3(opt_alpha, opt_lambda, U_max + 0.01, 'r*', 'MarkerSize', 13, 'LineWidth', 2);
plot3(alpha_sweep, 0.5 * ones(1, Na), U_surface(lw_half_idx, :), 'w-', 'LineWidth', 2.5);
plot3(knee_alpha, 0.5, knee_U, 'wp', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.2 0.2]);
set(ax6, 'Color', axBG, 'XColor', gcol, 'YColor', gcol, 'ZColor', gcol, 'FontSize', 7);
xlabel('\alpha', 'Color', 'white');
ylabel('\lambda', 'Color', 'white');
zlabel('Utility', 'Color', 'white');
title('(6) Pareto Surface', 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold');
view(40, 28);

ax7 = subplot(2, 4, 7); 
styleAx(ax7, axBG, gcol);
colors7 = {[0.3 0.8 1.0], [0.5 1.0 0.4], [1.0 0.5 0.2]};
for li_p = 1:numel(lw_plot)
    plot(se_slices(li_p, :), qs_slices(li_p, :), 'o-', 'Color', colors7{li_p}, ...
        'LineWidth', 1.5, 'MarkerFaceColor', colors7{li_p}, 'MarkerSize', 4, ...
        'DisplayName', sprintf('\\lambda=%.1f', lw_plot(li_p)));
end
[~, half_idx] = min(abs(lw_plot - 0.5));
plot(se_slices(half_idx, knee_idx), qs_slices(half_idx, knee_idx), 'p', 'MarkerSize', 13, ...
    'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', 'white', ...
    'DisplayName', sprintf('Knee \\alpha^*=%.2f', knee_alpha));
title('(7) Pareto Slices', 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold');
xlabel('SE (bps/Hz)', 'Color', 'white');
ylabel('Q_{sense}', 'Color', 'white');
legend('Location', 'best', 'FontSize', 6);

ax8 = subplot(2, 4, 8); 
styleAx(ax8, axBG, gcol);
rp_norm = range_profile / (max(range_profile) + 1e-30);
ct_norm = cfar_threshold / (max(range_profile) + 1e-30);
plot(range_axis, 10*log10(rp_norm + 1e-10), 'Color', [0.3 0.8 1.0], ...
    'LineWidth', 1.2, 'DisplayName', 'Range Profile');
hold on;
plot(range_axis, 10*log10(ct_norm + 1e-10), '--', 'Color', [1.0 0.5 0.1], ...
    'LineWidth', 1.5, 'DisplayName', 'CFAR Threshold');
if any(detections)
    stem(range_axis(detections), 10*log10(rp_norm(detections) + 1e-10), 'Color', [0.2 1.0 0.4], ...
        'Marker', 'v', 'MarkerFaceColor', [0.2 1.0 0.4], 'LineWidth', 1.5, 'DisplayName', 'Detections');
end
for t = 1:Ntargets
    xline(target_range(t), '--', 'Color', [1 0.9 0.1], 'LineWidth', 1);
end
xlim([0 Rmax]);
title('(8) Range Profile + CFAR', 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold');
xlabel('Range (m)', 'Color', 'white');
ylabel('Power (dB)', 'Color', 'white');
legend('Location', 'best', 'FontSize', 6);

fprintf('============================================================\n');
fprintf('  FINAL SUMMARY\n');
fprintf('============================================================\n');
fprintf('  Beamforming gains    : %s dB\n', num2str(round(BF_gain_single_dB, 1)));
fprintf('  Radar SINR           : %.2f dB\n', radar_SINR);
fprintf('  Range RMSE           : %.3f m\n', RMSE_range);
fprintf('  CFAR detections      : %d\n', sum(detections));
fprintf('  BER (sim)            : %.2e\n', BER_sim);
fprintf('  SE                   : %.2f bps/Hz\n', SE_op);
fprintf('  CRB Range            : %.4f m\n', crb_r_op);
fprintf('  CRB Velocity         : %.4f m/s\n', crb_v_op);
fprintf('  Pareto alpha*        : %.2f\n', opt_alpha);
fprintf('  Pareto lambda*       : %.2f\n', opt_lambda);
fprintf('============================================================\n');

function styleAx(ax, bgColor, gColor)
    set(ax, 'Color', bgColor, 'XColor', gColor, 'YColor', gColor, ...
        'GridColor', gColor, 'GridAlpha', 0.5, 'FontSize', 7.5, 'FontName', 'Helvetica');
    grid(ax, 'on');
    box(ax, 'on');
    hold(ax, 'on');
end
