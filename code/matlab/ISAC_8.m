clc; close all;
load('workspace_S1.mat');
load('workspace_S2.mat');
load('workspace_S7.mat');

fprintf('============================================================\n');
fprintf('  S8 - FULL PARETO SURFACE U(lambda, alpha)\n');
fprintf('============================================================\n\n');

alpha_sweep  = 0.05:0.05:0.95;
lambda_sweep = 0.00:0.05:1.00;
Na           = numel(alpha_sweep);
Nl           = numel(lambda_sweep);

U_surface  = zeros(Nl, Na);
SE_surface = zeros(Nl, Na);
Qs_surface = zeros(Nl, Na);

for li = 1:Nl
    lw = lambda_sweep(li);
    for ai = 1:Na
        a   = alpha_sweep(ai);
        Npi = max(round(a * Nactive), 8);
        Npi = min(Npi, Nactive - 8);
        Ndi = Nactive - Npi;

        pilot_start_a = floor((Nactive - Npi) / 2) + 1;
        pilot_idx_a   = pilot_start_a:(pilot_start_a + Npi - 1);
        pilot_freqs_a = ((pilot_idx_a - mean(pilot_idx_a)) * delta_f).';
        B_eff_a       = sqrt(mean(pilot_freqs_a.^2));

        SE_a    = (Ndi / Nactive) * log2(1 + beta * SNR_lin);
        SE_norm = SE_a / log2(1 + SNR_lin);

        snr_r   = 10^(SNR_dB/10) * (1-beta) * mean(target_rcs.^2) * BF_gain_mean_lin * (Npi / Nactive);
        crb_r   = 1 / (8*pi^2 * max(snr_r, 1e-12) * B_eff_a^2);
        Q_sense = 1 / (1 + sqrt(crb_r) * c/2);

        U_surface(li, ai)  = lw * SE_norm + (1-lw) * Q_sense;
        SE_surface(li, ai) = SE_a;
        Qs_surface(li, ai) = Q_sense;
    end
end

[U_max, U_max_idx] = max(U_surface(:));
[li_opt, ai_opt]   = ind2sub([Nl, Na], U_max_idx);
opt_lambda         = lambda_sweep(li_opt);
opt_alpha          = alpha_sweep(ai_opt);

[~, lw_half_idx] = min(abs(lambda_sweep - 0.5));
U_half           = U_surface(lw_half_idx, :);
[~, knee_idx]    = max(U_half);
knee_alpha       = alpha_sweep(knee_idx);
knee_U           = U_half(knee_idx);

lw_plot   = [0.2, 0.5, 0.8];
se_slices = zeros(numel(lw_plot), Na);
qs_slices = zeros(numel(lw_plot), Na);
for li_p = 1:numel(lw_plot)
    [~, idx] = min(abs(lambda_sweep - lw_plot(li_p)));
    se_slices(li_p, :) = SE_surface(idx, :);
    qs_slices(li_p, :) = Qs_surface(idx, :);
end

fprintf('  Global optimal:\n');
fprintf('    lambda* = %.2f\n', opt_lambda);
fprintf('    alpha*  = %.2f\n', opt_alpha);
fprintf('    U_max   = %.4f\n', U_max);
fprintf('\n  Knee point (lambda=0.5):\n');
fprintf('    alpha*  = %.2f\n', knee_alpha);
fprintf('    U_knee  = %.4f\n', knee_U);

figure('Name', 'S8 - Pareto Surface', 'Color', [0.12 0.12 0.15], ...
       'Position', [50 50 1200 520]);

subplot(1, 2, 1);
[A_grid, L_grid] = meshgrid(alpha_sweep, lambda_sweep);
surf(A_grid, L_grid, U_surface, 'EdgeColor', 'none', 'FaceAlpha', 0.88);
colormap(parula);
colorbar;
hold on;
plot3(opt_alpha, opt_lambda, U_max + 0.01, 'r*', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', 'Global Optimum');
plot3(alpha_sweep, 0.5 * ones(1, Na), U_surface(lw_half_idx, :), 'w-', 'LineWidth', 2.5, 'DisplayName', '\lambda=0.5 Slice');
plot3(knee_alpha, 0.5, knee_U, 'wp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.2 0.2], 'DisplayName', 'Knee Point');
set(gca, 'Color', [0.14 0.14 0.18], 'XColor', [0.7 0.7 0.8], 'YColor', [0.7 0.7 0.8], 'ZColor', [0.7 0.7 0.8], 'FontSize', 8);
xlabel('\alpha (Pilot Ratio)', 'Color', 'white');
ylabel('\lambda (ISAC Weight)', 'Color', 'white');
zlabel('Total Utility', 'Color', 'white');
title('3D Pareto Surface U(\lambda, \alpha)', 'Color', [0.9 0.9 1.0], 'FontSize', 10, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 7);
view(40, 28);

subplot(1, 2, 2);
ax2 = gca;
set(ax2, 'Color', [0.14 0.14 0.18], 'XColor', [0.7 0.7 0.8], ...
    'YColor', [0.7 0.7 0.8], 'GridColor', [0.35 0.35 0.45], 'GridAlpha', 0.5);
grid on; box on; hold on;

colors = {[0.3 0.8 1.0], [0.5 1.0 0.4], [1.0 0.5 0.2]};
for li_p = 1:numel(lw_plot)
    plot(se_slices(li_p, :), qs_slices(li_p, :), 'o-', 'Color', colors{li_p}, 'LineWidth', 1.8, ...
        'MarkerFaceColor', colors{li_p}, 'MarkerSize', 5, 'DisplayName', sprintf('\\lambda = %.1f', lw_plot(li_p)));
end

[~, half_idx] = min(abs(lw_plot - 0.5));
plot(se_slices(half_idx, knee_idx), qs_slices(half_idx, knee_idx), 'p', 'MarkerSize', 14, ...
    'MarkerFaceColor', [1.0 0.2 0.2], 'MarkerEdgeColor', 'white', ...
    'DisplayName', sprintf('Knee (\\alpha^*=%.2f)', knee_alpha));

title('Pareto Slices: SE vs Sensing Quality', 'Color', [0.9 0.9 1.0], 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Spectral Efficiency (bps/Hz)', 'Color', 'white');
ylabel('Sensing Quality Q_{sense}', 'Color', 'white');
legend('Location', 'best', 'FontSize', 8);

sgtitle(sprintf('S8 - Pareto Surface | \\alpha^*=%.2f | \\lambda^*=%.2f | U_{max}=%.3f', ...
    opt_alpha, opt_lambda, U_max), 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');

save('workspace_S8.mat', 'alpha_sweep', 'lambda_sweep', 'U_surface', ...
    'SE_surface', 'Qs_surface', 'U_max', 'opt_lambda', 'opt_alpha', ...
    'knee_alpha', 'knee_U', 'knee_idx', 'lw_half_idx', 'lw_plot', ...
    'se_slices', 'qs_slices');

fprintf('\n  [SAVED] workspace_S8.mat\n');
fprintf('  Next -> run S9_FinalPlots.m\n');
fprintf('============================================================\n');