clc; close all;
load('workspace_S1.mat');

fprintf('============================================================\n');
fprintf('  S2 - MIMO BEAMFORMING\n');
fprintf('============================================================\n\n');

steer = @(theta_deg, N) exp(1j * pi * (0:N-1)' * sind(theta_deg));

theta_user = 0;
w_comm     = steer(theta_user, Nt) / sqrt(Nt);

fprintf('  User direction      : %d degrees\n', theta_user);

W_sense = zeros(Nt, Ntargets);
for t = 1:Ntargets
    W_sense(:, t) = steer(target_angle(t), Nt) / sqrt(Nt);
    fprintf('  Target %d beam       : %+d degrees\n', t, target_angle(t));
end

W_isac = sqrt(beta) * (w_comm * ones(1, Ntargets)) + sqrt(1-beta) * W_sense;
W_isac = W_isac ./ max(vecnorm(W_isac), 1e-10);

w_sense_combined = sum(W_sense, 2);
w_sense_combined = w_sense_combined / max(norm(w_sense_combined), 1e-10);

w_isac_single = sqrt(beta) * w_comm + sqrt(1-beta) * w_sense_combined;
w_isac_single = w_isac_single / max(norm(w_isac_single), 1e-10);
w_tx          = w_isac_single;

theta_scan     = -90:0.5:90;
BP_comm        = zeros(size(theta_scan));
BP_isac        = zeros(size(theta_scan));
BP_isac_single = zeros(size(theta_scan));

for ii = 1:numel(theta_scan)
    sv = steer(theta_scan(ii), Nt);
    BP_comm(ii)        = abs(sv' * w_comm)^2;
    BP_isac(ii)        = sum(abs(sv' * W_isac).^2);
    BP_isac_single(ii) = abs(sv' * w_isac_single)^2;
end

BP_comm_dB        = 10*log10(BP_comm / max(BP_comm) + 1e-10);
BP_isac_dB        = 10*log10(BP_isac / max(BP_isac) + 1e-10);
BP_isac_single_dB = 10*log10(BP_isac_single / max(BP_isac_single) + 1e-10);

BF_gain_dB        = zeros(1, Ntargets);
BF_gain_single_dB = zeros(1, Ntargets);
for t = 1:Ntargets
    sv_t = steer(target_angle(t), Nt);
    BF_gain_dB(t)        = 20*log10(abs(sv_t' * W_isac(:, t)) + 1e-10);
    BF_gain_single_dB(t) = 20*log10(abs(sv_t' * w_isac_single) + 1e-10);
end
BF_gain_mean_lin = mean(10.^(BF_gain_single_dB/10));

fprintf('\n  Multi-beam BF gains : %s dB\n', num2str(round(BF_gain_dB, 2)));
fprintf('  Single-beam gains   : %s dB\n', num2str(round(BF_gain_single_dB, 2)));
fprintf('  Mean BF power gain  : %.3f dB\n', 10*log10(BF_gain_mean_lin));

figure('Name', 'S2 - ISAC Beampattern', 'Color', [0.12 0.12 0.15], ...
       'Position', [100 100 800 450]);

ax = axes;
set(ax, 'Color', [0.14 0.14 0.18], 'XColor', [0.7 0.7 0.8], ...
    'YColor', [0.7 0.7 0.8], 'GridColor', [0.35 0.35 0.45], 'GridAlpha', 0.5);
hold on; grid on; box on;

plot(theta_scan, BP_comm_dB, '--', 'Color', [1.0 0.6 0.2], 'LineWidth', 2, ...
    'DisplayName', 'Comm Only');
plot(theta_scan, BP_isac_dB, '-', 'Color', [0.3 0.8 1.0], 'LineWidth', 2, ...
    'DisplayName', 'ISAC Multi-beam');
plot(theta_scan, BP_isac_single_dB, '-', 'Color', [0.9 0.4 1.0], 'LineWidth', 1.8, ...
    'DisplayName', 'ISAC Single Beam');

for t = 1:Ntargets
    xline(target_angle(t), '--', 'Color', [0.4 1.0 0.5], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
end
xline(theta_user, ':', 'Color', [1.0 0.6 0.2], 'LineWidth', 1.5, 'HandleVisibility', 'off');

ylim([-50 3]); xlim([-90 90]);
title('ISAC ULA Beampattern', 'Color', [0.9 0.9 1.0], 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Angle (degrees)', 'Color', 'white');
ylabel('Normalized Array Gain (dB)', 'Color', 'white');
legend('Location', 'best', 'FontSize', 9);

save('workspace_S2.mat');
fprintf('\n  [SAVED] workspace_S2.mat\n');
fprintf('  Next -> run S3_Transmitter.m\n');
fprintf('============================================================\n');