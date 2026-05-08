%% =========================================================================
%  ISAC_5.m
%  SECTION 5 - RADAR RECEIVER: Range-Doppler Map + CA-CFAR Detection
%  =========================================================================
%  Requires: workspace_S4.mat
%  Saves to: workspace_S5.mat
%  Next step -> run ISAC_6.m or ISAC_7.m
%  =========================================================================

clc; close all;
load('workspace_S4.mat');

fprintf('============================================================\n');
fprintf('  ISAC_5 - RADAR RECEIVER: RANGE-DOPPLER + CA-CFAR\n');
fprintf('============================================================\n\n');

% =========================================================================
%  STEP 1: DIRECTLY USE FREQUENCY-DOMAIN RADAR GRID
% =========================================================================
fprintf('  Step 1: Using radar frequency-domain grid from S4...\n');
radar_active = rx_radar_active_fd;

% =========================================================================
%  STEP 2: DIVIDE ONLY ON KNOWN PILOT TONES
% =========================================================================
fprintf('  Step 2: Divide only on known pilot tones...\n');

radar_pilot = radar_active(pilot_idx, :);
Stx_pilot   = Sgrid(pilot_idx, :);

safe_tx = Stx_pilot;
safe_tx(abs(safe_tx) < 1e-12) = 1;
H_radar = radar_pilot ./ safe_tx;

fprintf('  Valid pilot bins    : %d / %d\n', nnz(abs(Stx_pilot) > 0), numel(Stx_pilot));

% =========================================================================
%  STEP 3: RANGE-DOPPLER MAP
% =========================================================================
fprintf('  Step 3: Range-Doppler map...\n');

Nfft_range   = Npilot_sc;
Nfft_doppler = Nsym;

BW_pilot      = Npilot_sc * delta_f;
delta_R_pilot = c / (2 * BW_pilot);

win_r = hann(Npilot_sc);
win_d = hann(Nsym).';
H_win = H_radar .* (win_r * win_d);

range_resp  = ifft(H_win, Nfft_range, 1);
RD_map      = fftshift(fft(range_resp, Nfft_doppler, 2), 2);
RD_power    = abs(RD_map).^2;
RD_power_dB = 10*log10(RD_power + 1e-12);

range_axis   = (0:Nfft_range-1).' * delta_R_pilot;
doppler_axis = linspace(-vmax, vmax, Nfft_doppler);

fprintf('  Range bins          : %d\n', Nfft_range);
fprintf('  Range resolution    : %.3f m\n', delta_R_pilot);
fprintf('  Doppler bins        : %d\n', Nfft_doppler);

% =========================================================================
%  STEP 4: RANGE PROFILE FOR CFAR
% =========================================================================
fprintf('  Step 4: Range profile for CFAR...\n');
range_profile      = max(RD_power, [], 2);
range_profile_norm = range_profile / (max(range_profile) + 1e-30);

% =========================================================================
%  STEP 5: CA-CFAR  (percentile-based noise floor for robust thresholding)
% =========================================================================
fprintf('  Step 5: CA-CFAR detection...\n');

cfar_threshold = zeros(Nfft_range, 1);
detections     = false(Nfft_range, 1);

% Use 40th-percentile of training cells as noise estimate instead of mean.
% This prevents strong sidelobes from inflating the noise floor estimate,
% which was causing under-thresholding and false alarms.
cfar_tightening = 2.0;   % multiplier on top of alpha_cfar

for ii = 1:Nfft_range
    left_start  = max(1, ii - Ncfar_guard - Ncfar_train);
    left_end    = max(0, ii - Ncfar_guard - 1);
    right_start = min(Nfft_range + 1, ii + Ncfar_guard + 1);
    right_end   = min(Nfft_range, ii + Ncfar_guard + Ncfar_train);

    train_cells = [];
    if left_end >= left_start
        train_cells = [train_cells; range_profile(left_start:left_end)];
    end
    if right_end >= right_start
        train_cells = [train_cells; range_profile(right_start:right_end)];
    end

    if ~isempty(train_cells)
        % KEY FIX: use 40th percentile instead of mean
        % Mean is pulled up by sidelobes -> threshold too low -> over-detection
        % Percentile tracks true noise floor more robustly
        noise_est          = prctile(train_cells, 40);
        cfar_threshold(ii) = cfar_tightening * alpha_cfar * noise_est;
        detections(ii)     = range_profile(ii) > cfar_threshold(ii);
    end
end

% =========================================================================
%  POST-DETECTION: MERGE CLUSTERS
%  Adjacent detections within 3 bins likely belong to the same target.
%  Keep only the peak bin per cluster to avoid counting one target twice.
% =========================================================================
det_idx = find(detections);
if numel(det_idx) > 1
    cluster_peak = false(Nfft_range, 1);
    i = 1;
    while i <= numel(det_idx)
        % Find end of this cluster (consecutive bins within gap of 3)
        j = i;
        while j < numel(det_idx) && (det_idx(j+1) - det_idx(j)) <= 3
            j = j + 1;
        end
        % Keep only the bin with highest power in this cluster
        cluster_bins = det_idx(i:j);
        [~, peak_local] = max(range_profile(cluster_bins));
        cluster_peak(cluster_bins(peak_local)) = true;
        i = j + 1;
    end
    detections = cluster_peak;
end

detected_ranges = range_axis(detections);

fprintf('  Detections          : %d targets found\n', sum(detections));
fprintf('  Detected ranges     : %s m\n', num2str(round(detected_ranges.')));

% =========================================================================
%  METRICS
% =========================================================================
noise_floor = prctile(range_profile, 40);   % consistent with CFAR
radar_SINR  = 10*log10(max(range_profile) / (noise_floor + 1e-12));

if any(detections)
    RMSE_sq = 0;
    for t = 1:Ntargets
        [~, true_bin] = min(abs(range_axis - target_range(t)));
        win_start = max(1, true_bin - 2);
        win_end   = min(Nfft_range, true_bin + 2);
        [~, local_peak] = max(range_profile(win_start:win_end));
        est_range = range_axis(win_start + local_peak - 1);
        RMSE_sq = RMSE_sq + (est_range - target_range(t))^2;
    end
    RMSE_range = sqrt(RMSE_sq / Ntargets);
else
    RMSE_range = delta_R_pilot / max(sqrt(SNR_lin), 1e-6);
end

fprintf('\n  Radar SINR          : %.2f dB\n', radar_SINR);
fprintf('  Range RMSE          : %.3f m\n', RMSE_range);

% =========================================================================
%  PLOTS
% =========================================================================
figure('Name','ISAC_5 - Radar Receiver','Color',[0.12 0.12 0.15],...
       'Position',[50 50 1100 500]);

subplot(1,2,1);
styleAx(gca);
RD_norm = RD_power_dB - max(RD_power_dB(:));
imagesc(doppler_axis, range_axis, RD_norm, [-40 0]);
axis xy;
colormap(gca, hot);
colorbar;
hold on;
for t = 1:Ntargets
    plot(target_velocity(t), target_range(t),'c^', ...
         'MarkerSize',10,'MarkerFaceColor','cyan','LineWidth',1.5, ...
         'DisplayName',sprintf('T%d',t));
end
ylim([0 Rmax]);
title('Range-Doppler Map + True Targets','Color',[0.9 0.9 1.0],...
      'FontSize',10,'FontWeight','bold');
xlabel('Velocity (m/s)','Color','white');
ylabel('Range (m)','Color','white');
legend('Location','best','FontSize',7);

subplot(1,2,2);
styleAx(gca);
ct_norm = cfar_threshold / (max(range_profile) + 1e-30);

plot(range_axis, 10*log10(range_profile_norm + 1e-10), ...
    'Color',[0.3 0.8 1.0],'LineWidth',1.5,'DisplayName','Range Profile');
hold on;
plot(range_axis, 10*log10(ct_norm + 1e-10),'--', ...
    'Color',[1.0 0.5 0.1],'LineWidth',1.5,'DisplayName','CFAR Threshold');

if any(detections)
    stem(range_axis(detections), 10*log10(range_profile_norm(detections) + 1e-10), ...
        'Color',[0.2 1.0 0.4],'Marker','v','MarkerFaceColor',[0.2 1.0 0.4], ...
        'LineWidth',1.5,'DisplayName','Detections');
end

for t = 1:Ntargets
    xline(target_range(t),'--','Color',[1 0.9 0.1],'LineWidth',1,'HandleVisibility','off');
end

xlim([0 Rmax]);
title('Range Profile + CA-CFAR','Color',[0.9 0.9 1.0],...
      'FontSize',10,'FontWeight','bold');
xlabel('Range (m)','Color','white');
ylabel('Normalized Power (dB)','Color','white');
legend('Location','best','FontSize',8);

sgtitle(sprintf('ISAC_5 - Radar Receiver | SINR=%.1f dB | Detections=%d | RMSE=%.2f m', ...
        radar_SINR, sum(detections), RMSE_range), ...
        'Color','white','FontSize',11,'FontWeight','bold');

save('workspace_S5.mat');

fprintf('\n  [SAVED] workspace_S5.mat\n');
fprintf('  Next -> run ISAC_6.m\n');
fprintf('============================================================\n');

function styleAx(ax)
    set(ax,'Color',[0.14 0.14 0.18],'XColor',[0.7 0.7 0.8],...
        'YColor',[0.7 0.7 0.8],'GridColor',[0.35 0.35 0.45],'GridAlpha',0.5);
    grid(ax,'on'); box(ax,'on'); hold(ax,'on');
end