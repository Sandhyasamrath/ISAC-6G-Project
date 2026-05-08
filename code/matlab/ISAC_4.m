%% =========================================================================
%  ISAC_4.m
%  SECTION 4 - CHANNEL MODELING (Comm + Radar Dual Path)
%  =========================================================================
%  Requires: workspace_S3.mat
%  Saves to: workspace_S4.mat
%  Next step -> run ISAC_5.m and/or ISAC_6.m
%  =========================================================================

clc; close all;
load('workspace_S3.mat');

fprintf('============================================================\n');
fprintf('  ISAC_4 - CHANNEL MODELING\n');
fprintf('============================================================\n\n');

% =========================================================================
%  PATH A: COMMUNICATION CHANNEL - Frequency-Selective Rayleigh Fading
% =========================================================================
fprintf('  [PATH A] Communication Channel - Rayleigh Fading\n');

Lch       = 6;
tap_power = [0, -3, -6, -9, -12, -15];
tap_amp   = 10.^(tap_power/20);
tap_coeff = tap_amp .* (randn(1, Lch) + 1j*randn(1, Lch)) / sqrt(2);

fprintf('  Taps                : %d\n', Lch);
fprintf('  PDP (dB)            : %s\n', num2str(tap_power));
fprintf('  CP length           : %d  (>= Lch: %s)\n', Ncp, mat2str(Ncp >= Lch));

rx_comm_td = zeros(Nsc + Ncp + Lch - 1, Nsym);
for s = 1:Nsym
    rx_comm_td(:, s) = conv(tx_td(:, s), tap_coeff);
end
rx_comm_td = rx_comm_td(1:Nsc+Ncp, :);

% Scale AWGN to the actual received signal power
sig_pow_comm   = mean(abs(rx_comm_td(:)).^2);
noise_var_comm = sig_pow_comm / SNR_lin;

noise_comm = sqrt(noise_var_comm/2) .* ...
    (randn(size(rx_comm_td)) + 1j*randn(size(rx_comm_td)));
rx_comm_td = rx_comm_td + noise_comm;

fprintf('  Comm SNR            : %d dB\n', SNR_dB);
fprintf('  Comm signal power   : %.4e\n', sig_pow_comm);
fprintf('  Comm noise variance : %.4e\n', noise_var_comm);
fprintf('  rx_comm_td size     : [%d x %d]\n\n', size(rx_comm_td,1), size(rx_comm_td,2));

% =========================================================================
%  PATH B: RADAR CHANNEL - Frequency-Domain Echo Model
% =========================================================================
fprintf('  [PATH B] Radar Channel - Frequency-Domain Echo Model\n');

k_active  = (-Nactive/2 : Nactive/2-1).';
f_vec     = k_active * delta_f;
t_sym_vec = (0:Nsym-1) * Tsym;

rx_radar_active_fd = zeros(Nactive, Nsym);

for t = 1:Ntargets
    tau_t = 2 * target_range(t) / c;
    fd_t  = 2 * target_velocity(t) / lambda;

    range_shift   = exp(-1j * 2*pi * f_vec * tau_t);
    doppler_shift = exp( 1j * 2*pi * fd_t * t_sym_vec);

    gain_lin = target_rcs(t) * 10^(BF_gain_single_dB(t)/20);

    target_echo = gain_lin * ((range_shift * doppler_shift) .* Sgrid);
    rx_radar_active_fd = rx_radar_active_fd + target_echo;

    fprintf('  Target %d  R=%.0fm  v=%+.0fm/s  tau=%.3fus  fd=%.1fHz  BF=%.1fdB\n', ...
        t, target_range(t), target_velocity(t), tau_t*1e6, fd_t, BF_gain_single_dB(t));
end

noise_radar = sqrt(noise_var/2) .* ...
    (randn(size(rx_radar_active_fd)) + 1j*randn(size(rx_radar_active_fd)));
rx_radar_active_fd = rx_radar_active_fd + noise_radar;

rx_radar_fd = zeros(Nsc, Nsym);
rx_radar_fd(guard_left+1 : guard_left+Nactive, :) = rx_radar_active_fd;

rx_radar_td = zeros(Nsc + Ncp, Nsym);
for s = 1:Nsym
    full_sym_fd = zeros(Nsc,1);
    full_sym_fd(guard_left+1 : guard_left+Nactive) = rx_radar_active_fd(:,s);
    sym_td = ifft(ifftshift(full_sym_fd), Nsc);
    rx_radar_td(:,s) = [sym_td(end-Ncp+1:end); sym_td];
end

fprintf('\n  rx_radar_active_fd size : [%d x %d]\n', size(rx_radar_active_fd,1), size(rx_radar_active_fd,2));
fprintf('  rx_radar_td size        : [%d x %d]\n', size(rx_radar_td,1), size(rx_radar_td,2));

% =========================================================================
%  PLOT
% =========================================================================
figure('Name','ISAC_4 - Channel Modeling','Color',[0.12 0.12 0.15],...
       'Position',[100 100 1000 420]);

subplot(1,2,1);
styleAx(gca);
plot(real(rx_comm_td(:,1)),'Color',[0.3 0.8 1.0],'LineWidth',0.8,'DisplayName','Re{rx comm}');
hold on;
plot(real(rx_radar_td(:,1)),'Color',[1.0 0.5 0.2],'LineWidth',0.8,'DisplayName','Re{rx radar}');
title('Received Signals - Symbol 1','Color',[0.9 0.9 1.0],'FontSize',10,'FontWeight','bold');
xlabel('Sample Index','Color','white');
ylabel('Amplitude','Color','white');
legend('Location','best','FontSize',8);

subplot(1,2,2);
styleAx(gca);
psd_comm  = abs(fftshift(fft(rx_comm_td(:,1), 1024))).^2;
psd_radar = abs(fftshift(fft(rx_radar_td(:,1), 1024))).^2;
f_axis = linspace(-BW/2, BW/2, 1024) / 1e6;
plot(f_axis, 10*log10(psd_comm  + 1e-10),'Color',[0.3 0.8 1.0],'LineWidth',1,'DisplayName','Comm Channel');
hold on;
plot(f_axis, 10*log10(psd_radar + 1e-10),'Color',[1.0 0.5 0.2],'LineWidth',1,'DisplayName','Radar Echo');
title('Power Spectral Density','Color',[0.9 0.9 1.0],'FontSize',10,'FontWeight','bold');
xlabel('Frequency (MHz)','Color','white');
ylabel('PSD (dB)','Color','white');
legend('Location','best','FontSize',8);

sgtitle('ISAC_4 - Dual-Path Channel Modeling','Color','white','FontSize',12,'FontWeight','bold');

save('workspace_S4.mat');

fprintf('\n  [SAVED] workspace_S4.mat\n');
fprintf('  Next -> run ISAC_5.m and/or ISAC_6.m\n');
fprintf('============================================================\n');

function styleAx(ax)
    set(ax,'Color',[0.14 0.14 0.18],'XColor',[0.7 0.7 0.8],...
        'YColor',[0.7 0.7 0.8],'GridColor',[0.35 0.35 0.45],'GridAlpha',0.5);
    grid(ax,'on'); box(ax,'on'); hold(ax,'on');
end
