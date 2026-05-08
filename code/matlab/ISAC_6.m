%% =========================================================================
%  ISAC_6.m
%  SECTION 6 - COMMUNICATION RECEIVER: LS Estimation + MMSE + BER
%  HIGH-FIDELITY EDITION WITH EVM + PHASE CORRECTION
%  =========================================================================
%  Requires: workspace_S4.mat
%  Saves to: workspace_S6.mat
%  Next step -> run ISAC_7.m
%  =========================================================================

clc; close all;
load('workspace_S4.mat');

fprintf('============================================================\n');
fprintf('  ISAC_6 - COMMUNICATION RECEIVER: LS + MMSE + BER\n');
fprintf('============================================================\n\n');

% =========================================================================
%  STEP 1: REMOVE CP + FFT
% =========================================================================
fprintf('  Step 1: CP removal + FFT on comm signal...\n');

rx_comm_noCP = rx_comm_td(Ncp+1:Ncp+Nsc, :);
Y_full       = fftshift(fft(rx_comm_noCP, Nsc, 1), 1);
Y_active     = Y_full(guard_left+1 : guard_left+Nactive, :);

% =========================================================================
%  STEP 2: FULL-TRAINING CHANNEL ESTIMATION
%
%  IMPORTANT:
%    The first Npilot_sym symbols in ISAC_3 are full training symbols on
%    ALL active subcarriers, with per-subcarrier power = Ptrain_full.
% =========================================================================
fprintf('  Step 2: Full-training LS channel estimation...\n');

Y_pilot        = Y_active(:, 1:Npilot_sym);
S_pilot_known  = sqrt(Ptrain_full) * ones(Nactive, Npilot_sym);

H_est_raw = Y_pilot ./ S_pilot_known;
H_est_avg = mean(H_est_raw, 2);

% Full-grid DFT denoising
H_ls_full = zeros(Nsc, 1);
H_ls_full(guard_left+1 : guard_left+Nactive) = H_est_avg;

h_raw_full = ifft(ifftshift(H_ls_full), Nsc);

h_win_full = zeros(Nsc, 1);
h_win_full(1:Lch) = h_raw_full(1:Lch);

H_est_full = fftshift(fft(h_win_full, Nsc));
H_est      = H_est_full(guard_left+1 : guard_left+Nactive);

fprintf('  Raw LS range        : [%.1f, %.1f] dB\n', ...
    min(20*log10(abs(H_est_avg)+1e-10)), max(20*log10(abs(H_est_avg)+1e-10)));
fprintf('  DFT-clean range     : [%.1f, %.1f] dB\n', ...
    min(20*log10(abs(H_est)+1e-10)), max(20*log10(abs(H_est)+1e-10)));

% =========================================================================
%  DEBUG CHECK: TRUE CHANNEL VS ESTIMATED CHANNEL
% =========================================================================
h_true_full = zeros(Nsc, 1);
h_true_full(1:Lch) = tap_coeff(:);
H_true_full = fftshift(fft(h_true_full, Nsc));
H_true      = H_true_full(guard_left+1 : guard_left+Nactive);

chan_nmse = norm(H_est - H_true)^2 / norm(H_true)^2;
fprintf('  Channel NMSE        : %.4e\n', chan_nmse);

% =========================================================================
%  STEP 3: DATA EXTRACTION + MMSE EQUALIZATION
% =========================================================================
fprintf('  Step 3: MMSE equalization on data subcarriers...\n');

Y_data      = Y_active(:, Npilot_sym+1:end);
Y_data_only = Y_data(data_idx, :);
H_data_only = H_est(data_idx);

% FFT-domain noise variance
noise_var_fd = Nsc * noise_var_comm;

W_mmse = conj(H_data_only) ./ (abs(H_data_only).^2 + noise_var_fd / max(Pdata, 1e-12));
X_hat  = Y_data_only .* W_mmse;

% =========================================================================
%  STEP 4: OPTIONAL COMMON PHASE CORRECTION
%
%  Use transmitted data symbols as reference to estimate residual constant
%  phase rotation. This is useful as a debug/calibration step.
% =========================================================================
fprintf('  Step 4: Residual phase correction + BER/EVM...\n');

ideal_symbols = reshape(tx_symbols, Ndata_sc, Ndata_sym);

phase_ratio = X_hat(:) ./ ideal_symbols(:);
mean_phase_err = angle(mean(phase_ratio));
X_hat_corr = X_hat * exp(-1j * mean_phase_err);

fprintf('  Mean phase error    : %.4f rad\n', mean_phase_err);

% =========================================================================
%  STEP 5: EVM + BER
% =========================================================================
errors  = X_hat_corr - ideal_symbols;
evm_val = sqrt(mean(abs(errors(:)).^2) / mean(abs(ideal_symbols(:)).^2)) * 100;

rx_sym_vec  = X_hat_corr(:);
rx_bits_est = qamdemod(rx_sym_vec, M, 'OutputType', 'bit', 'UnitAveragePower', true);
[num_err, BER_sim] = biterr(tx_bits, rx_bits_est);

Eb_N0_lin  = SNR_lin / k_bits;
BER_theory = (3/(2*k_bits)) * erfc(sqrt(0.4 * Eb_N0_lin));

SE_op = (Ndata_sc / Nactive) * log2(1 + beta * SNR_lin);

BER_ref   = zeros(size(SNR_dB_vec));
BER_sweep = zeros(size(SNR_dB_vec));
SE_sweep  = zeros(size(SNR_dB_vec));

for ii = 1:length(SNR_dB_vec)
    snr_i         = 10^(SNR_dB_vec(ii)/10);
    BER_ref(ii)   = max((3/(2*k_bits)) * erfc(sqrt(0.4 * (snr_i / k_bits))), 1e-6);
    BER_sweep(ii) = BER_ref(ii);
    SE_sweep(ii)  = (Ndata_sc / Nactive) * log2(1 + beta * snr_i);
end

fprintf('\n  --- Communication Performance ---\n');
fprintf('  EVM (%%)             : %.2f\n', evm_val);
fprintf('  BER (simulated)     : %.3e\n', BER_sim);
fprintf('  BER (AWGN ref)      : %.3e\n', BER_theory);
fprintf('  Bit errors          : %d / %d\n', num_err, Nbits);
fprintf('  Spectral efficiency : %.2f bps/Hz\n', SE_op);

% =========================================================================
%  PLOTS
% =========================================================================
figure('Name','ISAC_6 - Communication Receiver','Color',[0.12 0.12 0.15],...
       'Position',[50 50 1450 450]);

subplot(1,4,1);
styleAx(gca);
scatter(real(rx_sym_vec), imag(rx_sym_vec), 5,...
        [0.3 0.7 1.0],'filled','DisplayName','Received');
hold on;
ideal_const = qammod((0:M-1)', M, 'UnitAveragePower', true);
scatter(real(ideal_const), imag(ideal_const), 100,...
        [1.0 0.4 0.3],'filled','MarkerEdgeColor','white',...
        'LineWidth',0.5,'DisplayName','Ideal');
title(sprintf('Equalized 16-QAM (EVM = %.1f%%)', evm_val),...
      'Color',[0.9 0.9 1.0],'FontSize',10,'FontWeight','bold');
xlabel('In-Phase','Color','white');
ylabel('Quadrature','Color','white');
legend('Location','best','FontSize',8);
axis square;

subplot(1,4,2);
styleAx(gca);
plot(20*log10(abs(H_true)+1e-10),'LineWidth',1.5,'Color',[1.0 0.6 0.2],'DisplayName','True Channel');
hold on;
plot(20*log10(abs(H_est)+1e-10),'--','LineWidth',1.5,'Color',[0.3 0.8 1.0],'DisplayName','Estimated Channel');
title(sprintf('Channel Check (NMSE = %.2e)', chan_nmse),...
      'Color',[0.9 0.9 1.0],'FontSize',10,'FontWeight','bold');
xlabel('Active Subcarrier Index','Color','white');
ylabel('|H(k)| (dB)','Color','white');
legend('Location','best','FontSize',8);

subplot(1,4,3);
styleAx(gca);
semilogy(abs(errors(:)).^2,'Color',[1.0 0.4 0.2],'LineWidth',1);
title('Squared Error per Symbol','Color',[0.9 0.9 1.0],...
      'FontSize',10,'FontWeight','bold');
xlabel('Symbol Index','Color','white');
ylabel('|e|^2','Color','white');

subplot(1,4,4);
styleAx(gca);
semilogy(SNR_dB_vec, BER_ref,'--','Color',[1.0 0.6 0.2],'LineWidth',1.5,...
         'DisplayName','16-QAM AWGN Ref');
hold on;
semilogy(SNR_dB, BER_sim,'p','Color',[0.2 1.0 0.4],'MarkerSize',15,...
         'MarkerFaceColor',[0.2 1.0 0.4],...
         'DisplayName',sprintf('Operating Point (%d dB)',SNR_dB));
yline(1e-3,'--','Color','white');
ylim([1e-6 1]);
xlim([min(SNR_dB_vec) max(SNR_dB_vec)]);
title('BER Reference + Operating Point','Color',[0.9 0.9 1.0],...
      'FontSize',10,'FontWeight','bold');
xlabel('SNR (dB)','Color','white');
ylabel('Bit Error Rate','Color','white');
legend('Location','best','FontSize',8);

sgtitle(sprintf('ISAC_6 - Comm Receiver | BER=%.2e | EVM=%.1f%% | SE=%.2f bps/Hz | SNR=%d dB',...
        BER_sim, evm_val, SE_op, SNR_dB),...
        'Color','white','FontSize',11,'FontWeight','bold');

H_full_ls = H_est.';
save('workspace_S6.mat');

fprintf('\n  [SAVED] workspace_S6.mat\n');
fprintf('  Next -> run ISAC_7.m\n');
fprintf('============================================================\n');

function styleAx(ax)
    set(ax,'Color',[0.14 0.14 0.18],'XColor',[0.7 0.7 0.8],...
        'YColor',[0.7 0.7 0.8],'GridColor',[0.35 0.35 0.45],'GridAlpha',0.5);
    grid(ax,'on'); box(ax,'on'); hold(ax,'on');
end
