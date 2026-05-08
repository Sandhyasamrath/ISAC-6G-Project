clc; clear; close all;
rng(42);

fprintf('============================================================\n');
fprintf('  S1 - SYSTEM PARAMETERS\n');
fprintf('============================================================\n\n');

% RF / carrier
fc     = 28e9;
c      = 3e8;
lambda = c / fc;
BW     = 100e6;

% OFDM grid
Nsc        = 256;
Nguard     = 56;
Nactive    = Nsc - Nguard;
Nsym       = 128;
Ncp        = 32;
Npilot_sym = 4;

% Derived parameters
delta_f     = BW / Nsc;
Tsym_noCP   = Nsc / BW;
Tsym        = (Nsc + Ncp) / BW;
TCPI        = Nsym * Tsym;
Rmax        = c / (2 * delta_f);
delta_R     = c / (2 * BW);
BW_active   = Nactive * delta_f;
delta_R_eff = c / (2 * BW_active);
vmax        = lambda / (4 * Tsym);
delta_v     = lambda / (2 * TCPI);

% Modulation
M      = 4;   % QPSK
k_bits = log2(M);

% MIMO / beamforming
Nt    = 8;
Nr    = 4;
d_ant = lambda / 2;

% Optimization defaults
alpha    = 0.40;
beta     = 0.75;
Ptotal   = 1.0;
lambda_w = 0.5;

% SNR
SNR_dB_vec = -5:2:25;
SNR_dB     = 20;
SNR_lin    = 10^(SNR_dB/10);
noise_var  = Ptotal / SNR_lin;

% Radar targets
target_range    = [80, 150, 280];
target_velocity = [15, -20, 30];
target_rcs      = [1.0, 0.6, 0.8];
target_angle    = [10, -20, 35];
Ntargets        = numel(target_range);

% CA-CFAR (tighter settings)
Ncfar_guard  = 2;
Ncfar_train  = 8;
Pfa          = 1e-5;           % tightened from 1e-4
Ntrain_total = 2 * Ncfar_train;
alpha_cfar   = Ntrain_total * (Pfa^(-1/Ntrain_total) - 1);

% CP interpretation
Tcp     = Ncp / BW;
max_tau = 2 * max(target_range) / c;
Rcp     = c * Tcp / 2;

% Active / pilot allocation
guard_left = floor(Nguard / 2);
Npilot_sc  = round(alpha * Nactive);
Npilot_sc  = max(8, min(Npilot_sc, Nactive - 8));
Ndata_sc   = Nactive - Npilot_sc;

% Use contiguous pilot band so alpha maps cleanly to sensing bandwidth
pilot_start = floor((Nactive - Npilot_sc) / 2) + 1;
pilot_idx   = pilot_start : (pilot_start + Npilot_sc - 1);
data_idx    = setdiff(1:Nactive, pilot_idx);
Ndata_sym   = Nsym - Npilot_sym;

% Power allocation
Pcomm       = beta * Ptotal;
Psense      = (1 - beta) * Ptotal;
Pdata       = Pcomm / max(Ndata_sc, 1);
Ppilot      = Psense / max(Npilot_sc, 1);
Ptrain_full = Ptotal / Nactive;

fprintf('  Carrier frequency   : %.0f GHz\n', fc/1e9);
fprintf('  Bandwidth           : %.0f MHz\n', BW/1e6);
fprintf('  Wavelength          : %.4f mm\n', lambda*1e3);
fprintf('  Subcarriers total   : %d (Active: %d, Guard: %d)\n', Nsc, Nactive, Nguard);
fprintf('  OFDM symbols        : %d (Pilot sym: %d)\n', Nsym, Npilot_sym);
fprintf('  Subcarrier spacing  : %.3f kHz\n', delta_f/1e3);
fprintf('  Symbol dur (no CP)  : %.3f us\n', Tsym_noCP*1e6);
fprintf('  Symbol dur (w/ CP)  : %.3f us\n', Tsym*1e6);
fprintf('  CPI duration        : %.4f ms\n', TCPI*1e3);

fprintf('\n  --- Radar Parameters ---\n');
fprintf('  Max unambig. range  : %.1f m\n', Rmax);
fprintf('  Range res. (theory) : %.3f m\n', delta_R);
fprintf('  Range res. (active) : %.3f m\n', delta_R_eff);
fprintf('  Max velocity        : %.2f m/s\n', vmax);
fprintf('  Velocity resolution : %.4f m/s\n', delta_v);

fprintf('\n  --- Targets ---\n');
for t = 1:Ntargets
    fprintf('  Target %d            : R=%.0f m  v=%+.0f m/s  RCS=%.1f  angle=%+.0f deg\n', ...
        t, target_range(t), target_velocity(t), target_rcs(t), target_angle(t));
end

fprintf('\n  --- CA-CFAR ---\n');
fprintf('  Pfa                 : %.0e\n', Pfa);
fprintf('  Guard cells         : %d\n', Ncfar_guard);
fprintf('  Training cells      : %d\n', Ncfar_train);
fprintf('  alpha_cfar          : %.4f\n', alpha_cfar);

fprintf('\n  --- ISAC Waveform ---\n');
fprintf('  Modulation          : %d-QAM/QPSK equivalent\n', M);
fprintf('  Pilot ratio alpha   : %.2f (%d pilot SC, %d data SC)\n', alpha, Npilot_sc, Ndata_sc);
fprintf('  Power split beta    : %.2f (Pcomm=%.2f, Psense=%.2f)\n', beta, Pcomm, Psense);
fprintf('  Pilot band          : SC %d to %d within active band\n', pilot_idx(1), pilot_idx(end));
fprintf('  Ptrain_full / SC    : %.5f\n', Ptrain_full);
fprintf('  Ppilot / SC         : %.5f\n', Ppilot);
fprintf('  Pdata / SC          : %.5f\n', Pdata);

fprintf('\n  --- CP Condition Check ---\n');
fprintf('  CP length           : %d samples\n', Ncp);
fprintf('  CP duration         : %.3f us\n', Tcp*1e6);
fprintf('  CP covers range     : %.1f m\n', Rcp);
fprintf('  Max target delay    : %.3f us\n', max_tau*1e6);
if max_tau > Tcp
    fprintf('  WARNING             : Radar echoes extend beyond CP.\n');
    fprintf('                        This is acceptable for radar processing.\n');
else
    fprintf('  OK                  : All echoes fit within CP.\n');
end

save('workspace_S1.mat');

fprintf('\n  [SAVED] workspace_S1.mat\n');
fprintf('  Next -> run ISAC_2.m\n');
fprintf('============================================================\n');