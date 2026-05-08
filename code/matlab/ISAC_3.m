%% =========================================================================
%  ISAC_3.m
%  SECTION 3 - TRANSMITTER: OFDM-ISAC GRID GENERATION
%  =========================================================================
%  Requires: workspace_S2.mat
%  Saves to: workspace_S3.mat
%  Next step -> run ISAC_4.m
%  =========================================================================

clc; close all;
load('workspace_S2.mat');

fprintf('============================================================\n');
fprintf('  ISAC_3 - TRANSMITTER: OFDM-ISAC GRID\n');
fprintf('============================================================\n\n');

% -------------------------------------------------------------------------
%  STEP 1: GENERATE RANDOM COMMUNICATION BITS
% -------------------------------------------------------------------------
Nbits   = Ndata_sc * Ndata_sym * k_bits;
tx_bits = randi([0 1], Nbits, 1);

fprintf('  Total bits          : %d\n', Nbits);
fprintf('  Data subcarriers    : %d\n', Ndata_sc);
fprintf('  Data OFDM symbols   : %d\n', Ndata_sym);

% -------------------------------------------------------------------------
%  STEP 2: 16-QAM MODULATION
% -------------------------------------------------------------------------
tx_symbols = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
fprintf('  QAM symbols         : %d  (16-QAM)\n', length(tx_symbols));

% -------------------------------------------------------------------------
%  STEP 3: OFDM GRID DESIGN
%
%  FIX:
%    1. First Npilot_sym symbols are FULL training symbols on ALL active SC.
%       These are used for robust communication channel estimation.
%    2. Remaining symbols carry:
%         - pilot_idx  : known sensing/reference tones
%         - data_idx   : 16-QAM communication data
%
%  This keeps communication estimation stable and preserves radar pilots
%  across the rest of the frame.
% -------------------------------------------------------------------------
Sgrid = zeros(Nactive, Nsym);

% Full training symbol power spread across all active subcarriers
Ptrain_full = Ptotal / Nactive;

% First Npilot_sym OFDM symbols: all-active known training
for s = 1:Npilot_sym
    Sgrid(:, s) = sqrt(Ptrain_full) * ones(Nactive, 1);
end

% Remaining symbols: pilots on pilot_idx, data on data_idx
tx_sym_mat   = reshape(tx_symbols, Ndata_sc, Ndata_sym);
pilot_symbol = ones(Npilot_sc, 1);

for s = Npilot_sym+1:Nsym
    Sgrid(pilot_idx, s) = sqrt(Ppilot) * pilot_symbol;
    Sgrid(data_idx,  s) = sqrt(Pdata)  * tx_sym_mat(:, s - Npilot_sym);
end

fprintf('\n  Pilot / Training structure:\n');
fprintf('    Full training sym : %d  (all %d active SC known)\n', Npilot_sym, Nactive);
fprintf('    Pilot SC in data  : %d  (known sensing tones)\n', Npilot_sc);
fprintf('    Data  SC in data  : %d  (16-QAM)\n', Ndata_sc);
fprintf('    Ptrain_full / SC  : %.5f\n', Ptrain_full);
fprintf('    Ppilot / SC       : %.5f\n', Ppilot);
fprintf('    Pdata  / SC       : %.5f\n', Pdata);

% -------------------------------------------------------------------------
%  STEP 4: MAP TO FULL SUBCARRIER GRID
% -------------------------------------------------------------------------
Sfull = zeros(Nsc, Nsym);
Sfull(guard_left+1 : guard_left+Nactive, :) = Sgrid;

% -------------------------------------------------------------------------
%  STEP 5: IFFT + CP
% -------------------------------------------------------------------------
tx_td = zeros(Nsc + Ncp, Nsym);
for s = 1:Nsym
    ofdm_sym    = ifft(ifftshift(Sfull(:, s)), Nsc);
    tx_td(:, s) = [ofdm_sym(end-Ncp+1:end); ofdm_sym];
end

fprintf('\n  IFFT size           : %d\n', Nsc);
fprintf('  CP length           : %d samples\n', Ncp);
fprintf('  tx_td size          : [%d x %d]\n', size(tx_td,1), size(tx_td,2));

% -------------------------------------------------------------------------
%  PLOTS
% -------------------------------------------------------------------------
figure('Name','ISAC_3 - OFDM-ISAC Grid','Color',[0.12 0.12 0.15],...
       'Position',[100 100 900 420]);

subplot(1,2,1);
set(gca,'Color',[0.14 0.14 0.18],'XColor',[0.7 0.7 0.8],...
    'YColor',[0.7 0.7 0.8],'GridColor',[0.35 0.35 0.45],'GridAlpha',0.5);
hold on; grid on;

grid_power = abs(Sfull(:, Npilot_sym+1)).^2;
bar(1:Nsc, grid_power, 'FaceColor',[0.3 0.8 1.0],...
    'EdgeColor','none','DisplayName','Occupied SC');

pilot_full_idx = guard_left + pilot_idx;
bar(pilot_full_idx, abs(Sfull(pilot_full_idx, Npilot_sym+1)).^2,...
    'FaceColor',[1.0 0.6 0.2],'EdgeColor','none','DisplayName','Pilot SC');

title('Data-Symbol Subcarrier Occupancy','Color',[0.9 0.9 1.0],...
      'FontSize',10,'FontWeight','bold');
xlabel('Subcarrier Index','Color','white');
ylabel('Power','Color','white');
legend('Location','best','FontSize',8);

subplot(1,2,2);
set(gca,'Color',[0.14 0.14 0.18],'XColor',[0.7 0.7 0.8],...
    'YColor',[0.7 0.7 0.8],'GridColor',[0.35 0.35 0.45],'GridAlpha',0.5);
hold on; grid on;

t_axis = (0:Nsc+Ncp-1) / BW * 1e6;
plot(t_axis, real(tx_td(:,1)),'Color',[0.3 0.8 1.0],'LineWidth',1,...
     'DisplayName','Re{s(t)}');
plot(t_axis, imag(tx_td(:,1)),'Color',[1.0 0.5 0.2],'LineWidth',1,...
     'DisplayName','Im{s(t)}');
xline(Ncp/BW*1e6,'--','Color','white');

title('Time-Domain OFDM Symbol (with CP)','Color',[0.9 0.9 1.0],...
      'FontSize',10,'FontWeight','bold');
xlabel('Time (us)','Color','white');
ylabel('Amplitude','Color','white');
legend('Location','best','FontSize',8);

sgtitle('ISAC_3 - OFDM-ISAC Transmitter','Color','white',...
        'FontSize',12,'FontWeight','bold');

% -------------------------------------------------------------------------
%  SAVE
% -------------------------------------------------------------------------
save('workspace_S3.mat');

fprintf('\n  [SAVED] workspace_S3.mat\n');
fprintf('  Next -> run ISAC_4.m\n');
fprintf('============================================================\n');
