clear; close all; clc;

%% Simulation Parameters
B = 160e6;              % Bandwidth [Hz]
N_sub = 2048;           % Number of sub-carriers [samples]
L_cp = 256;             % Cycle Prefix length [samples]
L_preamble = 2;         % Number of repetition of preamble [OFDM symbols]
L_data = 30;            % Data length, Comm. Mode [OFDM symbols]
fc = 5e9;               % Carrier frequency [Hz]
SNR_dB = -5:5:30;       % SNR [dB]

%%  Transmitter
%   Bits Generation (Data)
M = 16;                             % Modulation order
k = log2(M);                        % Bits per symbol
num_bits = N_sub * L_data * k;      % Total number of bits
tx_bits = randi([0 1], num_bits, 1);

%   Map bits to QAM symbols (Modulation)
tx_syms = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
tx_syms_grid = reshape(tx_syms, N_sub, L_data);

%   Preamble Generation
preamble_bits = randi([0 1], N_sub, 1);                 % Random bits
preamble_sym = pskmod(preamble_bits, 2);                % BPSK
preamble_grid = repmat(preamble_sym, 1, L_preamble);    % Used to replicate the preamble L_preamble times

%   Concatenate in Frequency Domain
tx_grid_full = [preamble_grid, tx_syms_grid];           % Used to concatenate the preamble and the data in the frequency domain

%   OFDM Modulation (IFFT) - Full Frame
tx_time_full = ifft(tx_grid_full, N_sub);               % Convert frequency-domain symbols to time-domain waveforms discretized in time with N_sub samples

%   Add Cyclic Prefix
tx_cp = [tx_time_full(end-L_cp+1:end, :); tx_time_full];

%   Serialize
tx_signal = tx_cp(:);

%%  Channel
ber_sim = zeros(length(SNR_dB), 1);


%%  Plots
%   First 50 bits
%   Constellation of Transmitted Symbols
%   Time Domain Signal (First 200 samples)

figure('Name', 'Transmitter Visualization', 'Color', 'w', 'Position', [100, 100, 1200, 400]);

%   First 50 bits
subplot(1,3,1);
stem(tx_bits(1:50), 'filled', 'LineWidth', 1.5, 'Color', '#0072BD');
title('First 50 Transmitted Bits', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Bit Index', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Value', 'Interpreter', 'latex', 'FontSize', 12);
ylim([-0.2 1.2]);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

%   Symbols (Constellation)
subplot(1,3,2);
plot(real(tx_syms), imag(tx_syms), 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'none');
title('Transmitted Constellation (16-QAM)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('In-Phase (I)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Quadrature (Q)', 'Interpreter', 'latex', 'FontSize', 12);
grid on; axis square;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

%   Signal (Time Domain)
subplot(1,3,3);
t_plot = (0:199); % First 200 samples
plot(t_plot, real(tx_signal(1:200)), 'LineWidth', 1.2, 'Color', '#77AC30');
title('Transmitted Signal (Time Domain)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Sample Index', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12);
grid on; xlim([0 200]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

%   Save the plot
saveas(gcf, 'Transmitter_Viz.png');

