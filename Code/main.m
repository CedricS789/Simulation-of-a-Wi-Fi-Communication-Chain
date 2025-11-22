clear; close all; clc;

%% 1.2.1 - Simulation Parameters
B = 160e6;              % Bandwidth [Hz]
N_sub = 2048;           % Number of sub-carriers [samples]
L_cp = 256;             % Cycle Prefix length [samples]
L_preamble = 2;         % Number of repetition of preamble [OFDM symbols]
L_data = 30;            % Data length, Comm. Mode [OFDM symbols]
fc = 5e9;               % Carrier frequency [Hz]
SNR_dB = -5:5:30;       % SNR [dB]

%% 1.2.2 - Transmitter
% Bits Generation (Data)
M = 16;                             % Modulation order
k = log2(M);                        % Bits per symbol
num_bits = N_sub * L_data * k;      % Total number of bits
tx_bits = randi([0 1], num_bits, 1);

% Map bits to QAM symbols (Modulation)
tx_syms = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
tx_syms_matrix = reshape(tx_syms, N_sub, L_data);

% Preamble Generation
preamble_bits = randi([0 1], N_sub, 1);                 % Random bits
preamble_sym = pskmod(preamble_bits, 2);                % BPSK
preamble_matrix = repmat(preamble_sym, 1, L_preamble);  % Used to replicate the preamble L_preamble times (S/P Conversation)

% Concatenate in Frequency Domain
tx_matrix_full = [preamble_matrix, tx_syms_matrix];     % Used to concatenate the preamble and the data in the frequency domain

% OFDM Modulation (IFFT) - Full Frame
tx_time_full = ifft(tx_matrix_full, N_sub);             % Convert frequency-domain symbols to time-domain waveforms discretized in time with N_sub samples

% Add Cyclic Prefix
tx_cp = [tx_time_full(end-L_cp+1:end, :); tx_time_full];

% Serialize
tx_signal = tx_cp(:);                                   % (P/S Conversion)

%% 1.2.3 - Receiver Chain (ideal)
rx_signal = tx_signal;                                  % Assume ideal channel (rx_signal = tx_signal)  

% S/P Conversion (Reshape)
rx_matrix_cp = reshape(rx_signal, N_sub + L_cp, L_preamble + L_data);

% CP Removal
rx_matrix = rx_matrix_cp(L_cp+1:end, :);

% FFT
rx_matrix_full = fft(rx_matrix, N_sub);

% Extract Data (Discard Preamble)
rx_syms_matrix = rx_matrix_full(:, L_preamble+1:end);

% P/S Conversion (Serialize)
rx_syms = rx_syms_matrix(:);

% QAM Demapping
rx_bits = qamdemod(rx_syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);

% BER Calculation
[num_errors, ber] = biterr(tx_bits, rx_bits);
fprintf('Bit Error Rate (BER): %.4e\n', ber);

% Validation: Verify BER is 0
if ber == 0
    disp('Validation Successful: BER is exactly 0.');
else
    warning('Validation Failed: BER is %.4e (should be 0).', ber);
end

%% Plots
% First 50 bits
% Constellation of Transmitted Symbols
% Time Domain Signal (First 200 samples)

figure('Name', 'Transmitter Visualization', 'Color', 'w', 'Position', [100, 100, 1200, 400]);

% First 50 bits
subplot(1,3,1);
stem(tx_bits(1:50), 'filled', 'LineWidth', 1.5, 'Color', '#0072BD');
title('First 50 Transmitted Bits', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Bit Index', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Value', 'Interpreter', 'latex', 'FontSize', 12);
ylim([-0.2 1.2]);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

% Symbols (Constellation)
subplot(1,3,2);
plot(real(tx_syms), imag(tx_syms), 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'none');
title('Transmitted Constellation (16-QAM)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('In-Phase (I)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Quadrature (Q)', 'Interpreter', 'latex', 'FontSize', 12);
grid on; axis square;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

% Signal (Time Domain)
subplot(1,3,3);
t_plot = (0:199); % First 200 samples
plot(t_plot, real(tx_signal(1:200)), 'LineWidth', 1.2, 'Color', '#77AC30');
title('Transmitted Signal (Time Domain)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Sample Index', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12);
grid on; xlim([0 200]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

% Save the plot
saveas(gcf, 'Transmitter_Viz.png');

figure('Name', 'Receiver Visualization', 'Color', 'w', 'Position', [100, 550, 800, 400]);

% Received Constellation
subplot(1,2,1);
plot(real(rx_syms), imag(rx_syms), 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'none');
title('Received Constellation', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('In-Phase (I)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Quadrature (Q)', 'Interpreter', 'latex', 'FontSize', 12);
grid on; axis square;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

% First 50 Received Bits
subplot(1,2,2);
stem(rx_bits(1:50), 'filled', 'LineWidth', 1.5, 'Color', '#0072BD');
title('First 50 Received Bits', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Bit Index', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Value', 'Interpreter', 'latex', 'FontSize', 12);
ylim([-0.2 1.2]);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

% Save the plot
saveas(gcf, 'Receiver_Viz.png');

