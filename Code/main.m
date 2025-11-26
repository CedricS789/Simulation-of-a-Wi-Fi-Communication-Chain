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

% Check orthogonality between 2 sub carriers
sub_1_time = tx_time_full(1, :);   % Sub-carrier 1 (time domain)
sub_700_time = tx_time_full(700, :); % Sub-carrier 700 (time domain)
sub_1200_time = tx_time_full(1200, :); % Sub-carrier 1200 (time domain)
sub_1400_time = tx_time_full(1400, :); % Sub-carrier 1400 (time domain)
sub_1600_time = tx_time_full(1600, :); % Sub-carrier 1600 (time domain)

sub_1_freq = tx_matrix_full(1, :);   % Sub-carrier 1 (frequency domain)
sub_700_freq = tx_matrix_full(700, :); % Sub-carrier 700 (frequency domain)
sub_1200_freq = tx_matrix_full(1200, :); % Sub-carrier 1200 (frequency domain)
sub_1400_freq = tx_matrix_full(1400, :); % Sub-carrier 1400 (frequency domain)
sub_1600_freq = tx_matrix_full(1600, :); % Sub-carrier 1600 (frequency domain)

% Integral between sub_1 and sub_700
% Correct syntax: trapz(x) or trapz(x, y)
integral = trapz(sub_1_time .* conj(sub_700_time));
integral_1200 = trapz(sub_1_time .* conj(sub_1200_time));
integral_1400 = trapz(sub_1_time .* conj(sub_1400_time));
integral_1600 = trapz(sub_1_time .* conj(sub_1600_time));

fprintf('Integral between sub_1 and sub_700: %.2e\n', abs(integral));
fprintf('Integral between sub_1 and sub_1200: %.2e\n', abs(integral_1200));
fprintf('Integral between sub_1 and sub_1400: %.2e\n', abs(integral_1400));
fprintf('Integral between sub_1 and sub_1600: %.2e\n', abs(integral_1600));

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
fprintf('BER = %.4e\n', ber);


%%  Receiver Chain (AWGN Channel)
ber_sim = zeros(length(SNR_dB), 1);

for i = 1:length(SNR_dB)
    EbNo = SNR_dB(i);

    %   Convert Eb/N0 to SNR
    %   SNR = (Eb/N0) * (Bits per Symbol) * (Code Rate)
    %   Code Rate here accounts for CP overhead: N_sub / (N_sub + L_cp)
    snr = EbNo + 10*log10(k) + 10*log10(N_sub / (N_sub + L_cp));

    %   AWGN Channel
    rx_signal = awgn(tx_signal, snr, 'measured');

    %   S/P Conversion (Reshape)
    rx_matrix_cp = reshape(rx_signal, N_sub + L_cp, L_preamble + L_data);

    %   CP Removal
    rx_matrix = rx_matrix_cp(L_cp+1:end, :);

    %   FFT
    rx_grid_full = fft(rx_matrix, N_sub);

    %   Extract Data (Discard Preamble)
    rx_syms_grid = rx_grid_full(:, L_preamble+1:end);

    %   P/S Conversion (Serialize)
    rx_syms = rx_syms_grid(:);

    %   QAM Demapping
    rx_bits = qamdemod(rx_syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);

    %   BER Calculation
    [~, ber_sim(i)] = biterr(tx_bits, rx_bits);
end

%% Receiver Chain (Multipath Channel)
h = rand(1, 3); % Random multipath channel

rx_signal_multipath = filter(h, 1, tx_signal); % Apply channel

% S/P Conversion (Reshape)
rx_matrix_cp_multipath = reshape(rx_signal_multipath, N_sub + L_cp, L_preamble + L_data);

% CP Removal
rx_matrix_multipath = rx_matrix_cp_multipath(L_cp+1:end, :);

% FFT
rx_matrix_full_multipath = fft(rx_matrix_multipath, N_sub);

% Extract Data (Discard Preamble)
rx_syms_matrix_multipath = rx_matrix_full_multipath(:, L_preamble+1:end);

% P/S Conversion (Serialize)
rx_syms_multipath = rx_syms_matrix_multipath(:);

%%  Theoretical BER
ber_theo = berawgn(SNR_dB, 'qam', M);


%% Frequency Domain Channel Equalization with Zero Forcing
% Transform the time domain h into frequency domain
H_fft = fft(h, N_sub);

% Expand the channel to match data dimension
H_matrix = repmat(H_fft.', 1, L_data);

% Apply Zero Forcing Equalization
rx_syms_matrix_equalized = rx_syms_matrix_multipath ./ H_matrix;

% P/S Conversion (Serialize)
rx_syms_equalized = rx_syms_matrix_equalized(:);


%%  Plots
%   First 50 bits
%   Constellation of Transmitted Symbols
%   Time Domain Signal (First 200 samples)
%   BER vs Eb/N0

%   Transmitter Visualization
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

%   Receiver Visualization
figure('Name', 'Receiver Visualization', 'Color', 'w', 'Position', [100, 550, 800, 400]);

%   Received Constellation (Last SNR)
subplot(1,2,1);
plot(real(rx_syms), imag(rx_syms), 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'none');
title(['Received Constellation ($E_b/N_0$ = ' num2str(SNR_dB(end)) ' dB)'], 'Interpreter', 'latex', 'FontSize', 14);
xlabel('In-Phase (I)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Quadrature (Q)', 'Interpreter', 'latex', 'FontSize', 12);
grid on; axis square;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

%   BER vs Eb/N0
subplot(1,2,2);
semilogy(SNR_dB, ber_theo, 'k-', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB, ber_sim, 'bo--', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
title('BER vs $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$E_b/N_0$ [dB]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('BER', 'Interpreter', 'latex', 'FontSize', 12);
legend('Theoretical', 'Simulated', 'Interpreter', 'latex', 'FontSize', 11);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

%   Save the plot
saveas(gcf, 'Receiver_Viz.png');

%   Multipath Visualization
figure('Name', 'Multipath Visualization', 'Color', 'w', 'Position', [100, 100, 1000, 500]);

%   Before Equalization
subplot(1,2,1);
plot(real(rx_syms_multipath), imag(rx_syms_multipath), 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', '#7E2F8E', 'MarkerEdgeColor', 'none');
title('Before Equalization', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('In-Phase (I)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Quadrature (Q)', 'Interpreter', 'latex', 'FontSize', 12);
grid on; axis square;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

%   After Equalization
subplot(1,2,2);
plot(real(rx_syms_equalized), imag(rx_syms_equalized), 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', '#77AC30', 'MarkerEdgeColor', 'none');
title('After Zero-Forcing Equalization', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('In-Phase (I)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Quadrature (Q)', 'Interpreter', 'latex', 'FontSize', 12);
grid on; axis square;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

%   Save the plot
saveas(gcf, 'Multipath_Viz.png');

