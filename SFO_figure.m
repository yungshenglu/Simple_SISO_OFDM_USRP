% Plot two figure with SFO and without SFO correction
SFO = decode(1, 1);
nonSFO = decode(0, 0);

% Plot Tx data
file = fopen('tx_vec.bin', 'rb');
complex_array = fread(file, 'float');

tx_vec = complex_array(1 : 2 : end) + 1i * complex_array(2 : 2 : end);
tx_vec = abs(tx_vec) .^ 2;

figure;
plot(tx_vec);
title('Tx data');

% Plot Rx signals
file = fopen('rx_signals.bin', 'rb');
complex_array = fread(file, 'float');

rx_signals = complex_array(1 : 2 : end) + 1i * complex_array(2 : 2 : end);
rx_signals = abs(rx_signals) .^ 2;

figure;
plot(rx_signals);
title('Raw Rx signal');

% Plot filtered Rx signal
rx_vec_air = read_complex_binary('rx_signals.bin');
tx_offset = 108400;
raw_rx_dec = rx_vec_air(tx_offset : 12000 + tx_offset).';
raw_rx_dec = abs(raw_rx_dec) .^ 2;

figure;
plot(raw_rx_dec);
title('Filtered Rx signal');

% Plot phases of decoded signal
figure;
plot(SFO, 'o');
hold on;
plot(nonSFO, 'x');
ylim([0, 3.5]);
title('Phases of decoded signal');










