% Set the number of symbols and receivers
num_symbols = 1000;
num_receivers = 2;

% Generate a random sequence of QPSK symbols
qpsk_symbols = 1/sqrt(2) * (2*randi([0 1], 1, num_symbols)-1 + 1j*(2*randi([0 1], 1, num_symbols)-1));

% Add AWGN to the symbols
snr_db = 10; % Signal-to-noise ratio in dB
qpsk_symbols_noisy = awgn(qpsk_symbols, snr_db, 'measured');

% Split the symbols into two branches
qpsk_symbols_branch1 = qpsk_symbols_noisy(1:2:end);
qpsk_symbols_branch2 = qpsk_symbols_noisy(2:2:end);

% Perform selection combining
qpsk_symbols_combined = zeros(1, num_symbols/2);
for i = 1:num_symbols/2
    r1 = abs(qpsk_symbols_branch1(i));
    r2 = abs(qpsk_symbols_branch2(i));
    if r1 > r2
        qpsk_symbols_combined(i) = qpsk_symbols_branch1(i);
    else
        qpsk_symbols_combined(i) = qpsk_symbols_branch2(i);
    end
end

% Calculate the bit error rate (BER)
ber = sum(qpsk_symbols ~= qpsk_symbols_combined);
fprintf('Bit error rate: %g\n', ber);