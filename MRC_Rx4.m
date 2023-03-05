N = 2*100000;
bits = randi([0,1],N,1);
snr_db = 0:1:30;


% Map the bits to QPSK symbols
qpsk_symbols = 1/sqrt(2) * (2*bits(1:2:end)-1 + 1j*(2*bits(2:2:end)-1));

% Fading Coeff for 3 Rx
gauss1 = randn(N, 1);
gauss2 = randn(N, 1);
h_releigh_rx1 = 1 * sqrt(gauss1.^2 + gauss2.^2);

gauss_rx2_1 = randn(N, 1);
gauss_rx2_2 = randn(N, 1);
h_releigh_rx2 = 1 * sqrt(gauss_rx2_1.^2 + gauss_rx2_2.^2);

gauss_rx3_1 = randn(N, 1);
gauss_rx3_2 = randn(N, 1);
h_releigh_rx3 = 1 * sqrt(gauss_rx3_1.^2 + gauss_rx3_2.^2);

gauss_rx4_1 = randn(N, 1);
gauss_rx4_2 = randn(N, 1);
h_releigh_rx4 = 1 * sqrt(gauss_rx3_1.^2 + gauss_rx3_2.^2);

qpsk_symbols_rx2 = zeros(N/2,1);
qpsk_symbols_rx1 = zeros(N/2,1);
qpsk_symbols_rx3 = zeros(N/2,1);
qpsk_symbols_rx4 = zeros(N/2,1);
ber_ray = zeros(length(snr_db),1);
%% Modulation Rayleigh
for i = 1:length(snr_db)
    snr_lin = 10^(snr_db(i)/10);
    sd = 1/snr_lin;
    n_rx1 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx2 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx3 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx4 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    qpsk_symbols_combined = zeros(N/2,1);
    for j = 1:N/2
        qpsk_symbols_rx1(j) = h_releigh_rx1(j)*qpsk_symbols(j)+n_rx1(j);
        qpsk_symbols_rx2(j) = h_releigh_rx2(j)*qpsk_symbols(j)+n_rx2(j);
        qpsk_symbols_rx3(j) = h_releigh_rx3(j)*qpsk_symbols(j)+n_rx3(j);
        qpsk_symbols_rx4(j) = h_releigh_rx4(j)*qpsk_symbols(j)+n_rx4(j);
        % Perform maximum ratio combining  
        w1 = h_releigh_rx1(j)*conj(h_releigh_rx1(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w2 = h_releigh_rx2(j)*conj(h_releigh_rx2(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w3 = h_releigh_rx3(j)*conj(h_releigh_rx3(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w4 = h_releigh_rx4(j)*conj(h_releigh_rx4(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        qpsk_symbols_combined(j) = w1*qpsk_symbols_rx1(j)+w2*qpsk_symbols_rx2(j)+w3*qpsk_symbols_rx3(j)+w4*qpsk_symbols_rx4(j);
        
    end
    qpsk_demod_r = qpsk_demod(qpsk_symbols_combined,N/2);
    ber_ray(i) = sum(xor(qpsk_demod_r,bits));
end

% semilogy(snr_db,ber_ray)
%% Nakagami fading
h_naka_rx1 = sqrt(2).*h_releigh_rx1;
h_naka_rx2 = sqrt(2).*h_releigh_rx2;
h_naka_rx3 = sqrt(2).*h_releigh_rx3;
h_naka_rx4 = sqrt(2).*h_releigh_rx4;
qpsk_symbols_rx2 = zeros(N/2,1);
qpsk_symbols_rx3 = zeros(N/2,1);
qpsk_symbols_rx1 = zeros(N/2,1);
qpsk_symbols_rx4 = zeros(N/2,1);
ber_naka = zeros(length(snr_db),1);
qpsk_demod_r = zeros(N,1);
for i = 1:length(snr_db)
    snr_lin = 10^(snr_db(i)/10);
    sd = 1/snr_lin;
    n_rx1 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx2 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx3 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx4 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    qpsk_symbols_combined = zeros(N/2,1);
    for j = 1:N/2
        qpsk_symbols_rx1(j) = h_naka_rx1(j)*qpsk_symbols(j)+n_rx1(j);
        qpsk_symbols_rx2(j) = h_naka_rx2(j)*qpsk_symbols(j)+n_rx2(j);
        qpsk_symbols_rx3(j) = h_naka_rx3(j)*qpsk_symbols(j)+n_rx3(j);
        qpsk_symbols_rx4(j) = h_naka_rx3(j)*qpsk_symbols(j)+n_rx3(j);
        % Perform maximum ratio combining  
        w1 = h_releigh_rx1(j)*conj(h_releigh_rx1(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w2 = h_releigh_rx2(j)*conj(h_releigh_rx2(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w3 = h_releigh_rx3(j)*conj(h_releigh_rx3(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w4 = h_releigh_rx4(j)*conj(h_releigh_rx4(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        qpsk_symbols_combined(j) = w1*qpsk_symbols_rx1(j)+w2*qpsk_symbols_rx2(j)+w3*qpsk_symbols_rx3(j)+w4*qpsk_symbols_rx4(j);
    end
    qpsk_demod_r = qpsk_demod(qpsk_symbols_combined,N/2);
    ber_naka(i) = sum(xor(qpsk_demod_r,bits));
end

%% AWGN
qpsk_symbols_rx2 = zeros(N/2,1);
qpsk_symbols_rx1 = zeros(N/2,1);
qpsk_symbols_rx3 = zeros(N/2,1);
qpsk_symbols_rx4 = zeros(N/2,1);
ber_awgn = zeros(length(snr_db),1);
qpsk_demod_r = zeros(N,1);
for i = 1:length(snr_db)
    snr_lin = 10^(snr_db(i)/10);
    sd = 1/snr_lin;
    n_rx1 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx2 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx3 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    n_rx4 = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    qpsk_symbols_combined = zeros(N/2,1);
    for j = 1:N/2
        qpsk_symbols_rx1(j) = qpsk_symbols(j)+n_rx1(j);
        qpsk_symbols_rx2(j) = qpsk_symbols(j)+n_rx2(j);
        qpsk_symbols_rx3(j) = qpsk_symbols(j)+n_rx3(j);
        qpsk_symbols_rx4(j) = qpsk_symbols(j)+n_rx4(j);
        % Perform maximum ratio combining  
        w1 = h_releigh_rx1(j)*conj(h_releigh_rx1(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w2 = h_releigh_rx2(j)*conj(h_releigh_rx2(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w3 = h_releigh_rx3(j)*conj(h_releigh_rx3(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        w4 = h_releigh_rx4(j)*conj(h_releigh_rx4(j)) / (h_releigh_rx1(j)*conj(h_releigh_rx1(j)) + h_releigh_rx2(j)*conj(h_releigh_rx2(j))+h_releigh_rx3(j)*conj(h_releigh_rx3(j))+h_releigh_rx4(j)*conj(h_releigh_rx4(j)));
        qpsk_symbols_combined(j) = w1*qpsk_symbols_rx1(j)+w2*qpsk_symbols_rx2(j)+w3*qpsk_symbols_rx3(j)+w4*qpsk_symbols_rx4(j);
    end
    qpsk_demod_r = qpsk_demod(qpsk_symbols_combined,N/2);
    ber_awgn(i) = sum(xor(qpsk_demod_r,bits));
end

figure;
semilogy(snr_db, ber_ray, 'r-x');
xlabel('SNR in dB');
ylabel('BER');
hold on
semilogy(snr_db,ber_naka,'b-o');
semilogy(snr_db,ber_awgn,'g-s');
legend('Rayleigh','Nakagami','AWGN')

% Demodulate the QPSK symbols
function bits_demod = qpsk_demod(qpsk_symbols_noisy,num_symbols)
bits_demod = zeros(2*num_symbols,1);
for i = 1:num_symbols
    if real(qpsk_symbols_noisy(i)) >= 0
        bits_demod(2*i-1) = 1;
    else
        bits_demod(2*i-1) = 0;
    end
    if imag(qpsk_symbols_noisy(i)) >= 0
        bits_demod(2*i) = 1;
    else
        bits_demod(2*i) = 0;
    end
end
end

% function H=nak_m(m,Nr,Nt)
% n=zeros(Nr,Nt);
% for i=1:2*m
% n=n+randn(Nr,Nt).^2;
% end
% n=n/(2*m);
% phi=2*pi*rand(Nr,Nt);
% H=(n.^0.5).*cos(phi)+1i*(n.^0.5).*sin(phi);
% end