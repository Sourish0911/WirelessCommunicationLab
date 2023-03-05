N = 100000;
bits = randi([0,1],N,1);
symbol = zeros(N,1);
y_releigh = zeros(N,1);
y_releigh_dec = zeros(N,1);
for i=1:length(bits)
    if bits(i) == 0
        symbol(i)=-1;
    else
        symbol(i)=1;
    end
end
snr_db = -5:1:10;
% Generate two independent Gaussian random variables
gauss1 = randn(N, 1);
gauss2 = randn(N, 1);
h_releigh = 1 * sqrt(gauss1.^2 + gauss2.^2);
% BPSK
ber_gaus = zeros(length(snr_db),1);
for j = 1:length(snr_db)
    snr_lin = 10^(snr_db(j)/10);
    sd = 1/snr_lin;
    n = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    for i = 1:length(bits)
        y_releigh(i) = 1*symbol(i);
        y_releigh(i) = y_releigh(i)+n(i);
            if real(y_releigh(i))>=0
                y_releigh_dec(i)=1;
            else
                y_releigh_dec(i) = 0;
            end
            if bits(i)~=y_releigh_dec(i)
                ber_gaus(j)=ber_gaus(j)+1;
            else
                continue
            end
    end
end

ber_ray = zeros(length(snr_db),1);
for j = 1:length(snr_db)
    snr_lin = 10^(snr_db(j)/10);
    sd = 1/snr_lin;
    n = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    for i = 1:length(bits)
        y_releigh(i) = h_releigh(i)*symbol(i);
        y_releigh(i) = y_releigh(i)+n(i);
            if real(y_releigh(i))>=0
                y_releigh_dec(i)=1;
            else
                y_releigh_dec(i) = 0;
            end
            if bits(i)~=y_releigh_dec(i)
                ber_ray(j)=ber_ray(j)+1;
            else
                continue
            end
    end
end

ber_rice = zeros(length(snr_db),1);
h_rice = ricean(N,2);
for j = 1:length(snr_db)
    snr_lin = 10^(snr_db(j)/10);
    sd = 1/snr_lin;
    n = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    for i = 1:length(bits)
        y_releigh(i) = h_rice(i)*symbol(i);
        y_releigh(i) = y_releigh(i)+n(i);
            if real(y_releigh(i))>=0
                y_releigh_dec(i)=1;
            else
                y_releigh_dec(i) = 0;
            end
            if bits(i)~=y_releigh_dec(i)
                ber_rice(j)=ber_rice(j)+1;
            else
                continue
            end
    end
end

ber_nak = zeros(length(snr_db),1);
h_nak = nak_m(2,N,1);
for j = 1:length(snr_db)
    snr_lin = 10^(snr_db(j)/10);
    sd = 1/snr_lin;
    n = sqrt(sd)*(randn(N, 1)+1i*randn(N, 1));
    for i = 1:length(bits)
        y_releigh(i) = h_nak(i)*symbol(i);
        y_releigh(i) = y_releigh(i)+n(i);
            if real(y_releigh(i))>=0
                y_releigh_dec(i)=1;
            else
                y_releigh_dec(i) = 0;
            end
            if bits(i)~=y_releigh_dec(i)
                ber_nak(j)=ber_nak(j)+1;
            else
                continue
            end
    end
end
subplot(4,1,1)
semilogy(snr_db,ber_gaus)
legend('Gaussian');
xlabel('SNR in dB');
ylabel('BER');
hold on
subplot(4,1,2)
semilogy(snr_db,ber_ray)
legend('rayleigh');
xlabel('SNR in dB');
ylabel('BER');
hold on
subplot(4,1,3)
semilogy(snr_db,ber_rice)
legend('Ricean');
xlabel('SNR in dB');
ylabel('BER');
hold on
subplot(4,1,4)
semilogy(snr_db,ber_nak)
legend('Nakagami');
xlabel('SNR in dB');
ylabel('BER');
hold on
