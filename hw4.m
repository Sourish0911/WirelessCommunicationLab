clc
clear all;
N = 100000; % Number of bits
SNRdB = 0:1:10; % Range of SNR values in dB
SNR = 10.^(SNRdB/10); % SNR values in linear scale
data = randi([0 1], N, 1);
loop = 10;
BER_awgn = zeros(size(SNR));
BER_rayleigh = zeros(size(SNR));
BER_rician = zeros(size(SNR));
BER_nakagami = zeros(size(SNR));
BER_awgnth = zeros(size(SNR));
BER_rayleighth = zeros(size(SNR));
BER_ricianth = zeros(size(SNR));
qam16 = qammod(data, 16); % 16-QAM modulation
for i = 1:length(SNR)
    for lp = 1:loop
        
    N0=1./SNR(i);
    sigma(i)=sqrt(N0/2);
    noise=sigma(i)*(randn(N,1)+1i*randn(N,1));% AWGN
    h = sqrt(1/2)*(randn(N, 1) + 1j*randn(N, 1));
    %h = 1/sqrt(2)*(randn(N,1) + 1i*randn(N,1));
    Yk_gaussian = qam16 + noise;
    Yk_Rayleigh = h.*qam16 + noise;
    %% rician
    k1=10; %Rician factor
    mean=sqrt(k1/(k1+1));% mean
    sigma=sqrt(1/(2*(k1+1)));%variance
    Nr2=randn(N,1)*sigma+mean;
    Ni2=randn(N,1)*sigma;
    %To generate the Rician Random Variable
    h2=sqrt(Nr2.^2+Ni2.^2); %Rician fading coefficient
    Yk_Rician=qam16.*h2+noise; % s means transmitted signal...please take the value as u have taken
    
    %% nakagami -m
    m =2;
    h_naka = sqrt(m)*sqrt(1/2)*(randn(N, 1) + 1i*randn(N, 1));
    %h_naka = transpose(sqrt(m)*sqrt(1/gamrnd(m,1/m,[N,1])));
    YK_nakagami = qam16.*h_naka + noise;
    %% detections
    detected_awgn_qam16 = qamdemod(Yk_gaussian,16);
    detected_rayleigh_qam16 = qamdemod(Yk_Rayleigh./h,16);
    detected_rician_qam16 = qamdemod(Yk_Rician./h2,16);
    detected_nakagami_qam16 = qamdemod(YK_nakagami./h_naka,16);
    %% error
    errors_awgn = sum(detected_awgn_qam16 ~= data);
    errors_rayleigh = sum(detected_rayleigh_qam16 ~= data);
    errors_rician = sum(detected_rician_qam16 ~= data);
    errors_nakagami = sum(detected_nakagami_qam16 ~= data);
    BER_awgn(i) = BER_awgn(i) + errors_awgn/N;
    
    BER_rayleigh(i) = BER_rayleigh(i) + errors_rayleigh/N;
   
    %BER_rayleighth = 0.5*(1 - sqrt(SNR./(1 + SNR)));
    BER_rician(i) = BER_rician(i) + errors_rician/N;
    
    BER_nakagami(i) = BER_nakagami(i) + errors_nakagami/N;
    end
    BER_awgn(i) = BER_awgn(i)/loop;
    BER_rayleigh(i) = BER_rayleigh(i)/loop;
    BER_rician(i) = BER_rician(i)/loop;
    BER_nakagami(i) = BER_nakagami(i)/loop;
    BER_awgnth(i) = BER_awgnth(i) + qfunc(sqrt(2*SNR(i))); % AWGN theoritical
    BER_rayleighth(i) = BER_rayleighth(i) + 1./(4*SNR(i));% Rayleigh Theoretical
    BER_ricianth(i) = BER_ricianth(i) + 0.5*erfc(sqrt((k1*SNR(i))/(k1+SNR(i)))); % rician theoetical
end
figure;
semilogy(SNRdB, BER_awgn, 'r-x');
hold on
semilogy(SNRdB,BER_rayleigh,'b-o');
semilogy(SNRdB,BER_rician,'g-s');
semilogy(SNRdB,BER_nakagami,'m-d');
% semilogy(SNRdB,BER_ricianth,'or');
% semilogy(SNRdB,BER_rayleighth,'or');
% semilogy(SNRdB,BER_awgnth,'or');
legend('gaussian', 'Rayleigh Fading', 'rician','nakagami');