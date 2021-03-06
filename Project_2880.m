close all;
clear all;
load('CIR.mat', 'h');
IRchannel = h;
test_H = fft(h, 128);

global N;
N = 128; %number of subcarriers
global Lc;
Lc = 16; %length of cyclic prefix
global nbr_OFDM_symbols;
nbr_OFDM_symbols = 100;
global Pmax;
SNR = 10; %SNR in dB
H = fft(h, N);
[vecteur_ofdm_symbol, vector_data_brut]  = vecteur_ofdm_symbols();
L = length(vecteur_ofdm_symbol);
SNR = 10^(SNR/10);

Esym = sum(abs(vecteur_ofdm_symbol).^2)/(L); %
Pmax = (sum(abs(vecteur_ofdm_symbol).^2))/nbr_OFDM_symbols;
N0 = (Esym)/(SNR*2); %variance du bruit, /2 partie imaginaire et reelle

%%%%%%%%%%%%%%%%%%%%%%%%%%%
SER_SNR()



%give the symbol error rate function of the snr
function SER_SNR()
SNR_vect = 5:0.5:15;
BER =  zeros(1, length(SNR_vect));
for a=1:length(BER)
    BER(a) = OFDM_simulation(SNR_vect(a));
end
figure();
semilogy(SNR_vect, BER);

Pe = erfc(sqrt(3*10.^(SNR_vect/10)/(2*(4-1))));
hold on
semilogy(SNR_vect, Pe);
xlabel("SNR [dB]");
ylabel("Symbol error rate");
title("Symbol error rate in function of the SNR");
legend("experimental curve", "theoretical curve");
end


% simulation of the sending and the receiving of ofdm vector for one
% specific snr
%return the symbol error rate
function symbol_error_rate = OFDM_simulation(SNR)
[vecteur_ofdm_symbol, vector_data_brut] = vecteur_ofdm_symbols();
vecteur_after_chanel = add_awgn_noise(vecteur_ofdm_symbol, SNR);
receive_decode = decoder(vecteur_after_chanel);
symbol_error_rate = calcul_error(vector_data_brut, detection(receive_decode));
end


%produce an ofdm symbol and the ofdm symbol without cp corresponding
%cp is a vector of cyclix prefix of length Lc = 16
function [symbole_ifft, symbole_without_cp] = OFDM_symbol(cp)
global N;
%creation of the symbols to send
QAM = [1+1i, 1-1i, -1+1i, -1-1i];
random = randi([1,4], 1, N);
s = zeros(1,N);
for a = 1:N
    s(a) = QAM(random(a));
end
symbole_without_cp = s;
%computation of the IFT
inter = ifft(s(1:N));

%add of the cyclic prefix
symbole_ifft = [cp, inter];
end

%produce a vector of length nbr_OFDM_symbols
function [vecteur_ofdm_symbol, vector_data_brut]  = vecteur_ofdm_symbols()
global N;
global Lc;
global nbr_OFDM_symbols;

vecteur_ofdm_symbol = zeros(1,nbr_OFDM_symbols*(N+Lc));
%ofdm_without_cp = zeros(1, nbr_OFDM_symbols*N);
[vecteur_ofdm_symbol(1:Lc+N),data_brut] = OFDM_symbol(zeros(1,Lc));
vector_data_brut = data_brut;
%ofdm_without_cp(1 : N) = ofdm_symbol(Lc+1 : (N+Lc));
for a = 1:nbr_OFDM_symbols-1
    [vecteur_ofdm_symbol(a*(Lc+N)+1 : (a+1)*(Lc+N)),data_brut]  = OFDM_symbol(vecteur_ofdm_symbol(a*(Lc+N)-Lc+1 : a*(Lc+N)));
    vector_data_brut = [vector_data_brut, data_brut];
    %ofdm_without_cp(a*N+1 : (a+1)*N) = ofdm_symbol(a*(Lc+N)-N+1 : a*(N+Lc));
end
end


function y = add_awgn_noise(x, SNR)
L = length(x);
SNR = 10^(SNR/10);
Esym = sum(abs(x).^2)/(L);
N0 = Esym/SNR;
noiseSigma = sqrt(N0/2);
n = noiseSigma*(randn(1,L) + 1i*randn(1,L));
y = x + n;
end

%remove the cyclic prefix and compute the fft
function y = decoder(r)
global N;
global Lc;
global nbr_OFDM_symbols;
y = zeros(1, nbr_OFDM_symbols*N);
without_cp = zeros(1, N);
calcul_fft = zeros(1, N);
for a = 1:nbr_OFDM_symbols
    without_cp = r(a*(Lc+N)-N+1 : a*(N+Lc));
    calcul_fft = fft(without_cp);
    y((a-1)*N+1 : a*N) = calcul_fft;
end
end

%dectect in which quadrant is each symbol of vector_complexe
%The value of the return vector are : 1+1i, 1-1i, -1+1i, -1-1i
function adjust_vector = detection(vector_complexe)
adjust_vector = zeros(1, length(vector_complexe));
for a= 1:length(vector_complexe)
    value = vector_complexe(a);
    if(real(value)>0)
        if(imag(value)>0)
            adjust_vector(a) = 1+1i;
        else
            adjust_vector(a) = 1-1i;
        end
    else
        if(imag(value)>0)
            adjust_vector(a) = -1+1i;
        else
            adjust_vector(a) = -1-1i;
        end
    end
end
end

%compute the symbol error rate of the sent sequence
function symbol_error_rate = calcul_error(send, received)
nbr_error = sum(send ~= received);
symbol_error_rate = nbr_error/length(send);
end
