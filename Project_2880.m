

clear all;
load('CIR.mat', 'h');
IRchannel = h;


global N;
N = 128;
global Lc;
Lc = 16;
global nbr_OFDM_symbols;
nbr_OFDM_symbols = 10;
global Pmax;
Pmax = 1;
SNR = 12; %dB
H = fft(h, 128);
[vecteur_ofdm_symbol, vector_data_brut]  = vecteur_ofdm_symbols();
L = length(vecteur_ofdm_symbol);
SNR = 10^(SNR/10);
Esym = sum(abs(vecteur_ofdm_symbol).^2)/(L);
N0 = Esym/SNR; %variance du bruit
H_carre = abs(H).^2;
N0 = N0*ones(length(H_carre), 1);
bruit_sur_canal = N0./H_carre;
mu = water_level(bruit_sur_canal); %bruit_sur_canal = sigma_n_carre/|H_n|^2
sigma_x_carre = mu*ones(length(bruit_sur_canal),1) - bruit_sur_canal;
signe_sigma = sigma_x_carre > 0;
sigma_x_carre = sigma_x_carre .* signe_sigma; %met les valeurs negatives a zero
SNR = sigma_x_carre ./ N0;
Perror_target = 10^-5;
gamma = (2/3)*(erfcinv(Perror_target/2))^2;
nbr_bits = (1/2)*log2(1+SNR/gamma);




function mu = water_level(bruit_sur_canal)
global Pmax;
length_debut = length(bruit_sur_canal)
mu = (Pmax + sum(bruit_sur_canal))/length(bruit_sur_canal);
[Max, Indice] = max(bruit_sur_canal);
while(mu<Max)
    bruit_sur_canal(Indice) = [];
    mu = (Pmax + sum(bruit_sur_canal))/length(bruit_sur_canal);
    [Max, Indice] = max(bruit_sur_canal);
end
length_fin = length(bruit_sur_canal)
end





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
vecteur_after_chanel = channel(vecteur_ofdm_symbol, SNR);
receive_decode = decoder(vecteur_after_chanel);
symbol_error_rate = calcul_error(vector_data_brut, detection(receive_decode));
end


%produce an ofdm symbol and the ofdm symbol without cp corresponding
%cp est un vecteur de cyclic prefix de taille Lc = 16
function [symbole_ifft, symbole_without_cp] = OFDM_symbol(cp)
global N;
%cr�ation de la sequence de symboles � envoyer
QAM = [1+1i, 1-1i, -1+1i, -1-1i];
random = randi([1,4], 1, N);
s = zeros(1,N);
for a = 1:N
    s(a) = QAM(random(a)); 
end
symbole_without_cp = s;
%calcul de la TF inverse
inter = ifft(s(1:N));

%ajoute le cyclic prefix
symbole_ifft = [cp, inter];
end

%produce a vector of nbr_OFDM_symbols
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


function vector_corrupted = channel(vector_clean, SNR)
vector_corrupted = add_awgn_noise(vector_clean, SNR);
% SNR_naturel = 10^(SNR/20);
% amplitude = sqrt(2);%l'amplitude d un symbole
% noise_variance = amplitude^2/SNR_naturel;
% noise_variance = noise_variance/2;
% vector_noise_real = normrnd(0, noise_variance, [1, length(vector_clean)]);
% vector_noise_imaginary = 1i*normrnd(0, noise_variance, [1, length(vector_clean)]);
% vector_noise = vector_noise_real + vector_noise_imaginary;
% vector_corrupted = vector_clean + vector_noise;
end

function y = add_awgn_noise(x, SNR)
L = length(x);
SNR = 10^(SNR/10);
Esym = sum(abs(x).^2)/(L);
N0 = Esym/SNR;
% if(isreal(x))
%     noiseSigma = sqrt(N0);
%     n = noiseSigma*randn(1,L);
% else
    noiseSigma = sqrt(N0/2);
    n = noiseSigma*(randn(1,L) + 1i*randn(1,L));
%end
y = x + n;
end


function y = decoder(r)
global N;
global Lc;
global nbr_OFDM_symbols;
y = zeros(1, nbr_OFDM_symbols*N);
without_cp = zeros(1, N);
calcul_fft = zeros(1, N);
%retire le cyclic prefix et calcule la fft
for a = 1:nbr_OFDM_symbols
    without_cp = r(a*(Lc+N)-N+1 : a*(N+Lc));
    calcul_fft = fft(without_cp);
    y((a-1)*N+1 : a*N) = calcul_fft;
end
end

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

function symbol_error_rate = calcul_error(send, received)
nbr_error = sum(send ~= received); 
symbol_error_rate = nbr_error/length(send);
end






