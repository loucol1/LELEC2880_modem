close all;
clear all;
load('CIR.mat', 'h');
IRchannel = h;
test_H = fft(h, 128);

global N;
N = 128;%number of subcarriers
global Lc;
Lc = 16; % length of the cyclic prefix
global nbr_OFDM_symbols;
nbr_OFDM_symbols = 100;
global Pmax;
SNR = 10; % SNR in dB
H = fft(h, 128);
[vecteur_ofdm_symbol, vector_data_brut]  = vecteur_ofdm_symbols();
L = length(vecteur_ofdm_symbol);
SNR = 10^(SNR/10);

Esym = sum(abs(vecteur_ofdm_symbol).^2)/(L); %
Pmax = (sum(abs(vecteur_ofdm_symbol).^2))/nbr_OFDM_symbols;
N0 = (Esym)/(SNR*2); %variance du bruit, /2 partie imaginaire et reelle


H_carre = abs(H).^2;
Perror_target = 10^-5;
gamma = (2/3)*(erfcinv(Perror_target/2))^2; %SNR gap
N0 = N0*ones(length(H_carre), 1);
bruit_sur_canal = N0.*gamma./H_carre; %influence of the channel, the noise and the gamma needed for the computation of mu
mu = water_level(bruit_sur_canal); %bruit_sur_canal = sigma_n_carre/|H_n|^2
sigma_x_carre = mu*ones(length(bruit_sur_canal),1) - bruit_sur_canal; %Power on each subcarrier
signe_sigma = sigma_x_carre > 0;
sigma_x_carre = sigma_x_carre .* signe_sigma; %put the negative value to zero
SNR_n = sigma_x_carre ./ bruit_sur_canal; %SNR per subcarrier
nbr_bits = (1/2)*log2(1+SNR_n/gamma);
figure
bar(mu*ones(1, length(sigma_x_carre)), 'r');
hold on
bar(mu*ones(1, length(sigma_x_carre))-sigma_x_carre, 'b'); %effect of the channel, the noise and the gamma on a subcarrier
bit_rate = sum(nbr_bits); %total number of bit on all the subcarriers 

%uniform distribution 
P_uniform = ones(1,128)*Pmax/128;
SNR_uniform = P_uniform./bruit_sur_canal.';
nbr_bit_uniform = 0.5*log2(1+SNR_uniform/gamma);
bit_rate_uniform = sum(nbr_bit_uniform);

%recursive function that return the water level
function mu = water_level(bruit_sur_canal)
global Pmax;
mu = (Pmax + sum(bruit_sur_canal))/length(bruit_sur_canal);
[Max, Indice] = max(bruit_sur_canal);
while(mu<Max)
    bruit_sur_canal(Indice) = [];
    mu = (Pmax + sum(bruit_sur_canal))/length(bruit_sur_canal);
    [Max, Indice] = max(bruit_sur_canal);
end
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


