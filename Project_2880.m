

clear all;
global N;
N = 128;
global Lc;
Lc = 16;
global nbr_OFDM_symbols;
nbr_OFDM_symbols = 5;
%b= [1+1i, -1-1i];
%a = awgn(b);

SNR_vect = 0:0.005:40;
BER =  zeros(1, length(SNR_vect));
for a=1:length(BER)
    BER(a) = OFDM_simulation(SNR_vect(a));
end
semilogy(SNR_vect, BER);





function bit_error_rate = OFDM_simulation(SNR)
[vecteur_ofdm_symbol, vector_data_brut] = vecteur_ofdm_symbols();
vecteur_after_chanel = channel(vecteur_ofdm_symbol, SNR);
receive_decode = decoder(vecteur_after_chanel);
symbol_error_rate = calcul_error(vector_data_brut, detection(receive_decode));
bit_error_rate = symbol_error_rate*2;
end



%cp est un vecteur de cyclic prefix de taille Lc = 16
function [symbole_ifft, symbole_without_cp] = OFDM_symbol(cp)
global N;
%création de la sequence de symboles à envoyer
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
SNR_naturel = 10^(SNR/20);
amplitude = sqrt(2);%l'amplitude d un symbole
noise_variance = amplitude^2/SNR_naturel;
noise_variance = noise_variance/2;
vector_noise_real = normrnd(0, noise_variance, [1, length(vector_clean)]);
vector_noise_imaginary = 1i*normrnd(0, noise_variance, [1, length(vector_clean)]);
vector_noise = vector_noise_real + vector_noise_imaginary;
vector_corrupted = vector_clean + vector_noise;
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






