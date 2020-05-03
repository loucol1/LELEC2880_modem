
close all;
clear all;
load('CIR.mat', 'h');
IRchannel = h;


global N;
N = 128;
global Lc;
Lc = 16;
global nbr_OFDM_symbols;
nbr_OFDM_symbols = 100;
global Pmax;
SNR = 10; %dB
H = fft(h, 128);
[vecteur_ofdm_symbol, vector_data_brut]  = vecteur_ofdm_symbols();
L = length(vecteur_ofdm_symbol);
SNR = 10^(SNR/10);

Esym = sum(abs(vecteur_ofdm_symbol).^2)/(L); %
Pmax = (sum(abs(vecteur_ofdm_symbol).^2))/nbr_OFDM_symbols;
N0 = (Esym)/(SNR*2); %variance du bruit, /2 partie imaginaire et rï¿½elle


H_carre = abs(H).^2;
Perror_target = 10^-5;
gamma = (2/3)*(erfcinv(Perror_target/2))^2; %SNR gap
N0 = N0*ones(length(H_carre), 1);
bruit_sur_canal = N0.*gamma./H_carre;
mu = water_level(bruit_sur_canal); %bruit_sur_canal = sigma_n_carre/|H_n|^2
sigma_x_carre = mu*ones(length(bruit_sur_canal),1) - bruit_sur_canal; %Puissance de signal par channel
signe_sigma = sigma_x_carre > 0;
sigma_x_carre = sigma_x_carre .* signe_sigma; %met les valeurs negatives a zero
SNR_n = sigma_x_carre ./ bruit_sur_canal; %SNR par channel
nbr_bits = (1/2)*log2(1+SNR_n/gamma);
figure
bar(mu*ones(1, length(sigma_x_carre)));
hold on
bar(mu*ones(1, length(sigma_x_carre))-sigma_x_carre); %bruit
bit_rate = sum(nbr_bits); %nombre total de bit sur toutes les porteuses

%distribution uniforme
P_uniform = ones(1,128)*Pmax/128;
SNR_uniform = P_uniform./bruit_sur_canal.';
nbr_bit_uniform = 0.5*log2(1+SNR_uniform/gamma);
bit_rate_uniform = sum(nbr_bit_uniform);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Part3
SNR = 20;

%creation training sequence
training_seq = zeros(1,2*N+2*Lc);
first = -1*ones(1,N);
exposant = (0:N-1);
preambule = first.^exposant;
Preambule = ifft(preambule, 128);
training_seq(1:Lc) = Preambule(N-Lc+1:N);
training_seq(Lc+1:N+Lc) = Preambule;
training_seq(N+Lc+1:N+2*Lc) = Preambule(N-Lc+1:N);
training_seq(N+2*Lc+1:2*N+2*Lc) = Preambule;
%%%

%%%
y = conv(training_seq,h);
y = add_awgn_noise(y, SNR);
Y = fft(y);
% estimate_H = Y(Lc+1:N+Lc)./preambule;
% estimate_h = ifft(estimate_H,length(estimate_H));
% figure()
% plot((1:length(estimate_h)), abs(estimate_h))


%%%
h2=zeros(1,128); % canal avec 8 taps et des 0 derriere
h2(1:length(h))=h;
test_send = training_seq(Lc+1:N+Lc); %1 symbole de la training sequence sans cyclic prefix
H = fft(h2,128);
Test_send = fft(test_send, 128);
Y=H.*Test_send;
estimate_H = Y./preambule;
estimate_h = ifft(estimate_H,128);
estimate_h_prime = estimate_h.';
figure()
plot((1:length(estimate_h)), abs(estimate_h))
title('Channel estimation')
forme = abs(h.'-estimate_h(1:length(h))).^2;
MSE_h = sum(abs(h.'-estimate_h(1:length(h))).^2)/length(h);

Nbr_trial = 21;
MSE_vec = zeros(1,Nbr_trial);


for a=(0:Nbr_trial-1)
    MSE_moyenne = 0;
    nbr_rep=250;
    for b=(0:nbr_rep)
        SNR = a;
        h2=zeros(1,128);
        h2(1:length(h))=h;
        test_send = add_awgn_noise(training_seq(Lc+1:N+Lc),SNR);
        H = fft(h2,128);
        Test_send = fft(test_send, 128);
        Y=H.*Test_send;
        estimate_H = Y./preambule;
        estimate_h = ifft(estimate_H,128);
        estimate_H_8=fft(estimate_h(1:8),128);%On sait que h est de longeur 8 donc on ne prend
        %que les 8 premiers taps
        MSE_moyenne = MSE_moyenne + sum(abs(H-estimate_H_8).^2)/length(H);
    end
    
    MSE_vec(a+1) = MSE_moyenne/(nbr_rep+1);
end
figure()
plot((0:Nbr_trial-1), MSE_vec)
title('MSE function of SNR')

%%%% Part 4 - Viterbi decoding
%%


SNR = 13;
u = [0,1,1,0,0,1,1,0,1,0,1];
[x, symbols] = viterbi_encode(u);
u_receive = viterbi_decode(symbols);

% %calcul de la TF inverse
% inter = ifft(symbol);
% 
% %ajoute le cyclic prefix
% cp = zeros(1, Lc);
% symbole_ofdm = [cp, inter];
% 
% %modulation avec le channel + bruit
% y = conv(symbol_ofdm,h);
% y = add_awgn_noise(y, SNR);
% 
% %retire le cyclic prefix
% y = y(Lc+1:end);
% 
% %receiver
% u_receive = viterbi_decode(y)


function [x, symbols] = viterbi_encode(u)
x = zeros(1, 2*length(u));
symbols = zeros(1,length(u));
u_k2 = 0;
u_k1 = 0;
counter = 1;
while(counter<=length(u))
    u_k = u(counter);
    x(1, [2*counter-1, 2*counter]) = [bitxor(u_k, u_k1), bitxor(bitxor(u_k, u_k1), u_k2)];
    x1_x2 = [bitxor(u_k, u_k1), bitxor(bitxor(u_k, u_k1), u_k2)];
    if(x1_x2(1) == 0)
        x1_x2(1) = -1;
    end
    if(x1_x2(2) == 0)
        x1_x2(2) = -1;
    end
    symbols(counter) = x1_x2(1) + 1j*x1_x2(2);
    u_k2 = u_k1;
    u_k1 = u(counter);
    counter = counter + 1;
end
end

function u = viterbi_decode(y)
matrix = 9999*ones(4, length(y)+1, 2);
matrix(1,1,:) = [0, -1]; %[error, chemin(en decimal)]
for column = 1:length(y)
    for ligne = 1:4
        if(matrix(ligne, column, 1)~= -1)
            if ligne == 1
                matrix(1,column+1,:) = [matrix(1,column,1)+abs(y(column)-(-1-1j)), 0];
                matrix(2,column+1,:) = [matrix(1,column,1)+abs(y(column)-(1+1j)), 3];
            end
            if ligne == 2
                error = matrix(2,column,1)+abs(y(column)-(1+1j));
                if(error < matrix(3,column+1,1))
                    matrix(3,column+1,:) = [error, 3];
                end
                error = matrix(2,column,1)+abs(y(column)-(-1-1j));
                if(error < matrix(4,column+1,1))
                    matrix(4,column+1,:) = [error, 0];
                end
            end
            if ligne == 3
                error = matrix(3,column,1)+abs(y(column)-(-1+1j));
                if(error < matrix(1,column+1,1))
                    matrix(1,column+1,:) = [error,1];
                end
                error = matrix(3,column,1)+abs(y(column)-(1-1j));
                if(error < matrix(2,column+1,1))
                    matrix(2,column+1,:) = [error, 2];
                end
            end
            
            if ligne == 4
                error = matrix(4,column,1)+abs(y(column)-(1-1j));
                if(error < matrix(3,column+1,1))
                    matrix(3,column+1,:) = [error,2];
                end
                error = matrix(4,column,1)+abs(y(column)-(-1+1j));
                if(error < matrix(4,column+1,1))
                    matrix(4,column+1,:) = [error, 1];
                end
            end
        end
    end
end
x = zeros(1,length(y)); %chemin avec somme des erreur minimal
[Min, Indice] = min(matrix(:,end,1));
for column = length(y):1
    if Indice == 1
        if matrix(1,column,2)==0
            Indice = 1;
            x(1,column) = 0;
        else
            Indice = 3;
            x(1,column) = 1;
        end
    else if Indice == 2
            if matrix(2,column,2)==3
                Indice = 1;
                x(1,column) = 3;
            else
                Indice = 3;
                x(1,column) = 2;
            end
        else if Indice == 3
                if matrix(3,column,2) == 3
                    Indice = 2;
                    x(1,column) = 3;
                else
                    Indice = 4;
                    x(1,column) = 2;
                end
            else
                if matrix(4,column,2) == 0
                    Indice = 2;
                    x(1,column) = 0;
                else
                    Indice = 4;
                    x(1,column) = 1;
                end
            end
        end
    end
end

%retrouver u à partir de x
u = zeros(1,length(x));
u_k1 = 0;
for k = 1:length(x)
    u(k) = bitxor(u_k1, floor(x(k)/2));
    u_k1 = u(k);
end
end




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
%crï¿½ation de la sequence de symboles ï¿½ envoyer
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
