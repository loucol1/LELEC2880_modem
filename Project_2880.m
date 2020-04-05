

function main()
global N;
N = 128;
global Lc;
Lc = 16;
global SNR;
SNR = 50; %[dB]
global nbr_OFDM_symbols;
nbr_OFDM_symbols = 2;
%b= [1+1i, -1-1i];
%a = awgn(b);
s = vecteur_ofdm_symbols();
r = decoder(s);
%test = calcul_error(s, r)

end

%cp est un vecteur de cyclic prefix de taille Lc = 16
function w = OFDM_symbol(cp)
global N;
%cr�ation de la sequence de symboles � envoyer
QAM = [1+1i, 1-1i, -1+1i, -1-1i];
random = randi([1,4], 1, 128);
s = zeros(1,N);
for a = 1:N
    s(a) = QAM(random(a)); 
end

%calcul de la TF inverse
w = ifft(s(1:N));

%ajoute le cyclic prefix
w = [cp, w];
end

function y = vecteur_ofdm_symbols()
global N;
global Lc;
global nbr_OFDM_symbols;
y = zeros(1,nbr_OFDM_symbols*(N+Lc));
y(1:Lc+N) = OFDM_symbol(zeros(1,Lc));
for a = 1:nbr_OFDM_symbols
    y(a*(Lc+N)+1 : (a+1)*(Lc+N)) = OFDM_symbol(y(a*(Lc+N)-Lc+1 : a*(Lc+N)));
end
end

function r = channel(vecteur_ofdm_symbols)
global SNR;
r = awgn(vecteur_ofdm_symbols, SNR);
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

function symbol_error_rate = calcul_error(send, received)
nbr_error = sum(send ~= received);
symbol_error_rate = nbr_error/length(send);
end





