%%coucou sacha

%cp est un vecteur de cyclic prefix de taille Lc = 16
function y = OFDM_symbol(cp)
N = 128;
Lc = 16;
SNR = 10; %[dB]

%création de la sequence de symboles à envoyer
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

%ajoute le AWGN
r = awgn(w, SNR);

%retire le cyclic prefix
y = r(Lc:N+Lc);
end
