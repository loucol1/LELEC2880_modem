close all;
clear all;
load('CIR.mat', 'h');
IRchannel = h;
H = fft(h, 128);
H = H.';

global N;
N = 128;
global Lc;
Lc = 16;
global nbr_OFDM_symbols;
nbr_OFDM_symbols = 100;
global Pmax;
global L;
L=8;%longeur du canal

%%Part3
SNR = 20;

%creation training sequence
training_seq = zeros(1,2*N);
first = -1*ones(1,2*N);
exposant = (0:2*N-1);
preambule = first.^exposant;
%%%

bruit = add_awgn_noise(preambule, SNR)-preambule;
y_without_noise=conv(preambule, h);
y_without_noise=y_without_noise(1:end-(length(h)-1));
y = y_without_noise + bruit;
T = zeros(N-L+1,L);
count=1;
for a=(L:-1:1)
    T(:,count) = preambule(a:N-count+1);
    count = count+1;
end
y=y(L:N);
y=y.';
%estimate_h = inv(T)*y;
test = T.\T






% Y_1 = fft(y(1:N),128);
% Y_2 = fft(y(N+1:end),128);
% Preambule1 = fft(preambule(1:N),N);
% Preambule2 = fft(preambule(N+1:end),N);
% H_1 = Y_1./Preambule1;
% H_2 = Y_2./Preambule2;
% H_moy = (H_1+H_2)/2;
% 
% 
% 
% 
% 
% 
% 
% 
% Nbr_trial = 21;
% save_estimate_H = zeros(Nbr_trial, 128);
% MSE_vec = zeros(1,Nbr_trial);
% 
% 
% for a=(0:Nbr_trial-1)
%     MSE_moyenne = 0;
%     nbr_rep=250;
%     for b=(0:nbr_rep)
%         SNR = a;
%         bruit = add_awgn_noise(preambule, SNR)-preambule;
%         y_without_noise=conv(preambule, h);
%         y_without_noise=y_without_noise(1:end-(length(h)-1));
%         y = y_without_noise + bruit;
%         Y_1 = fft(y(1:N),128);
%         Y_2 = fft(y(N+1:end),128);
%         H_1 = Y_1./fft(preambule(1:N),N);
%         H_2 = Y_2./fft(preambule(N+1:end),N);
%         H_moy = (H_1+H_2)/2;
%         MSE_moyenne = MSE_moyenne + sum(abs(H-H_moy).^2)/length(H);
%     end
%     MSE_vec(a+1) = MSE_moyenne/(nbr_rep+1);
%     
%     %rempli la matrice save
%     save_estimate_H(a+1, :) = H_moy;
% end
% save('save_estimate_H');
% 
% figure()
% plot((0:Nbr_trial-1), MSE_vec)
% title('MSE function of SNR')

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

