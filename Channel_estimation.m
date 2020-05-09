close all;
clear all;
load('CIR.mat', 'h');
IRchannel = h;
test_H = fft(h, 128);

global N;
N = 128;
global Lc;
Lc = 16;
global nbr_OFDM_symbols;
nbr_OFDM_symbols = 100;
global Pmax;

%%Part3
SNR = 20;

%creation training sequence
training_seq = zeros(1,2*N+2*Lc);
first = -1*ones(1,N);
exposant = (0:N-1);
preambule = first.^exposant
Preambule = ifft(preambule,128)
training_seq(1:Lc) = Preambule(N-Lc+1:N);
training_seq(Lc+1:N+Lc) = Preambule;
training_seq(N+Lc+1:N+2*Lc) = Preambule(N-Lc+1:N);
training_seq(N+2*Lc+1:2*N+2*Lc) = Preambule;
%%%

%%%
test_send = training_seq(1:N+Lc); %1 symbole de la training sequence avec cyclic prefix
y = conv(test_send,h);
%y_second = conv(test_send_second, h);
%y = add_awgn_noise(y, SNR);
Y = fft(y(Lc+1:end-7),128);
%Y_second = fft(y_second(N+2*Lc+1:end-7),128); %retire le cyclic prefix du deuxieme symbole_ofdm
H_e_test = Y./preambule;
estimate_h_test = ifft(H_e_test, 128);
estimate_h_test = estimate_h_test.';
% estimate_H = Y(Lc+1:N+Lc)./preambule;
% estimate_h = ifft(estimate_H,length(estimate_H));
% figure()
% plot((1:length(estimate_h)), abs(estimate_h))


%%%
h2=zeros(1,128); % canal avec 8 taps et des 0 derriere
h2(1:length(h))=h;
test_send = training_seq(Lc+1:N+Lc); %1 symbole de la training sequence sans cyclic prefix
test_send_second = training_seq(N+Lc+1:end);
H = fft(h2,128);
Test_send = fft(test_send, 128);
Test_send_second = fft(test_send_second, 128);
Y=H.*Test_send;
Y_second = H.*Test_send_second;
estimate_H = Y./preambule;
estimate_H_second = Y_second./preambule;
estimate_h = ifft(estimate_H,128);
estimate_h_second = ifft(estimate_H_second,128);
estimate_h_prime = estimate_h.';
% figure()
% plot((1:length(estimate_h)), abs(estimate_h))
% title('Channel estimation')
forme = abs(h.'-estimate_h(1:length(h))).^2;
MSE_h = sum(abs(h.'-estimate_h(1:length(h))).^2)/length(h) + sum(abs(h.'-estimate_h_second(1:length(h))).^2)/length(h);
MSE_h = MSE_h/2;

Nbr_trial = 21;
save_estimate_H = zeros(Nbr_trial, 128);
MSE_vec = zeros(1,Nbr_trial);


for a=(0:Nbr_trial-1)
    MSE_moyenne = 0;
    nbr_rep=250;
    for b=(0:nbr_rep)
        SNR = a;
        h2=zeros(1,128);
        h2(1:length(h))=h;
        test_send = add_awgn_noise(training_seq(Lc+1:N+Lc),SNR);
        test_send_second = add_awgn_noise(training_seq(N+2*Lc+1:end),SNR);
        H = fft(h2,128);
        Test_send = fft(test_send, 128);
        Test_send_second = fft(test_send_second, 128);
        Y=H.*Test_send;
        Y_second = H.*Test_send_second;
        estimate_H = Y./preambule;
        estimate_H_second = Y_second./preambule;
        estimate_h = ifft(estimate_H,128);
        estimate_h_second = ifft(estimate_H_second,128);
        estimate_H_8=fft(estimate_h(1:8),128);%On sait que h est de longeur 8 donc on ne prend
        %que les 8 premiers taps
        estimate_H_8_second=fft(estimate_h_second(1:8),128);
        MSE_moyenne = MSE_moyenne + sum(abs(H-estimate_H_8).^2)/length(H) +sum(abs(H-estimate_H_8_second).^2)/length(H);
    end
    MSE_moyenne = MSE_moyenne/2;
    MSE_vec(a+1) = MSE_moyenne/(nbr_rep+1);
    
    %rempli la matrice save
    save_estimate_H(a+1, :) = (estimate_H_8+estimate_H_8_second)/2;
end
save('save_estimate_H');

figure()
plot((0:Nbr_trial-1), MSE_vec)
title('MSE function of SNR')

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

