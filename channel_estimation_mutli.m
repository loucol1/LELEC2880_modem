%% channel estimation, multiple training sequence

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
snr = 20;
nbr_iter = 100;
Nbr_symbol = 10;
save_snr_H = zeros(snr, N, Nbr_symbol);
MSE = zeros(snr,Nbr_symbol); 




for nbr_symbol = 1:Nbr_symbol
    for SNR = 0:snr
        estimate_H_moyen = zeros(1,N);
        for iter = 1:nbr_iter
            matrix = ones(Nbr_symbol, N);
            
            %creation du training sequence
            exposant = (0:N-1);
            Ik = (-1*ones(1,N)).^exposant;
            matrix = matrix.*Ik;
            
            %convolution avec le channel + bruit
            matrix_conv = matrix.*(repmat(H, Nbr_symbol,1));
            for a = [1:Nbr_symbol]
                matrix_conv(a,:) = add_awgn_noise(matrix_conv(a,:), snr);
            end
            
            estimate_H = matrix_conv./matrix;
            estimate_H_sum = sum(estimate_H,1)/Nbr_symbol;
            
            estimate_H_moyen = estimate_H_moyen + estimate_H_sum;
        end
        %pour chaque snr on a un estimé de channel
        estimate_H_moyen = estimate_H_moyen/nbr_iter;
        MSE(snr+1,nbr_symbol) = sum(abs(estimate_H_moyen-H).^2);
    end
end
plot((0:SNR), MSE(:,1))

    








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