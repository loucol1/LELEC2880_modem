%%
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




%%%% Part 4 - Viterbi decoding
nbr_diff_vect = zeros(1,16);
nbr_diff_vect2 = zeros(1,16);
nbr_diff_vect3 = zeros(1,16);
Nbr_trial_max = 1000;
for SNR = 0:15
    nbr_diff = 0;
    nbr_diff2 = 0;
    nbr_diff3 = 0;
    for nbr_trial = 1:Nbr_trial_max
        u = random_vector(128);
        
        [x, symbols] = viterbi_encode(u);
        
        
        %calcul de la TF inverse
        inter = ifft(symbols, 128);
        
        %ajoute le cyclic prefix
        cp = zeros(1, Lc);
        symbole_ofdm = [cp, inter];
        
        %modulation avec le channel + bruit
        y = conv(symbole_ofdm,h);
        
        y_without_viterbi = add_awgn_noise(y, SNR+10*log10(2));
        y = add_awgn_noise(y, SNR);
        
        
        %retire le cyclic prefix
        %Y = fft(y(Lc+1:end),128);
        %y = y(Lc+1:end-7);
        H_1=fft(h,128);
        H_1 = H_1.';
        %receiver
        
        Y_first = fft(y(Lc+1:end-7), 128);
        %Y_first = Y_first.';
        Y_second = H_1.*fft(inter, 128);
        %Y_second = Y_second.';
        Y_without_viterbi = fft(y_without_viterbi(Lc+1:end-7), 128);
        
        %Viterbi
        u_receive = viterbi_decode(Y_first, H_1);
        nbr_diff = nbr_diff+sum(abs(u_receive-u));
        
        %Viterbi sans channel
        u_receive3 = viterbi_decode(Y_first./H_1, ones(1,128));
        nbr_diff3 = nbr_diff3+sum(abs(u_receive3-u));
        
        %sans Viterbi
        Y_equalize = Y_without_viterbi./H_1;
        Y_decode = detection(Y_equalize);
        x_receive = symbol_to_bit(Y_decode);
        nbr_diff2 = nbr_diff2+sum(abs(x_receive-x));
    end
    nbr_diff_vect(SNR+1) = nbr_diff/(128*Nbr_trial_max);
    nbr_diff_vect2(SNR+1) = nbr_diff2/(128*2*Nbr_trial_max);
    nbr_diff_vect3(SNR+1) = nbr_diff3/(128*Nbr_trial_max);
end
figure
semilogy((0:15),nbr_diff_vect);
hold on
semilogy((0:15),nbr_diff_vect2);
hold on
semilogy((0:15),nbr_diff_vect3);
legend('Viterbi', 'QAM 4','Viterbi without channel');
xlabel('E_s/N_0 (dB)');
ylabel('BER');
title('BER function of SNR with and without Viterbi coding');





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

function u = viterbi_decode(y, H)
matrix = 99999*ones(4, length(y)+1, 2);
matrix(1,1,:) = [0, -1]; %[error, chemin(en decimal)]
for column = 1:length(y)
    for ligne = 1:4
        if(matrix(ligne, column, 1)~= 99999)
            if ligne == 1
                err = abs(y(column)-(-1-1j)*H(column))^2;
                %err = err*division(h(column));
                matrix(1,column+1,:) = [matrix(1,column,1)+err, 0];
                err = abs(y(column)-(1+1j)*H(column))^2;
                %err = err*division(h(column));
                matrix(2,column+1,:) = [matrix(1,column,1)+err, 3];
            end
            if ligne == 2
                err = abs(y(column)-(1+1j)*H(column))^2;
                %err = err*division(h(column));
                error = matrix(2,column,1)+err;
                if(error < matrix(3,column+1,1))
                    matrix(3,column+1,:) = [error, 3];
                end
                err = abs(y(column)-(-1-1j)*H(column))^2;
                %err = err*division(h(column));
                error = matrix(2,column,1)+err;
                if(error < matrix(4,column+1,1))
                    matrix(4,column+1,:) = [error, 0];
                end
            end
            if ligne == 3
                err = abs(y(column)-(-1+1j)*H(column))^2;
                %err = err*division(h(column));
                error = matrix(3,column,1)+err;
                if(error < matrix(1,column+1,1))
                    matrix(1,column+1,:) = [error,1];
                end
                err = abs(y(column)-(1-1j)*H(column))^2;
                %err = err*division(h(column));
                error = matrix(3,column,1)+err;
                if(error < matrix(2,column+1,1))
                    matrix(2,column+1,:) = [error, 2];
                end
            end
            
            if ligne == 4
                err = abs(y(column)-(1-1j)*H(column))^2;
                %err = err*division(h(column));
                error = matrix(4,column,1)+err;
                if(error < matrix(3,column+1,1))
                    matrix(3,column+1,:) = [error,2];
                end
                err = abs(y(column)-(-1+1j)*H(column))^2;
                %err = err*division(h(column));
                error = matrix(4,column,1)+err;
                if(error < matrix(4,column+1,1))
                    matrix(4,column+1,:) = [error, 1];
                end
            end
        end
    end
end

x = zeros(1,length(y)); %chemin avec somme des erreur minimal
[Min, Indice] = min(matrix(:,end,1));
for column = linspace(length(y)+1,2,length(y))
    if Indice == 1
        if matrix(1,column,2)==0
            Indice = 1;
            x(1,column-1) = 0;
        else
            Indice = 3;
            x(1,column-1) = 1;
        end
    else if Indice == 2
            if matrix(2,column,2)==3
                Indice = 1;
                x(1,column-1) = 3;
            else
                Indice = 3;
                x(1,column-1) = 2;
            end
        else if Indice == 3
                if matrix(3,column,2) == 3
                    Indice = 2;
                    x(1,column-1) = 3;
                else
                    Indice = 4;
                    x(1,column-1) = 2;
                end
            else
                if matrix(4,column,2) == 0
                    Indice = 2;
                    x(1,column-1) = 0;
                else
                    Indice = 4;
                    x(1,column-1) = 1;
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

function un_sur_h = division(h)
renv = abs(h-1);
if renv<0.01
    renv = 0.01;
end
un_sur_h = 1/renv;
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

function x = random_vector(Dim)
A = [0,1];
random = randi([1,2],1,Dim);
x = zeros(1,Dim);
for b= 1:Dim
    x(b) = A(random(b));
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

function vector_bit = symbol_to_bit(vector_symbol)
    vector_bit = zeros(1,2*length(vector_symbol));
    for a = 1:length(vector_symbol)
        if vector_symbol(a)==1+1j
            vector_bit(2*a-1:2*a) = [1, 1];
        end
        if vector_symbol(a)==1-1j
            vector_bit(2*a-1:2*a) = [1, 0];
        end
        if vector_symbol(a)==-1+1j
            vector_bit(2*a-1:2*a) = [0, 1];
        end
        if vector_symbol(a)==-1-1j
            vector_bit(2*a-1:2*a) = [0, 0];
        end
    end
end