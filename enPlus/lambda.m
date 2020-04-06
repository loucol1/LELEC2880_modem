clear all
load('CIR.mat', 'h');
IRchannel = h;
first = TF_h(0, 128, h);
second = TF_h(1, 128, h);
trois = TF_h(2, 128, h);
quatre = TF_h(4, 128, h);
l = lam(10^-5);
ldb = 10*log10(l);
function lambda = lam(Proba_erreur)
first = erfcinv(Proba_erreur/2)^2;
lambda = 2*first/3;
end

function H = TF_h(k,N, h)
H=0;
for a=1:length(h)
    H = H+h(a)*exp(-1i*2*pi*k*(a-1)/N);
end
end