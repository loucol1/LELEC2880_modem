%% test
clear all;
a=[1+2i,5-2i,6+1i];
b = [1+2i,5-2i,6+1i,0,0,0];
Test = fft(a);
Test2 = fft(b,6);