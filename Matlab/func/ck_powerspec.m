function [P,freqs] = ck_powerspec(X,samp)

% Quick Power spectrum based on FFT
% 
% function [Power,faxis] = powerspec(data,samp)
% computes the power spectrum of the signal X
% sampled at rate 'samp' [Hz]

L = length(X);
if rem(L,2)==1
  L = L-1;
end

freqs = samp*([0:L/2-1])/L;
F = fft(X([1:L]));
P = F.*conj(F)/L;
P = P([1:L/2]);