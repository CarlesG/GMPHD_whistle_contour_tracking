function [Fw,Bw] = Fw_estimate( x,fs,T )
%
%   Formant Frequency and Bandwidth Short Time Estimates
%
%   [Fw,Bw] = Fw_estimate(  x,fs,T )
%
%   Ramon Miralles (c) 2020
%   Based in "Speech Formant frequency and bandwidth tracking using
%   multiband energy demodulation, JASA, A.Potamianos, 1996.
%

%T=10e-3;
TN=round(T*fs);
round_length=floor(length(x)/TN)*TN;
x=x(1:round_length);  % Adjust the lenght of x so that it is a multiple of the block size TN

% % ESA estimate (Poor performance for the low first formant -> makes more
% % whistle harmonics appear
%  psi  = teager_keiser(x);
%  psidot  = teager_keiser( gradient(x) );
%  fe=real(1/(2*pi)*sqrt(psidot./psi))*fs;  
%  a=real(psi./sqrt(psidot));   

 
% Hilbert Trasnsform Demodulation (HTD) Estimate
fe=gradient(unwrap(angle(hilbert(x))))*fs/(2*pi);  
a=abs(hilbert(x));


% Estimate the derivative of a
adot=gradient(a);


% Estimate of the Formant Frequency (using second moments)
aux1=fe.*(a.^2);
AUX1=sum(reshape(aux1,TN,length(aux1)/TN));
aux2=a.^2;aux2=aux2(1:round_length);
AUX2=sum(reshape(aux2,TN,length(aux1)/TN));
Fw=AUX1./AUX2;

% Estimate of the Bandwidth (using second moments)
Fw_aux=repmat(Fw,TN,1);Fw_rep=Fw_aux(:);
aux1=(adot/(2*pi)).^2+(fe-Fw_rep).^2.*a.^2;
AUX1=sum(reshape(aux1,TN,length(aux1)/TN));
Bw=sqrt(AUX1./AUX2);

end

