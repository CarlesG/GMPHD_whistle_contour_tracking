function [ Pyk,BW, n ] = pyknogram_freqdomain( x,fs,flow,fhigh, BW, BWoverlap, T )
%
%
% FUNTION    : Pyknogram estiamtion (my own implementation)
%              Inspired in the manuscript :
%              "Speech formant frequency and bandwidth tracking using multiband energy demodulation
%               Alexandros Potamianos and Petros Maragos,  JASA, 1996."
%
%            (c) Ramon Miralles - February 2020.
%
% SYNTAX     : [ Pyk,BW, n ] = pyknogram( x,fs,flow,fhigh, BW, BWoverlap )
%
% INPUT      :
%            x         : 1D signal
%            fs        : Sample frequency
%            flow      : Start frequency of the Gabor filterbank
%            fhigh     : End frequency of the Gabor filterbank
%            BW        : Bandwidth of the Gabor filters
%            BWoverlap : Overlaping in % of the Gabor filters when creating
%                        the filterbank (typically 0%)
%            T         : Duration of the analysis frame (typically 10ms)
%
% OUTPUT     :
%            Pyk       : Matrix containing the Pyknogram
%            BW        : Matrix containing the formants Bandwidth                      
%            n         : Time vector of both of the output matrixes (same
%                        as for vector x)
%                        
%            
%
%---------------------------------------------------------------------


% Creating the Gabor filters of the filter bank
f0=flow:round(BW*(100-BWoverlap)/100):fhigh+BW/2;
n=(-.01:1/fs:.01)';
alfa=BW*sqrt(2*pi);

Gab_bf1=repmat(exp(-alfa^2*n.^2),1,length(f0));
Gab_bf2=cos(2*pi*n*f0);
h=Gab_bf1.*Gab_bf2;

N=length(x)+length(n)-1;
X=fft(x,N);
H=fft(h,N);

% Compute the pyknogram...
%T=10e-3;
TN=round(T*fs);
% FW=zeros(floor((length(x)-1)/TN),length(f0));
% BW=zeros(floor((length(x)-1)/TN),length(f0));
% FW=zeros(floor((N-1)/TN),length(f0));
% BW=zeros(floor((N-1)/TN),length(f0));
FW=zeros(floor(N/TN),length(f0));
BW=zeros(floor(N/TN),length(f0));
for l=1:length(f0);
    y=real(ifft(X.*H(:,l),N));
    [FW(:,l),BW(:,l)]=Fw_estimate( y,fs,T );
end


n=linspace(0,(length(x)-1)/fs,size(FW,1));
Pyk=FW;

end

