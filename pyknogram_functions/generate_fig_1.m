%
% The idea of this script is to test how it works the pyknogram for
% extracting cetacean whistles when compared to the thresholded spectrogram
% (D. Guillespie and other similar techniques).
%

clear all;
addpath('../generic_functions'); % pink2
addpath('/Users/rmiralle/rmr/math/export_fig/');

font_size=15;

% TEST SIGNALS...

[x,fs]=audioread('../example_sounds/GVI_N_1_20190701_034724_fragment.wav');tn=(0:length(x)-1)'/fs;

% Add Pink Noise according to SNR
% SNR=12;  % Signal to noise ratio in dB
% sig_pow=mean(x.^2);
% pnoise=pink2(length(x))';
% noise_pow=mean(pnoise.^2);
% scalef=sqrt( (sig_pow*10^(-SNR/10))/noise_pow);
% pnoise=pnoise*scalef;  % Apply scale factor to get the desired SNR
% x=(x+pnoise);

% Settings for the definition of the filter Bank
BW=1000; % BW in Hz
BWoverlap=50; % BW in %
flow=3000;fhigh=min([22000 fs/2-BW/2]); % MINIMUM flow -> flow=round(BW/2)


% Compute & draw the Spectrogram

figure(1);
vent=512;solape=400;
[P,f,n]=spectrogram(x,vent,solape,[],fs);
PdB=20*log10(abs(P));
fidx= find((f>=flow)&(f<=fhigh));
P_fzoom=P(fidx,:);
PdB_fzoom=20*log10(abs(P_fzoom));

figure(1);clf;set(gcf,'color','w');
htfr=imagesc(n,f(fidx),PdB_fzoom);set(gca,'YDir','normal');ylabel('Frequency [kHz]');xlabel('Time [sec.]');  
title(['Spectrogram [',num2str(flow),' Hz - ',num2str(fhigh),' Hz ]']);
axis tight;xlabel({'Time [sec.]';'(a)'});
yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
set(gca,'FontSize',font_size);

saveas(gcf,'spectrogram_example.png');
% And now compute the pyknogram and formants bandwidth...

tic
[ FW,BW_est, ndraw ] = pyknogram_freqdomain( x,fs,flow,fhigh, BW, BWoverlap,10e-3 );
toc

X=repmat(ndraw,1,size(FW,2));
figure(2);clf;set(gcf,'color','w');
sh=scatter(X(:),FW(:),'filled');box on;
sh.SizeData=15;
title(['Pyknogram [',num2str(flow),' Hz - ',num2str(fhigh),' Hz ]']);
ylabel('Frequency [kHz]');xlabel({'Time [sec.]';'(b)'});
yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
set(gca,'FontSize',font_size);axis([0 max(tn) flow fhigh]);

print2eps('pyknogram_example')