%
% The idea of this script is to test how it works the pyknogram for
% extracting cetacean whistles when compared to the thresholded spectrogram
% (D. Guillespie and other similar techniques).
%

clear all;
addpath('../generic_functions'); % pink2

font_size=15;

% TEST SIGNALS...

%[x,fs]=audioread('../example_sounds/delf_mular.wav');tn=(0:length(x)-1)'/fs;
[x,fs]=audioread('../example_sounds/gvi_sample.wav');tn=(0:length(x)-1)'/fs;
%[x,fs]=audioread('../example_sounds/GOZ_S_1_20171225_032016_frag.wav');tn=(0:length(x)-1)'/fs;
%[x,fs]=audioread('../example_sounds/party1.wav');tn=(0:length(x)-1)'/fs;
%[x,fs]=audioread('../example_sounds/GVI_N_1_20190701_034724_fragment.wav');tn=(0:length(x)-1)'/fs;

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
subplot(3,2,1:4)
htfr=imagesc(n,f(fidx),PdB_fzoom);set(gca,'YDir','normal');ylabel('Frequency [kHz]');xlabel('Time [sec.]');  
title(['Spectrogram [',num2str(flow),' Hz - ',num2str(fhigh),' Hz ]']);
yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);set(gca,'FontSize',font_size);
subplot(3,2,5:6);plot(tn,x);axis tight;grid;xlabel('Time [sec.]');set(gca,'FontSize',font_size);

% And now compute the pyknogram and formants bandwidth...
% tic
% [ FW,BW_est, ndraw ] = pyknogram( x,fs,flow,fhigh, BW, BWoverlap,10e-3 );
% toc

tic
[ FW,BW_est, ndraw ] = pyknogram_freqdomain( x,fs,flow,fhigh, BW, BWoverlap,10e-3 );
toc
% Obtain the colors for each pyknogram dot according to the spectrogram...

for k=1:size(FW,1)
    for l=1:size(FW,2)     
            [~,idx_n]=min(abs(n-ndraw(k)));
            [~,idx_k]=min(abs(f-FW(k,l)));
            pyk_color(k,l)=PdB(idx_k,idx_n);       
    end
end


X=repmat(ndraw,1,size(FW,2));
figure(2);
sh=scatter(X(:),FW(:),[],pyk_color(:),'filled');
sh.SizeData=15;
title(['Pyknogram [',num2str(flow),' Hz - ',num2str(fhigh),' Hz ]']);
ylabel('Frequency [kHz]');xlabel('Time [sec.]');
yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
set(gca,'FontSize',font_size);axis([0 max(tn) flow fhigh]);

%
%  TWO DIFFERENT APPROACHES:
%
%       (1) Locate whistles from BW by locating formants with BW lower than a
%       given threshold. -> Higer resolution but slower and more affected
%       by impulsive noise
%
%       (2) Locate whistles from FW filtering the pyknogram when a density
%       is higher than a given threshold -> More robust to impulsive noise
%

% ----------

% OPTION (1) Locate the elements with a Bandwidth lower than a given value (BW_thre)

%BW_thre=200;
BW_thre=BW/4;
idx_BW=find(BW_est<BW_thre);

% Create a struct and store possible whistles
Pyk1(1:length(idx_BW))=struct('time',[],'freq',[],'ampl',[],'label',[],'done',[]); 
kk=num2cell(X(idx_BW));[Pyk1.time]=kk{:};
kk=num2cell(FW(idx_BW));[Pyk1.freq]=kk{:};
kk=num2cell(pyk_color(idx_BW));[Pyk1.ampl]=kk{:};
[Pyk1.label]=deal(0);
[Pyk1.done]=deal(0);

figure(3);
sh=scatter([Pyk1.time],[Pyk1.freq],[],[Pyk1.ampl],'filled');
sh.SizeData=15;
title(['Formants bandwidth <',num2str(BW_thre),' Hz']);
ylabel('Frequency [kHz]');xlabel('Time [sec.]');
yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
set(gca,'FontSize',font_size);axis([0 max(tn) flow fhigh]);

% ----------
% OPTION (2) Remove those points from the pyknogram with a densisty lower than a given
% ammount (frequency separation greater than X Hz)




% X=repmat(ndraw,1,size(FW,2));
% figure(4);
% sh=scatter(X(idx_PYKf),FW(idx_PYKf),[],pyk_color(idx_PYKf),'filled');
% sh.SizeData=15;
% title('Pyknogram density < BW/4 Hz');
% ylabel('Frequency [kHz]');xlabel('Time [sec.]');  
% yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
% set(gca,'FontSize',font_size);axis([0 max(tn) flow fhigh]);
% 
% idx_combined=intersect(idx_BW,idx_PYKf');
% figure(5);
% sh=scatter(X(idx_combined),FW(idx_combined),[],pyk_color(idx_combined),'filled');
% sh.SizeData=15;
% title(['Combined']);
% ylabel('Frequency [kHz]');xlabel('Time [sec.]');  
% yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
% set(gca,'FontSize',font_size);axis([0 max(tn) flow fhigh]);


% Pyk2=struct('time',[],'freq',[],'ampl',[],'label',[],'done',[]); 
% count=1;
% for k=1:size(FW,1)
%     for l=2:size(FW,2)-1
%         if min([abs(FW(k,l)-FW(k,l-1)) abs(FW(k,l)-FW(k,l+1))]) <(BW*0.5*(100-BWoverlap)/100) 
%             Pyk2(count).time=ndraw(k);
%             Pyk2(count).freq=FW(k,l);
%             [~,idx_n]=min(abs(n-ndraw(k)));
%             [~,idx_k]=min(abs(f-FW(k,l)));
%             Pyk2(count).ampl=PdB(idx_k,idx_n);
%             Pyk2(count).label=0;
%             Pyk2(count).done=0;
%             count=count+1;
%         end
%     end
% end

idx_PYKf=[];
for k=1:size(FW,1)
    for l=2:size(FW,2)-1
        if min([abs(FW(k,l)-FW(k,l-1)) abs(FW(k,l)-FW(k,l+1))]) <(BW*0.5*(100-BWoverlap)/100)
            idx_PYKf=[idx_PYKf sub2ind(size(FW),k,l)];
        end
    end
end

% Create a struct and store possible whistles
Pyk2(1:length(idx_PYKf))=struct('time',[],'freq',[],'ampl',[],'label',[],'done',[]); 
kk=num2cell(X(idx_PYKf));[Pyk2.time]=kk{:};
kk=num2cell(FW(idx_PYKf));[Pyk2.freq]=kk{:};
kk=num2cell(pyk_color(idx_PYKf));[Pyk2.ampl]=kk{:};
[Pyk2.label]=deal(0);
[Pyk2.done]=deal(0);

figure(4); sh=scatter([Pyk2.time],[Pyk2.freq],[],[Pyk2.ampl],'filled');
sh.SizeData=15;
title('Pyknogram density < (BW*0.5*(100-BWoverlap)/100 Hz');
ylabel('Frequency [kHz]');xlabel('Time [sec.]');
yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
set(gca,'FontSize',font_size);axis([0 max(tn) flow fhigh]);


idx_combined=intersect(idx_BW,idx_PYKf');

% Create a struct and store possible whistles
Pyk(1:length(idx_combined))=struct('time',[],'freq',[],'ampl',[],'label',[],'done',[]); 
kk=num2cell(X(idx_combined));[Pyk.time]=kk{:};
kk=num2cell(FW(idx_combined));[Pyk.freq]=kk{:};
kk=num2cell(pyk_color(idx_combined));[Pyk.ampl]=kk{:};
[Pyk.label]=deal(0);
[Pyk.done]=deal(0);

figure(5);clf;
%sh=scatter([Pyk.time],[Pyk.freq],[],[Pyk.ampl],'filled');
sh=scatter([Pyk.time],[Pyk.freq],'filled');
sh.SizeData=15;
title(['Scatter plot of the Pyknogram Combination']);
ylabel('Frequency [kHz]');xlabel('Time [sec.]');
yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
set(gca,'FontSize',font_size);axis([0 max(tn) flow fhigh]);
set(gcf,'color','w');box on;grid;


%
%  NOW LETS COMPARE IT WITH THE SPECTROGRAM APPROACH....
%  (SAMARUC inspired functions)
%

[Pbin,f,t]= whistle_extraction_from_spectrogram(x,fs,flow,fhigh);

% imagesc(t,f,Pbin);set(gca,'Ydir','normal');
% title('Spectrogram binarized');
% ylabel('Frequency [kHz]');xlabel('Time [sec.]');
% yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
% set(gca,'FontSize',font_size);axis([0 max(t) flow fhigh]);

[i,j]=find(Pbin==1);

% Create a struct and store possible whistles
Spe(1:length(i))=struct('time',[],'freq',[],'ampl',[],'label',[],'done',[]); 
kk=num2cell(t(j));[Spe.time]=kk{:};
kk=num2cell(f(i)');[Spe.freq]=kk{:};
[Spe.ampl]=deal(1);
[Spe.label]=deal(0);
[Spe.done]=deal(0);

figure(6);clf;
sh=scatter([Spe.time],[Spe.freq],'filled');
sh.SizeData=15;
title('Scatter plot of the binarized spectrogram');
ylabel('Frequency [kHz]');xlabel('Time [sec.]');
yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
set(gca,'FontSize',font_size);axis([0 max(tn) flow fhigh]);
set(gcf,'color','w');box on;grid;


[Pyk_label, Pyk_label2]=whistle_segmentation(Pyk);




