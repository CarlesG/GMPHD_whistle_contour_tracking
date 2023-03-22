function [Pbin,f,t]= whistle_extraction_from_spectrogram(x,fs,flow,fhigh)
%
%
% FUNTION    : whistle_extraction_from_spectrogram (whistle similarity)
%              Whistle extraction from the thresholded espectrogram after some processing inspired in Gillespie work
%              "Automatic detection and classification of odontocete
%              whistles, Douglas Gillespie, Marjolaine Caillat, and
%              Jonathan Gordon, JASA, 2013."
%
% SYNTAX     : [Pbin]= whistle_extraction_from_spectrogram(x,fs,flow,fhigh)
%
% INPUT      :
%            x         : 1D signal
%            fs        : Sample frequency
%            flow      : low bound of the frequency to analyze
%            fhigh     : high bound of the frequency to analyze
%   
%
% OUTPUT     :
%            
%
%---------------------------------------------------------------------


%Tdeclick=0.0107; % Window time to declicking
Tdeclick=0.010; % Window time: 10 ms to make it equal to the Pyknogram
window=round(Tdeclick*fs);
shift=round(window/2);
x=x(:);


xi_t=zeros(length(x),1);
power=6;
threshold=5;
indices=1:shift:size(x);

for i=1:numel(indices)-2
    xi=x(indices(i):indices(i)+window-1);
    SD=std(xi);
    m=mean(xi); 
    wi=1./(1+((xi-m)/(threshold*SD)).^power);
    if ~isnan(wi) xi_t(indices(i):indices(i)+window-1)=xi.*wi; else xi_t(indices(i):indices(i)+window-1)=xi;end % If clipping wi=NaN
end
    xi=x(indices(i):end);
    SD=std(xi);
    m=mean(xi); 
    wi=1./(1+((xi-m)/(threshold*SD)).^power);
    if ~isnan(wi) xi_t(indices(i):end)=xi.*wi; else xi_t(indices(i):end)=xi;end


%% 2nd ------ Spectrogram Calculation ------

%noverlap=round(75*window/100);  %  Overlap in samples of (25% recomends Guillespie)
noverlap=round(0*window/100);  %  Overlap 0% to make it comaprable with the Pyknogram
xi_t=xi_t-mean(xi_t);
[P,f,t]=spectrogram(xi_t,window,noverlap,2*window,fs);dt=t(2)-t(1);  % Faster than my function
fidx= find((f>=flow)&(f<=fhigh));
f=f(fidx);
P=P(fidx,:);
n=round(t*fs-(window/2-1));  % For compatibility with my tfd function
Pl=10*log10(abs(P));



%% 3rd ------- Median filtered - Denoisisng  ------

num=30;
Pl_aux1=zeros(size(Pl,1)+num*2,size(Pl,2)*(num*2+1));
Pl_aux=[nan(num,size(Pl,2));Pl;nan(num,size(Pl,2))];

for i=0:2*num
    Pl_aux1(:,1+(size(Pl_aux,2)*i):size(Pl_aux,2)*(i+1)) = circshift(Pl_aux, [i-num, 0]);
end

% mm2=median(Pl_aux1,2,'omitnan');
mm2=nanmedian(Pl_aux1,2);  % The same as median(P1,2,'omitnan') for older versions of matlab
matrix_mm2=repmat(mm2(num+1:end-num),1,size(Pl,2));
Pl_median=Pl-matrix_mm2;

 
%% ------- Constant tone removal (Gillespie) ----- %

alfa=0.02;
bf(:,1)=alfa.*Pl_median(:,1);

for i=2:size(Pl_median,2)
    bf(:,i)=alfa*Pl_median(:,i)+(1-alfa)*bf(:,i-1);
end
% 
bf=circshift(bf,[0,1]);
bf(:,1)=zeros(size(Pl_median,1),1);
Pl_denoise=Pl_median-bf;


%%  ------- Gaussian smoothing Kernel (Gillespie) ------ 

% h=[1 2 1;2 4 2;1 2 1]; % original Gillespie work
h=[1 2 1;2 4 2;1 2 1]/4; % modified
Pl_smoothed=filter2(h,Pl_denoise);
  
%% ------- End computing Spectrogram -------- 

P=Pl_smoothed;

% --------------- drawing spectrograms -------------------------- %

% figure(2);imagesc(Pl);set(gca,'Ydir','normal');
% figure(3);imagesc(Pl_median);set(gca,'Ydir','normal');
% figure(4);imagesc(Pl_denoise);set(gca,'Ydir','normal');
% figure(5);imagesc(Pl_smoothed);set(gca,'Ydir','normal');

% ------------------------------------------------------------------- %

%% Whistle similarity with other method

% ---------- whistle_similarity with method 2 to binarize ----------  %

% ------ method 1 ------- %
% Pbin=(P==ones(vent,1)*max(P)); % Get the binary spectrogram with the maximum in each column
% figure(5),imagesc(Pbin);
% Pbin = bwmorph(Pbin,'bridge');
% figure(9),imagesc(Pbin);
% Pbin = bwmorph(Pbin,'bridge');

% ------ method 2 ------- %
%Pbin=(P>ones(size(P,1),1)*max(P)*9/10); % Get the binary spectrogram with a value extracted for each column (9/10*max)
Pbin=(P>ones(size(P))*max(P(:))*4/10); % A different apporach with a maximum value for the whole matrix

%figure(6),imagesc(Pbin);
Pbin = bwmorph(Pbin,'bridge');