%% GMPHD DETECTOR
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This script runs a GMPHD detector to track frequency contours of whistles.
% The output of the detector is stored in a struct called DT which has 
% 3 fields for each detected whistle contour: freq (frequency), time, label
% (a unique ID of a given contour).

% Pina Gruden, Institute of Sound and Vibration Research (ISVR), 2014-2018

%-------------------------------------------------------------------------
% For detailed explanation see: Gruden, P. and White, P. (2016). Automated
% tracking of dolphin whistles using Gaussian mixture probability
% hypothesis density filters. Journal of the Acoustical Society of America,
% 140(3),1981-1991.
%-------------------------------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear, close all
addpath('pyknogram_functions'); % All pyknogram functions

%% ////////////// Specify PARAMETERS //////////////
%NOTE: The filter was developed based on using window length (win_width_s) 
% of 10.7 ms (which resulted in 93.75 Hz frequency bin resolution). 
% The time increment (dt) between the windows was 5.35 ms (50 % overlap).
% This is set as default.
% If you decide to change this, then you will have to change variances of
% the GM-PHD filter, specifically the model.Q (system noise covariance will 
% have to be retrained.

% ------------ Parameters for obtaining the measurements ------------------
win_width_s = 0.0053; %window length in seconds (default)
freqrange = [5000,50000]; % lower and higher frequency range limits in
% Hz in which your signals occur 
peak_thr =8; %threshold in dB for finding spectral peaks 
%(these become your measurements that will go to the detector)

% ---------  PARAMETERS for the GMPHD detector ----------------------
%--------------Parameters for System and observation Models --------------
% System model: xk =F*xk_1 + v; 
% Measurement model: zk = H*xk + w;
dt=0.00535; %time increment in s (default)
% dt=0.0100; %time increment in s (default)
models.F = [1, dt; 0, 1]; %state transition matrix
models.Q =[5013,166770;166770, 53973000]; %system noise covariance matrix
models.H = [1, 0]; % measurement matrix
models.R= round((1/win_width_s)^2/12) ;%measurement noise covariance matrix

%-------------Parameters for birth,clutter and other--------------------
nClutter=10; %number of clutter points per time step
parameters.wth=0.009; %threshold used in state estimation
parameters.pdet = 0.85; %probability of detection (used in Update calculation)
parameters.psurv = 0.994;%probability of survival
parameters.clutter = nClutter/(freqrange(2)-freqrange(1)); % Clutter intensity
parameters.Jmax = 100; %max number of Gaussian componenets used in pruning
parameters.U = 10; %pruning threshold used in Pruning & merging
parameters.Tr = 0.001; %truncation threshold used in Pruning & merging 
load('birthpdf.mat');
models.birthpdf=birthpdf;

t_min = 0.053; %default value.
% t_min = 0.07; % current value: min track length in seconds. 
tl = t_min / dt;
% tl=10; %specify minimal track length in time steps - this is how long a 
%detection needs to be in order to become a true detection- tl=10 gives you
%tl*dt = 0.053 s (dt=time increment of your spectrogram). Using shorter tl
%will give you more false alarms, but potentially better recall (more
%detected sounds), and using higer tl will give you high precision (less
%false alrams) but consequently worse recall (less sounds that you expected
%to retrieve).
%% /////////////// RUN GMPHD DETECTOR ///////////////////
% ~~~~~~~~~~~~~~~~~~~~~~ Read audio data ~~~~~~~~~~~~~~~~~~~~~~~~~~~
[file, path] = uigetfile('*.wav','Select File to open');
[x,fs] = audioread(fullfile(path, file));
x=x(:,1); %select first channel
%Note: the following only handles one channel data. Adjust as necessary. 
secs_analisys = 1; % 4 seconds of analysis.
N_max = length(x);
N_analysis = round(secs_analisys * fs);
E_sp = [];
E_pic = [];
% Create detected events structure
idx = 1:N_analysis:N_max;


% Settings for the definition of the filter Bank of piknogram
BW=1000; % BW in Hz
BWoverlap=50; % BW in %
flow=5000;fhigh=min([50000 fs/2-BW/2]); % MINIMUM flow -> flow=round(BW/2)
Tl = 1;
%flow=2000;fhigh=min([25000 fs/2-BW/2]); % MINIMUM flow -> flow=round(BW/2)

for i = 1:length(idx)
%                i
%                 idx(i)
%                 idx(i) + N_analysis - 1
    if i==length(idx)
        y_part = x(idx(i):end); % Last chunk
    else
        y_part = x(idx(i) : idx(i) + N_analysis - 1); % Any chunk of N_analysis points
    end
    %% ~~~~~~~~~~~~~ Make MEASUREMENTS (Zsets) ~~~~~~~~~~~~~~~~~~~
    win_width = round(win_width_s * fs); 
    slide_incr = round(dt*fs);
    numstps = ceil((length(y_part) - win_width) / slide_incr); %Number of overlapping windows
    %% Piknogram implementation
    
    [FW,BW_est, ndraw] = pyknogram_freqdomain(y_part, fs, flow, fhigh, BW, BWoverlap, 10e-3);
    X=repmat(ndraw,1,size(FW,2)); 
    %  Get points with a Bandwidth lower than a given value (BW_thre)
    BW_thre=BW/4;
    BW_est(1,:)=BW; % To get rid of some problems in the BW estimation of the first and last
    BW_est(end,:)=BW; % To get rid of some problems in the BW estimation of the first and last
    idx_BW=find(BW_est<BW_thre);
    
    % ----------
    %  Point density based filter
    [P,idx_new]  = kernel_density( FW, ndraw, BW, BWoverlap );
    idx_combined = intersect(idx_BW,idx_new');   
    
    % Create a struct and store possible whistles
    Pyk=struct('time',cell(1,length(idx_combined)),'freq',[],'ampl',[],'label',[],'done',[]);
    kk=num2cell(X(idx_combined));[Pyk.time]=kk{:};
    kk=num2cell(FW(idx_combined));[Pyk.freq]=kk{:};
    
    % Represent candidates
%     figure(2);clf;
%     %sh=scatter([Pyk.time],[Pyk.freq],[],[Pyk.ampl],'filled');
%     sh=scatter([Pyk.time],[Pyk.freq],[],'b','filled');
%     sh.SizeData=15; 
%     title(['PRWE method']);
%     ylabel('Frequency [kHz]');xlabel('Time [sec.]');
%     yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     axis([0 Tl flow fhigh]);
%     set(gcf,'color','w');box on;grid;

    % Ordering candidates;
    tiempo = [Pyk.time];
    frecuencia = [Pyk.freq];
    [tiempo_ordenado, ind_ordenado] = sort(tiempo);
    frecuencia_ordenado = frecuencia(ind_ordenado);
    [tiempo_unico,ia,ic] = unique(tiempo_ordenado);
    Zset = cell([1 length(tiempo_unico)]);
    ic_unico = unique(ic);

    for ii = 1:length(Zset)
        Zset{ii} = sort(frecuencia_ordenado(ic == ic_unico(ii)));
    end
    
    % Represent ordered candidates
%     figure(2),clf
%     for j = 1:length(tiempo_unico)
%         plot(tiempo_unico(j),Zset{j},'b.'),hold on
%     end
%     title('Candidates Picknogram')
%     ylabel('Frequency [kHz]');xlabel('Time [sec.]');
%     yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     axis([0 Tl flow fhigh]);
%     set(gcf,'color','w');box on;grid;
% %     
% Conformation of all complet Zset candidates
    Zset_all = cell(1, numstps);
    t=(0:(size(Zset_all,2))).*dt;  
    
    for z = 1:length(Zset)
        [~,ix] = min(abs(t - tiempo_unico(z)));
        Zset_all{ix} = Zset{z};
    end
    Zset = Zset_all;
    % count the number of elements different of empty in the cell array
%     count = length(Zset_all) - nnz(cellfun(@isempty,Zset_all))
    %% Spectrogram implementation

    [Zset_sp] = preprocess_getZset(win_width_s,dt,fs,y_part,freqrange,peak_thr);
    %Zset is a cell array where each cell contains spectral peaks (frequencies)
    % from a particular time step that were above threshold (peak_thr - in dB).
% Represent Spectrogram candidates
%     figure(1),clf
%     t=(0:(size(Zset_sp,2))).*dt;
%     for m=1:size(Zset_sp,2)
%         if ~isempty(Zset_sp{m})
% %             plot(tiempo_unico(m),Zset{m},'k.'),hold on
%             plot(t(m),Zset_sp{m},'k.'),hold on % Uncomment to plot spectrogram candidates
%         end
%     end
%     title('Candidates Spectrogram')
%     ylabel('Frequency [kHz]');xlabel('Time [sec.]');
%     yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     axis([0 Tl flow fhigh]);
%     set(gcf,'color','w');box on;grid;
        
    %% ~~~~~~~~~~~~~~~~~ Run GMPHD detection ~~~~~~~~~~~~~~~~~~~
    
     [Xk_m_sp,XTag_sp] = gmphd_freqonly_adaptive(Zset_sp,parameters,models);
     [Xk_m,XTag] = gmphd_freqonly_adaptive(Zset,parameters,models);
    
    %% ~~~~~~~~~~~~~~~~ Extract Tracks (whistle contours) ~~~~~~~~~~~~~~~~~~
    [DT_sp,~,~] = tracktarget(XTag_sp,Xk_m_sp,dt,tl);
    [DT,~,~] = tracktarget(XTag,Xk_m,dt,tl);
    %DT is a structure of detected signals with 3 fields for each (freq x time x label)
    
    %save([folder,'GMPHD_',file(1:end-4),'.mat'],'DT')
    
    %% ~~~~~~~~~~~~~~~~~~ PLOT Detections ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %Plot detections against measurements Picnogram
    figure(3),clf
    t=(0:(size(Zset,2))).*dt + (i-1)*Tl;
    for m=1:size(Zset,2)
        if ~isempty(Zset{m})
%             plot(tiempo_unico(m),Zset{m},'k.'),hold on
            plot(t(m),Zset{m},'k.'),hold on % Uncomment to plot spectrogram candidates
        end
    end

    for m=1:size(DT,2)
	    DT(m).time = DT(m).time - 0.01 + idx(i)./fs;
        plot(DT(m).time, DT(m).freq,'LineWidth',1.5),hold on
    end

    title('Candidates vs detections PICKNOGRAM','Interpreter','latex')
    xlabel('time (s)','Interpreter','latex')
    ylabel('frequency (kHz)','Interpreter','latex')
    yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
    yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
    axis([min(t) max(t) flow fhigh]);
    set(gcf,'color','w');box on;grid;

    figure(4),clf
    t=(0:(size(Zset_sp,2))).*dt + (i-1)*Tl;
    for m=1:size(Zset,2)
        if ~isempty(Zset_sp{m})
%             plot(tiempo_unico(m),Zset{m},'k.'),hold on
            plot(t(m),Zset_sp{m},'k.'),hold on % Uncomment to plot spectrogram candidates
        end
    end

    for m=1:size(DT_sp,2)
	    DT_sp(m).time = (DT_sp(m).time * fs + idx(i))./fs;
        plot(DT_sp(m).time, DT_sp(m).freq,'LineWidth',1.5),hold on
    end

    title('Candidates vs detections SPECTROGRAM','Interpreter','latex')
    xlabel('time (s)','Interpreter','latex')
    ylabel('frequency (kHz)','Interpreter','latex')
    yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
    yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
    axis([min(t) max(t) flow fhigh]);
    set(gcf,'color','w');box on;grid;
    E_pic = [E_pic DT];
    E_sp = [E_sp DT_sp];
    DT = [];
    DT_sp = [];

    % Representation of the spectrogram
    figure(5),spectrogram(y_part,2048,[],freqrange(1):100:freqrange(2),fs,'yaxis'); %remember that the units are dB/Hz 
%     [~,f,t,ps_teorical] = spectrogram(y_part,2048,[],freqrange(1):100:freqrange(2),fs,'yaxis'); %remember that the units are dB/Hz 
%     figure(6), clf, pcolor(t,f,10*log10(abs(ps_teorical))),shading interp,colorbar
%     title('Spectrogram (psd) dB/Hz','Interpreter','latex')
%     xlabel('time (s)','Interpreter','latex')
%     ylabel('frequency (kHz)','Interpreter','latex')
%     yt=flow:2000:fhigh;set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt/1000);
%     axis([min(t) max(t) flow fhigh]);
%     set(gcf,'color','w');box on;grid;
end

