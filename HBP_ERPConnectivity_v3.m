%%
clear all
close all
clc

% NB: Must have the same parameters of convolution, dB normalization 
%     and ITPC of Frequency_Analysis

%% LOADING

load j23_r4_Keta2_PREPROCESSED_ok;

clear eeg_good eeg_diff sw_start sw_end ch_minus ch_plus


%% WAVELET CONVOLUTION 

% ISPC - INTER SITE PHASE CLUSTERING:
%
% CREATE A FAMILY OF WAVELETS:
%

wavelet_time = -2:1/fs:2; 
n = 3;                    % number of cycles 

lowest_frequency  =   1;  % lowest pk [Hz]
highest_frequency =  14;  % highest pk [Hz]
num_wavelets      =  (highest_frequency-lowest_frequency)+1;  % number of pk (frequency bands) 

% create a vector of peak frequencies:
frequencies = linspace(lowest_frequency,highest_frequency,num_wavelets);

% initialize the matrix of wavelet family:
wavelet_family = zeros(num_wavelets,length(wavelet_time));
 
% Loop through frequencies and make a family of wavelets.
for fi=1:num_wavelets
    % create a sine wave at this frequency
    sinewave = exp(2*1i*pi*frequencies(fi).*wavelet_time); % the "1i" makes it a complex wavelet
    % create a Gaussian window
    gaus_win = exp(-wavelet_time.^2./(2*(n/(2*pi*frequencies(fi)))^2));
    % create wavelet via element-by-element multiplication of the sinewave and gaussian window
    wavelet_family(fi,:) = sinewave.*gaus_win;
end
clear fi


%
% CONVOLUTION 
%

prestim = abs(time(:,1)); % in [s], define the max abs time pre stimulus
X1= -1; % in [s], define the lower time for TIME FREQUENCY ANALYSIS
X2= 1; % in [s], define the upper time for TIME FREQUENCY ANALYSIS

% New time window for Burst Suppression Sudy (ONLY!)
%X1= -2; % in [s], define the lower time for TIME FREQUENCY ANALYSIS
%X2= 1; % in [s], define the upper time for TIME FREQUENCY ANALYSIS

% define subset data for TIME FREQUENCY ANALYSIS
data = eeg (((fs.*X1)+(fs.*prestim)+1):((fs.*X2)+(fs.*prestim)+1),:,:);

% time_data
time_data = X1:1/fs:X2;

% CONVOLUTION: loop through frequencies and compute synchronization
CONV = zeros (size(wavelet_family,1),size(data,1),size(data,3),size(data,2));

for fr = 1:size(wavelet_family,1);
for sw = 1:size(data,3);
for ch = 1:size(data,2);
    
    CONV(fr,:,sw,ch) = conv (squeeze(data(:,ch,sw)),squeeze(wavelet_family(fr,:)), 'same'); 

end 
end
end

clear fr sw ch
clear gaus_win
clear highest_frequency
clear lowest_frequency
clear n 
clear num_wavelets
clear prestim
clear sinewave
clear wavelet_family
clear wavelet_time


% EXTRACTION OF CONVOLUTION PRODUCTS

%
% EXTRACT PHASES AND APPLY EULER's FORMULA:
%

phases = angle(CONV(:,:,:,:));    % exctract phases from convolution

% Compute the average vector of phases of each trials at all time points and peak frequencies
% by Euler's formula (vector lengths set to 1)
% NB: av_vector_phase = average phase among trials [radiants]
%     av_vector_length= average phase locking among trials [between 0 and 1]

av_vector_phase = squeeze (angle(mean(exp(1i.*(phases(:,:,:,:))),3))); 
av_vector_length = squeeze (abs(mean(exp(1i.*(phases(:,:,:,:))),3))); 


clear CONV

% BASELINE THRESHOLD OF ISPC (Bootstrap)

% Define data of Baseline

% baseline window
p1 = -.5;    %[s]
p2 = -.2;

% Take the phases of baseline
Phases_baseline = phases(:,((fs.*p1)+(fs.*abs(X1))+1):((fs.*p2)+(fs.*abs(X1))+1),:,:);


% Connectome of Phase differences of baseline
connect_phasediff_baseline = zeros(size (phases,1),size (Phases_baseline,2),size (phases,3),size (phases,4),size (phases,4));

for chii = 1:size(phases,4)
for chjj = 1:size(phases,4)
    
    connect_phasediff_baseline(:,:,:,chii,chjj) = Phases_baseline(:,:,:,chii) - Phases_baseline(:,:,:,chjj);

end
end

clear chii chjj sw fr sp
% clear Phases_baseline 


% Euler's Formula:
Euler_baseline = (exp(1i*(connect_phasediff_baseline(:,:,:,:,:))));

% Real Baseline:
ISPC_baseline = squeeze(abs(mean(Euler_baseline(:,:,:,:,:),3)));
ISPC_baseline = squeeze(mean(ISPC_baseline,2));

% ISPC THRESHOLD (Boostrap Thresholding)

% BOOSTRAP CORRECTED (about 10-12min for 200 straps)

% Define Boostrap Parameters:

baseline = Euler_baseline; % matrice da cui si calcola la distribuzione (sample X sweep)
n_straps = 500; % total number of resamplings (bootstraps; usually between 500 and 1000 is enough)
alpha = 0.05;   % significant level (probability of rejecting the null hypothesis when it's true) can be 0,05 or less
n_samples = size (baseline,2); % number of samples
n_sweeps = size (baseline,3);  % number of sweeps
n_freq = size(baseline,1);     % number of frequencies
n_chii = size(baseline,4);  % number of channels
n_chjj = size(baseline,5);  % number of channels

% loop over bootstrap samples

bootstrap_max = zeros(n_freq,n_chii,n_chjj,n_straps);
bootstrap_min = zeros(n_freq,n_chii,n_chjj,n_straps);

for bb =   1:n_straps

    
        % Take randomly (n = n_sweeps) samples from data (indexes of sweeps, 1:nsweeps), 
        % with replacement:
        resampled_sweeps = datasample(1:n_sweeps,n_sweeps); 
        
        % Take the actual sweeps corresponding to the randomized indexes
        % from the baseline matrix, for 1 channel (ch) and 1 frequency (fr) at the time:
        %resampled_baseline = baseline(fr,:,resampled_sweeps,chii,chjj);
        resampled_baseline = baseline(:,:,resampled_sweeps,:,:);
        
        % Make the ISPC  of the new resempled baseline matrix:
        
   %    % Extract Clustering  (resampled_ISPC_baseline: vector n of points x 1)
        resampled_ISPC_baseline = squeeze(abs(mean(resampled_baseline,3)));
   %    % Baseline Correction (connect_baseline:  single value)
        connect_baseline = squeeze(mean(resampled_ISPC_baseline,2));
  
   %    % Baseline Corrected (no offset) (connect_ISPCnorm_baseline: vector n of points x1)
        connect_ISPCnorm_baseline = zeros(n_freq, n_samples, n_chii, n_chjj);
  for fr =   1:n_freq
  for chii = 1:n_chii
  for chjj = 1:n_chjj        
        
   %    % Baseline Corrected (no offset) (connect_ISPCnorm_baseline: vector n of points x1)
        connect_ISPCnorm_baseline(fr,:,chii,chjj) = resampled_ISPC_baseline(fr,:,chii,chjj) - connect_baseline(fr,chii,chjj);
  end
  end
  end
  clear fr chii chjj
  
        
        % Find the maximum and the minimum value of all the new ISPC of
        % running baseline and create a matrix with MAX and MIN for each
        % frequecy and channel x channel:
        bootstrap_max(:,:,:,bb) = max(connect_ISPCnorm_baseline,[],2);
        bootstrap_min(:,:,:,bb) = min(connect_ISPCnorm_baseline,[],2);
end
clear bb

% Find thresholds

max_percentile = zeros(n_freq, n_chii, n_chjj);
min_percentile = zeros(n_freq, n_chii, n_chjj);

for chii = 1:n_chii
for chjj = 1:n_chjj
for fr = 1:n_freq
    
    max_percentile(fr,chii,chjj) = prctile(bootstrap_max(fr,chii,chjj,:),100.*(1-alpha));
    min_percentile(fr,chii,chjj) = prctile(bootstrap_min(fr,chii,chjj,:),100.*alpha);

end
end
end
clear ii chii chjj fr

% 

clear baseline Euler_baseline n_straps alpha n_samples n_sweeps n_freq n_chii n_chjj
clear bootstrap_max bootstrap_min resampled_sweeps resampled_baseline resampled_ISPC_baseline
clear resampled_ISPC_baseline connect_baseline connect_ISPCnorm_baseline connect_phasediff_baseline




% INTER SITE PHASE CLUSTERING-trial (ISPCtr) across all channels 


% Define data 

% signal window
w1= -0.1; % in [s], define the time window for mean ISPCtrials (valore usato x paper )
w2= 0.4; 

% Take the phases of signal of interest
Phases_signal = phases(:,((fs.*w1)+(fs.*abs(X1))+1):((fs.*w2)+(fs.*abs(X1))+1),:,:);
% Take the signal time
time_signal = time_data(:,((fs.*w1)+(fs.*abs(X1))+1):((fs.*w2)+(fs.*abs(X1))+1));



% Connectome of Phase differences of signal
connect_phasediff_signal = zeros(size (phases,1),size (Phases_signal,2),size (phases,3),size (phases,4),size (phases,4));

for chii = 1:size(phases,4)
for chjj = 1:size(phases,4)
    
    connect_phasediff_signal(:,:,:,chii,chjj) = Phases_signal(:,:,:,chii) - Phases_signal(:,:,:,chjj);

end
end

clear chii chjj sw fr sp
clear Phases_baseline Phases_signal


% On Signal:
% Extract Clustering
ISPC_signal_raw = squeeze(abs(mean(exp(1i*(connect_phasediff_signal(:,:,:,:,:))),3)));
 % Extract Preferred angle of the phase difference
angle_signal = squeeze(angle(mean(exp(1i*(connect_phasediff_signal(:,:,:,:,:))),3)));

clear connect_phasediff_signal

% Baseline Correction:

ISPC_signal= zeros(size(ISPC_signal_raw,1), size(ISPC_signal_raw,2),size(ISPC_signal_raw,3),size(ISPC_signal_raw,4));
  for fr =   1:size(ISPC_signal_raw,1)
  for chii = 1:size(ISPC_signal_raw,3)
  for chjj = 1:size(ISPC_signal_raw,4)        
        
        ISPC_signal(fr,:,chii,chjj) = ISPC_signal_raw(fr,:,chii,chjj) - ISPC_baseline(fr,chii,chjj);
  end
  end
  end
  clear fr chii chjj
  clear ISPC_signal_raw
  


% Thresholding
ISPC_signalsig = zeros (size(ISPC_signal,1), size(ISPC_signal,2), size(ISPC_signal,3), size(ISPC_signal,4));

for fr   = 1:size(ISPC_signal,1)
for sp   = 1:size(ISPC_signal,2)
for chii = 1:size(ISPC_signal,3) 
for chjj = 1:size(ISPC_signal,4)
    
    if ISPC_signal(fr,sp,chii, chjj) > max_percentile(fr,chii,chjj)
        ISPC_signalsig (fr,sp,chii,chjj) = ISPC_signal(fr,sp,chii, chjj);
    else
        
    if ISPC_signal(fr,sp,chii, chjj) < min_percentile(fr,chii,chjj)
        ISPC_signalsig (fr,sp,chii,chjj) = ISPC_signal(fr,sp,chii, chjj);
    else
        
        ISPC_signalsig (fr,sp,chii,chjj) = 0;
    end
    end
    
end
end
end
end
clear fr chii chjj sp
clear ISPC_signal


% TEST FOR VOLUME CONDUCTION:
 
 % Evaluates if the preferred angle of the clustering (mean phase angle for ITPC or mean phase angle difference for ISPC)
 % is significantly different from another specific phase angle.
 % Since volume conduction is istantaneous, it will generate a ISPC
 % at the preferred phase angle difference 0 or Pi greco (depending on wich side of the dipole) 
 % --> the preferred angle difference is compared to 0 or Pi 
 
 % NB1: all values must be in radiants
 % NB2: H0 (null hypothesis) = the phase angle difference is different from 0
 %      -> if p<0,05 -> H0 is rejected -> the phase angle difference = 0, therefore explained by volume conduction
 %      -> if p>0,05 -> H0 can't be rejected -> no evidence for volume conduction
 % NB3: the statistic of the test (u) is distributed as a gauss under H0, with mean of about 0 and varianve of about 1
 %      -> u = 1.96 correspond to p=0,05 (two-tailed)
 %      -> if u < 1.96 and > -1.96 -> H0 can't be rejected -> no evidence for volume conduction 
 %

 % Original test (Durand&Greenwood,1958; Zar,1999):
 % vtest  = @(icpcmag,n,val) n.*icpcmag*cos(val).*sqrt(2./n);
 
 % Modified by Choen: 
 % resistant to high number of samples and produces less false positive 
 % the original cosine component is sobstituted by a gaussian component
 % gvtest = @(icpcmag,n,val) n.*(icpcmag*exp((-(val).^2)./(4.*pi./n)).*(sqrt(2./n)));
 
 % Variables:
 % n = number of samples or trials 
 % icpcmag = clustering value (ITPC or ISPC)
 % val = phase angle difference (preferred angle of ISPC - 0 or Pi for tesing volume conduction)


angle_signal_s = angle_signal;

% Make gv test
connect_gvtest_0  = zeros(size(angle_signal,1), size(angle_signal,2), size(angle_signal,3), size(angle_signal,4));
connect_gvtest_pi = zeros(size(angle_signal,1), size(angle_signal,2), size(angle_signal,3), size(angle_signal,4));
connect_gvtest_pin = zeros(size(angle_signal,1), size(angle_signal,2), size(angle_signal,3), size(angle_signal,4));

for fr = 1:size(angle_signal,1)
for sp = 1:size(angle_signal,2)
for chii = 1:size(angle_signal,3)
for chjj = 1:size(angle_signal,4)
    
    connect_gvtest_0(fr, sp, chii, chjj) =  size(data,3).*(ISPC_signalsig(fr, sp, chii, chjj).*exp((-(angle_signal_s(fr, sp, chii, chjj)-0).^2)./(4.*pi./size(data,3))).*(sqrt(2./size(data,3))));
    connect_gvtest_pi(fr, sp, chii, chjj)=  size(data,3).*(ISPC_signalsig(fr, sp, chii, chjj).*exp((-(angle_signal_s(fr, sp, chii, chjj)-pi).^2)./(4.*pi./size(data,3))).*(sqrt(2./size(data,3))));
    connect_gvtest_pin(fr, sp, chii, chjj)= size(data,3).*(ISPC_signalsig(fr, sp, chii, chjj).*exp((-((angle_signal_s(fr, sp, chii, chjj).*-1)-pi).^2)./(4.*pi./size(data,3))).*(sqrt(2./size(data,3))));
   
end 
end 
end
end
clear fr sp chii chjj


% Obtain p value matrix and select only the not significant points
pconnect_gvtest_0 = 1-normcdf(connect_gvtest_0);
pconnect_gvtest_pi = 1-normcdf(connect_gvtest_pi);
pconnect_gvtest_pin = 1-normcdf(connect_gvtest_pin);

ISPC_signalf = zeros (size(ISPC_signalsig,1), size(ISPC_signalsig,2), size(ISPC_signalsig,3), size(ISPC_signalsig,4));

for fr = 1:size(ISPC_signalsig,1)
for sp = 1:size(ISPC_signalsig,2)
for chii = 1:size(ISPC_signalsig,3) 
for chjj = 1:size(ISPC_signalsig,4)
    
    if pconnect_gvtest_0(fr,sp,chii,chjj) < 0.05
    ISPC_signalf (fr,sp,chii, chjj) = 0;
    else

    if pconnect_gvtest_pi(fr,sp,chii,chjj) < 0.05
    ISPC_signalf (fr,sp,chii, chjj) = 0;
    else      
        
    if pconnect_gvtest_pin(fr,sp,chii,chjj) < 0.05
    ISPC_signalf (fr,sp,chii, chjj) = 0;
    else    
        
    ISPC_signalf (fr,sp,chii, chjj) = ISPC_signalsig (fr,sp,chii, chjj);
    
    end
    end
    end
end
end  
end
end

clear fr chii chjj sp 
clear connect_gvtest_0 connect_gvtest_pi pconnect_gvtest_0 pconnect_gvtest_pi pconnect_gvtest_pin pconnect_gvtest_0n
clear phases max_percentile min_percentile connect_gvtest_0n connect_gvtest_pin

%%
% PLOT exploratory
% Select and mean in the window of interest:
% Defie the window:

% Time window:
tw1 = 0.18; % [s]
tw2 = 0.4;
twx1 = find(round(time_signal,4) == tw1);
twx2 = find(round(time_signal,4) == tw2); 

% High Frequency Range:
fH_start = 5; % [Hz] bottom limit
fH_x1 = find(frequencies == fH_start);
fH_end = 14;  % [Hz] upper limit
fH_x2 = find(frequencies == fH_end);
% Low Frequency Range:
fL_start = 1; % [Hz] bottom limit
fL_x1 = find(frequencies == fL_start);
fL_end = 14;  % [Hz] upper limit
fL_x2 = find(frequencies == fL_end);

CONNECT = ISPC_signalf;   % ISPC without volume conduction
%CONNECT = ISPC_signalsig; % ISPC with volume conduction

connect_ISPC_t1 = squeeze (mean(CONNECT(:,twx1:twx2,:,:),2));

connect_ISPC_t1_fH = squeeze (mean(connect_ISPC_t1(fH_x1:fH_x2,:,:),1));
connect_ISPC_t1_fL = squeeze (mean(connect_ISPC_t1(fL_x1:fL_x2,:,:),1));
 

x1 = -0.1;
x2 = .4;
y1 = -100;
y2 = 200;

s1= -.3;
s2= .3;

f1=5;
f2=14;

cmap = 'parula';
c1 = 0;
c2 = 0.1;

seed1 = 9;
ch1 = 14;

seed2 =12;
ch2 =15;


clf
 % Plot voltage signal
figure (1)
subplot(3,2,1)
plot(time_data,squeeze(mean(data(:,ch1,:),3)),'k')
hold on
plot(time_data,squeeze(mean(data(:,seed1,:),3)),'b')
xlabel('Time [s]'), ylabel('Amplitude [uV]')
axis([x1,x2, y1,y2])
title(['sweep; ch' num2str(seed1) 'vs' num2str(ch1) ''])

subplot(3,2,2)
plot(time_data,squeeze(mean(data(:,ch2,:),3)),'k')
hold on
plot(time_data,squeeze(mean(data(:,seed2,:),3)),'b')
xlabel('Time [s]'), ylabel('Amplitude [uV]')
axis([x1,x2, y1,y2])
title(['sweep; ch' num2str(seed2) 'vs' num2str(ch2) ''])
%

subplot(3,2,3)
contourf(time_signal,frequencies,CONNECT(:,:,seed1,ch1),100,'linecolor','none')
colormap (cmap)
axis([x1,x2, f1,f2])
set(gca,'clim',[s1 s2])
xlabel('Time [s]'), ylabel('Frequency [Hz]')
title(['ISPC; ch' num2str(seed1) 'vs' num2str(ch1) ''])

subplot(3,2,4)
contourf(time_signal,frequencies,CONNECT(:,:,seed2,ch2),100,'linecolor','none')
colormap (cmap)
axis([x1,x2, f1,f2])
set(gca,'clim',[s1 s2])
xlabel('Time [s]'), ylabel('Frequency [Hz]')
title(['ISPC; ch' num2str(seed2) 'vs' num2str(ch2) ''])
%

subplot(3,2,5)
imagesc(connect_ISPC_t1_fH)
colormap(cmap)
set(gca,'clim',[c1 c2])
xlabel('Channels'), ylabel('Channels')
title(['ISPC ' num2str(fH_start) '-' num2str(fH_end) 'Hz ' num2str(tw1) '-' num2str(tw2) 's'])


subplot(3,2,6)
imagesc(connect_ISPC_t1_fL)
colormap(cmap)
set(gca,'clim',[c1 c2])
xlabel('Channels'), ylabel('Channels')
title(['ISPC ' num2str(fL_start) '-' num2str(fL_end) 'Hz ' num2str(tw1) '-' num2str(tw2) 's'])



clear x1 x2 y1 y2 z1 z2 f1 f2 ch fr1 fr2 cmap ch1 ch2 seed1 seed2 
clear connect_ISPC_t1_fH connect_ISPC_t1_fL connect_ISPC_t1 CONNECT
clear connect_ISPC_t1 twx1 twx2 tw1 tw2 w1 w2 X1 X2
clear fL_x2 fL_end fL_x1 fL_start fH_x2 fH_end fH_x1 fH_start 




% ANALYSIS

% DEFINE MATRIX

% Define time and frequency window

% Time window:
tw1 = 0.18; % [s]
tw2 = 0.40;
twx1 = find(round(time_signal,4) == tw1);
twx2 = find(round(time_signal,4) == tw2); 

% High Frequency Range: (or theta-alfa)
fH_start = 5; % [Hz] bottom limit
fH_x1 = find(frequencies == fH_start);
fH_end = 14;  % [Hz] upper limit
fH_x2 = find(frequencies == fH_end);

% Take the connectivity matrix and print ranges:
connect_ISPC_t1 = squeeze (mean(ISPC_signalf(:,twx1:twx2,:,:),2));

connect_ISPC_FH = squeeze (mean(connect_ISPC_t1(fH_x1:fH_x2,:,:),1));

range_T  = [tw1,tw2];
range_FH = [fH_start,fH_end];

clear connect_ISPC_t1 twx1 twx2 tw1 tw2
clear fH_x2 fH_end fH_x1 fH_start 
clear connect_ISPC_t1_fH connect_ISPC_t1_fL


% CONNECTIVITY DEGREE:
% Number of supra or below threshold connections for each electrode

% Convert all positive values in 1
positive_ISPC_FH = connect_ISPC_FH > 0;

% % Convert all negative values in 1
% negative_ISPC_FH = connect_ISPC_FH < 0;

% Positive Connectivity Degrees in single channels 
% normalized for the total possible connections:
 ConnDeg_FH = (sum(positive_ISPC_FH, 2))/(size(positive_ISPC_FH,1)-1);

% Positive Connectivity Degrees in global matrix 
% normalized for the total possible connections:
 ConnDegGl_FH = (sum((sum(positive_ISPC_FH, 2)), 1))/((size(positive_ISPC_FH,1).^2)-(size(positive_ISPC_FH,1)));


 
% CLUSTERING COEFFICIENT:

% initialize for high frequencies
clustcoef_FH = zeros(size(data,2),1);
        
for chani = 1:size(data,2)
    % find neighbors (suprathreshold connections with electrode chani)
    neighbors = find(positive_ISPC_FH(chani,:));
    n = length(neighbors); % count how many they are (connectivity degree))
            
    % cluster coefficient not computed for islands
    % (consider only the electrodes that have more then 1
    % suprashold connection, all the other mantein the clustcoef = 0 of the initialized matrix)
    if n>1
        % "local" network of neighbors
        % create a new matrix with the connections among only the
        % suprashold connections of electrode chani
        localnetwork = positive_ISPC_FH(neighbors,neighbors);
        % localnetwork is symmetric; remove redundant values by
        % replacing with NaN (tril function: take only the upper part of the matrix)
        localnetwork = localnetwork + tril(nan(n));
                
        % compute cluster coefficient (neighbor connectivity scaled)
        clustcoef_FH(chani) = 2*nansum(localnetwork(:)) / ((n-1)*n);
    end
end

clear chani localnetwork neighbors n


clustcoefGl_FH = mean(clustcoef_FH,1);


clear positive_ISPC_FH


%%
% PLOT exploratory
% Select and mean in the window of interest:
% Defie the window:

% Time window:
tw1 = 0.18; % [s]
tw2 = 0.4;
twx1 = find(round(time_signal,4) == tw1);
twx2 = find(round(time_signal,4) == tw2); 

% High Frequency Range:
fH_start = 5; % [Hz] bottom limit
fH_x1 = find(frequencies == fH_start);
fH_end = 14;  % [Hz] upper limit
fH_x2 = find(frequencies == fH_end);
% Low Frequency Range:
fL_start = 1; % [Hz] bottom limit
fL_x1 = find(frequencies == fL_start);
fL_end = 14;  % [Hz] upper limit
fL_x2 = find(frequencies == fL_end);

CONNECT = ISPC_signalf;   % ISPC without volume conduction
%CONNECT = ISPC_signalsig; % ISPC with volume conduction

connect_ISPC_t1 = squeeze (mean(CONNECT(:,twx1:twx2,:,:),2));

connect_ISPC_t1_fH = squeeze (mean(connect_ISPC_t1(fH_x1:fH_x2,:,:),1));
connect_ISPC_t1_fL = squeeze (mean(connect_ISPC_t1(fL_x1:fL_x2,:,:),1));
 

x1 = -0.05;
x2 = .4;
y1 = -100;
y2 = 150;

s1= -.3;
s2= .3;

f1=5;
f2=20;

cmap = 'parula';
c1 = 0;
c2 = 0.1;

 seed1 = 10;
 ch1 = 13;
%seed1 = 3;
%ch1 = 13;

seed2 =9;
ch2 =14;
sw=60;

clf
 % Plot voltage signal
figure (1)
subplot(3,2,1)
plot(time_data,squeeze(mean(data(:,ch1,:),3)),'LineWidth',1,'Color',[0 0 0])
hold on
plot(time_data,squeeze(mean(data(:,seed1,:),3)),'LineWidth',1,'Color',[0.5 0.5 0.5])
xlabel('Time [s]'), ylabel('Amplitude [uV]')
axis([x1,x2, y1,y2])
title(['sweep; ch' num2str(seed1) 'vs' num2str(ch1) ''])
legend ('V2','PA')

subplot(3,2,2)
% plot(time_data,squeeze(mean(data(:,ch2,:),3)),'LineWidth',1,'Color',[0 0.1 0.8])
% hold on
% plot(time_data,squeeze(mean(data(:,seed2,:),3)),'LineWidth',1,'Color',[0 0.5 0.8])

plot(time_data,squeeze(mean(data(:,ch2,:),3)),'LineWidth',1,'Color',[0 0.5 0.2])
hold on
plot(time_data,squeeze(mean(data(:,seed2,:),3)),'LineWidth',1,'Color',[0 0.7 0.4])

xlabel('Time [s]'), ylabel('Amplitude [uV]')
axis([x1,x2, y1,y2])
title(['sweep; ch' num2str(seed2) 'vs' num2str(ch2) ''])
legend ('V2','PA')
%
% % Sevoflurane colour    
% subplot (4,2,2);
%     plot (time_data, MEDIA2 (:, :), 'LineWidth',0.5,'Color',[0 0.5 0.8])
%     axis([x1,x2, y1,y2])
%     xlabel ('(s)')
%     ylabel ('(uV)')
%     title  ('Mean ERPs for all channels')
% hold on
%     plot (time_data, MEDIA2 (:, ch), 'LineWidth',1.5,'Color',[0 0.1 0.8])
%     axis([x1,x2, y1,y2])

%Ketamine colour    
% subplot (4,2,2);
%     plot (time_data, MEDIA2 (:, :), 'LineWidth',0.5,'Color',[0 0.7 0.4])
%     axis([x1,x2, y1,y2])
%     xlabel ('(s)')
%     ylabel ('(uV)')
%     title  ('Mean ERPs for all channels')
% hold on
%     plot (time_data, MEDIA2 (:, ch), 'LineWidth',1.5,'Color',[0 0.5 0.3])
%     axis([x1,x2, y1,y2])




subplot(3,2,3)
contourf(time_signal,frequencies,CONNECT(:,:,seed1,ch1),100,'linecolor','none')
colormap (cmap)
axis([x1,x2, f1,f2])
set(gca,'clim',[s1 s2])
xlabel('Time [s]'), ylabel('Frequency [Hz]')
title(['ISPC; ch' num2str(seed1) 'vs' num2str(ch1) ''])

subplot(3,2,4)
contourf(time_signal,frequencies,CONNECT(:,:,seed2,ch2),100,'linecolor','none')
colormap (cmap)
axis([x1,x2, f1,f2])
set(gca,'clim',[s1 s2])
xlabel('Time [s]'), ylabel('Frequency [Hz]')
title(['ISPC; ch' num2str(seed2) 'vs' num2str(ch2) ''])
%

subplot(3,2,5)
imagesc(connect_ISPC_t1_fH)
colormap(cmap)
set(gca,'clim',[c1 c2])
xlabel('Channels'), ylabel('Channels')
title(['ISPC ' num2str(fH_start) '-' num2str(fH_end) 'Hz ' num2str(tw1) '-' num2str(tw2) 's'])


subplot(3,2,6)
imagesc(connect_ISPC_t1_fL)
colormap(cmap)
set(gca,'clim',[c1 c2])
xlabel('Channels'), ylabel('Channels')
title(['ISPC ' num2str(fL_start) '-' num2str(fL_end) 'Hz ' num2str(tw1) '-' num2str(tw2) 's'])



clear x1 x2 y1 y2 z1 z2 f1 f2 ch fr1 fr2 cmap ch1 ch2 seed1 seed2 
clear connect_ISPC_t1_fH connect_ISPC_t1_fL connect_ISPC_t1 CONNECT
clear connect_ISPC_t1 twx1 twx2 tw1 tw2 w1 w2 X1 X2
clear fL_x2 fL_end fL_x1 fL_start fH_x2 fH_end fH_x1 fH_start 


 
%% Saving 

save ('j23_r4_Keta2_CONNECTIVITY_HBP_ok') 
