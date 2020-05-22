
close all
clear all
%% reading the data 
fq=25;
data1= xlsread('Anu - Copy.xlsx','Simple Data');
figure;
plot(data1(:,1), data1(:,2)-mean(data1(:,2))); 
xlabel('Time/s','fontsize', 14); ylabel('Signal magnitude', 'fontsize', 14); 
title('Raw Data from EMG - emg1', 'fontsize', 14);set(gca,'FontSize',14);

%% Step2: Perform Fast Fourier Transform and Shift Zero-component to center
% emg1
sEMG1 = data1(:,2); 
fft_sEMG1 = fft(sEMG1); 
fft_sEMG1 = fftshift(fft_sEMG1);
n = length(sEMG1);  
fq_axis1 = (-n/2:n/2-1)*fq/n; % zero-centered frequency range
abs_fft_sEMG1 = abs(fft_sEMG1);
figure; 
plot(fq_axis1, abs_fft_sEMG1 ); axis([fq_axis1(1) fq_axis1(end) 0 1000]); 
xlabel('Frequency/Hz','fontsize', 14); ylabel('X(jw) magnitude','fontsize', 14); 
title('Frequency Domain -emg1','fontsize', 14);
set(gca,'FontSize',14);

%% Step 3: Generate Bandpass filter 
% can change to your desired cutoff frequency
highpass = 5;   %%i am keeping this as 5 and 10 respectively
lowpass = 10;  
% get the index of -10 to -5 and 5 to 10Hz. 
cutoff1 = ceil((12.5-highpass)/(fq/length(sEMG1))); cutoff2 = ceil((12.5-lowpass)/(fq/length(sEMG1)));
cutoff3 = ceil((highpass+12.5)/(fq/length(sEMG1))); cutoff4 = ceil((lowpass+12.5)/(fq/length(sEMG1)));
sEMG1 = data1(:,2); 
H = zeros(length(sEMG1),1);
H(cutoff2:cutoff1) = 1; % take only the -10 to -5Hz
H(cutoff3:cutoff4) = 1; % take only the 5 to 10Hz
figure; 
plot(fq_axis1, H); 
set(gca,'YLim',[0 2]); 
xlabel('Freqeuncy/Hz','fontsize', 14); ylabel('Amplitude','fontsize', 14);
title('Bandpass filter','fontsize', 14);set(gca,'FontSize',14);

%% Step 4: Perform Bandpass filter, inverse shifting and inverse Fourier transform
% Biceps 1
cutoff1 = ceil((12.5-highpass)/(fq/length(sEMG1))); cutoff2 = ceil((12.5-lowpass)/(fq/length(sEMG1)));
cutoff3 = ceil((highpass+12.5)/(fq/length(sEMG1))); cutoff4 = ceil((lowpass+12.5)/(fq/length(sEMG1)));
H = zeros(length(sEMG1),1); H(cutoff2:cutoff1) = 1; % take only the -10 to -5Hz
H(cutoff3:cutoff4) = 1; % take only the 5 to 10Hz

yt1 = ifftshift(fft_sEMG1.*H); 
yt1 = ifft(yt1);
t = 1/fq*(1:length(yt1));
figure;
plot(t, real(yt1),'r');xlabel('Time/s','fontsize', 14); ylabel('Amplitude','fontsize', 14); 
title('Signal after bandpass filter - emg 1','fontsize', 14);set(gca,'FontSize',14);

%% Step 5: Feature extraction - MAV using convolution. 
dt = 2*fq; % hand flex for 2 sec; 
dt_5 = 5*fq;   % subject flex at 5s, 10s,15s, 20s, 25s, 30s.
wind = 2*fq;
filter = ones(20,1);
filter = filter/length(filter);

% Biceps 1 
MAV1 = zeros(length(yt1),1);
var1 = zeros(length(yt1),1);
MAV1 = conv(abs(real(yt1)),filter,'same');

%% Step 6: K-mean clustering
n = 10; %clustering window

for i = 1:length(MAV1)-n
    X1(i,:) =  MAV1(i:i+n);
end
[idx1,C1] = kmeans(X1,2);

% theorectical response and adjusting
dt_r = 0.04;
t_response = [0:dt_r:35]';
response = zeros(size(t_response));
START = 5/dt_r;
dt_r2 = 2/dt_r;
dt_r5 = 5/dt_r;
response(START+1:START+dt_r2,1) = 1;
response(2*START+1:2*START+dt_r2,1) = 1;
response(3*START+1:3*START+dt_r2,1) = 1;
response(4*START+1:4*START+dt_r2,1) = 1;
response(5*START+1:5*START+dt_r2,1) = 1;
response(6*START+1:6*START+dt_r2,1) = 1;

figure;
plot(data1(:,1),MAV1,'linewidth',2); hold on; 
plot(data1(:,1),real(yt1),'r'); 
plot(data1(:,1),[zeros(ceil(n/2),1);idx1;zeros(floor(n/2),1)],'k', 'linewidth', 2); 
plot(t_response(:,1),response,'g--','linewidth',3);
legend('MAV','EMG signal','Contract(Min); Relax (Max)', 'True response', 'location', 'Northeast');
xlabel('Time/s'); ylabel('Signal Amplitude'); title('emg 1'); grid on;
axis([0 round(max(data1(:,1))) min(real(yt1)) max(real(yt1))]);
set(gca,'FontSize',14);


num=[0 1];
den=[3.51 1.3 0];
fs=1/200;
syg=tf(num,den,fs);
syg;
[y,t]=lsim(syg,yt1);
figure;
stem(t,y);
ylabel('Amplitude of signal');
xlabel('Time');
title('System Plot');

[c,lags] = xcorr(yt1,y);
r=xcorr(yt1,y);
r;
figure;
stem(lags,c);

r1=corrcoef([yt1,y]);
r1;
c2=xcorr(yt1,y,'coeff');
c2;