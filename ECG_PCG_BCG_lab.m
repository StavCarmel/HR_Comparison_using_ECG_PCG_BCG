clear all; clc; close all;

%% Parameters

cut=1; % turn on (cut=1) if cutting the signal is needed 
%the seconds to cut from, change if needed
cut_rest_start=5; %[sec]
cut_rest_end=35; %[sec]
cut_exe_start=9.5; %[sec]
cut_exe_end=36; %[sec]

%change the channels according to the specific signals's channels (and add 1 because the first column is time)
channel_BCG=2;
channel_ECG=3;
channel_PCG=4;

Fs = 1000;          % [samples/sec]     - estimated sample rate according to data screen shot [Hz

rest_file_name='rest.txt';
exercise_file_name='exercise.txt';

%% Reading data
%change the text file name if needed
rest=load(rest_file_name);
exercise=load(exercise_file_name);
 
if cut==1
    %cut the start and end:
    rest=rest(cut_rest_start*Fs:cut_rest_end*Fs,:);
    exercise=exercise(cut_exe_start*Fs:cut_exe_end*Fs,:);
end

%%  1. ECG signal                                        

%% reading the data and defiyng the parameters

% taking the data from the current channel
raw_rest_ECG=rest(:,channel_ECG);
raw_exe_ECG=exercise(:,channel_ECG);

% making the time vectors
dt = 1/Fs; % [sec]
N_rest=length(raw_rest_ECG);  % [#] - number of samples
N_exe=length(raw_exe_ECG); % [#] - number of samples
t_rest=(0:N_rest-1)*dt; % [sec] - resting time vector
t_exe=(0:N_exe-1)*dt; % [sec] - active time vector

% plot the raw signals
%rest
figure;
subplot(3,1,1);
plot(t_rest,rest(:,channel_ECG),'b');
title({'Signals-rest','ECG'});
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_rest(end)]);
subplot(3,1,2);
plot(t_rest,rest(:,channel_PCG),'b');
title('PCG');
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_rest(end)]);
subplot(3,1,3);
plot(t_rest,rest(:,channel_BCG),'b');
title('BCG');
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_rest(end)]);

%after exercise
figure;
subplot(3,1,1);
plot(t_exe,exercise(:,channel_ECG),'b');
title({'Signals- after exercise','ECG'});
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_exe(end)]);
subplot(3,1,2);
plot(t_exe,exercise(:,channel_PCG),'b');
title('PCG');
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_exe(end)]);
subplot(3,1,3);
plot(t_exe,exercise(:,channel_BCG),'b');
title('BCG');
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_exe(end)]);

%% Filtering + R peaks detection

%filtering the ECG signal
filt_rest_ECG=filter_the_signal(raw_rest_ECG,Fs);
filt_exe_ECG=filter_the_signal(raw_exe_ECG,Fs);

%R detection
t_peak_rest_ECG=find_R_func(1,raw_rest_ECG,Fs);
t_peak_exe_ECG=find_R_func(1,raw_exe_ECG,Fs);


%% plot the signals with the R peaks

%plots unfiltered and filtered rest signal
figure;
subplot(2,1,1);
plot(t_rest,raw_rest_ECG,'b');
xlim([0 t_rest(end)]);
xlabel('time [sec]');
ylabel('ECG signal [mv]');
title('Unfiltered ECG - rest');
subplot(2,1,2);
plot(t_rest,filt_rest_ECG,'b');
hold on;
scatter(t_rest(t_peak_rest_ECG),filt_rest_ECG(t_peak_rest_ECG),10,'filled','r');
xlim([0 t_rest(end)]);
ylim([min(filt_rest_ECG)-0.05 max(filt_rest_ECG)+0.1])
xlabel('time [sec]');
ylabel('ECG signal [mv]');
title('Peak detection - rest');
legend('Filtered ECG','R peaks')
hold off;

%plot unfiltered and filtered after exercise signal
figure;
subplot(2,1,1);
plot(t_exe,raw_exe_ECG,'b');
xlim([0 t_exe(end)]);
xlabel('time [sec]');
ylabel('ECG signal [mv]');
title('Unfiltered ECG- after exercise');
subplot(2,1,2);
plot(t_exe,filt_exe_ECG,'b');
hold on;
scatter(t_exe(t_peak_exe_ECG),filt_exe_ECG(t_peak_exe_ECG),10,'filled','r');
xlim([0 t_exe(end)]);
ylim([min(filt_exe_ECG)-0.05 max(filt_exe_ECG)+0.1])
xlabel('time [sec]');
ylabel('ECG signal [mv]');
legend('Filtered ECG','R peaks')
title('Peak detection- after exercise');
hold off;


%% ECG Heart Rate Calculation

%rest
ind_peaks_diff=diff(t_peak_rest_ECG);
time_diff=ind_peaks_diff./Fs;
HR_rest_ECG = 60./time_diff;
HR_rest_ECG_mean = mean(HR_rest_ECG);
HR_rest_ECG_std = std(HR_rest_ECG);

%after exersice
ind_peaks_diff=diff(t_peak_exe_ECG);
time_diff=ind_peaks_diff./Fs;
HR_exe_ECG = 60./time_diff;
HR_exe_ECG_mean = mean(HR_exe_ECG);
HR_exe_ECG_std = std(HR_exe_ECG);

%%  2. PCG signal
%% filtering the signal

% taking the data from the current channel
raw_rest_PCG = rest(:, channel_PCG);
raw_exe_PCG = exercise(:, channel_PCG);

% normalization
raw_rest_PCG = raw_rest_PCG / max(abs(raw_rest_PCG)); 
raw_exe_PCG = raw_exe_PCG / max(abs(raw_exe_PCG));

% filtering with butterworth

[b,a]=butter(5,2*[60 150]/Fs,'bandpass');  %creates butterworth BPF of n=5 order
%apply BPF on signal
PCG_rest_filtered=filter(b,a,raw_rest_PCG); 
PCG_exe_filtered=filter(b,a,raw_exe_PCG);   


%% calculation of Shannon energy envelope 
%rest
Shanon_rest_PCG = zeros(1, length(PCG_rest_filtered));
for i = 1:length(PCG_rest_filtered)
    if PCG_rest_filtered(i) == 0
        Shanon_rest_PCG(i) = 0;
    else
        Shanon_rest_PCG(i) = ((abs(PCG_rest_filtered(i))^2) * log(abs(PCG_rest_filtered(i))^2));
    end
end

%smoothing 
Shanon_rest_PCG_smooth = zeros(1, length(PCG_rest_filtered));
for i = 1:length(Shanon_rest_PCG)
    if i <= 100
        Shanon_rest_PCG_smooth(i) = (-(1/i))*sum(Shanon_rest_PCG(1:i));
    else
        Shanon_rest_PCG_smooth(i) = (-(1/100))*sum(Shanon_rest_PCG(i-100:i));
    end
end

Shanon_clean_rest_PCG = zeros(1, length(raw_rest_PCG));
for i = 1:length(Shanon_rest_PCG_smooth)
    Shanon_clean_rest_PCG(i) = (Shanon_rest_PCG_smooth(i) - mean(Shanon_rest_PCG_smooth)) / (std(Shanon_rest_PCG_smooth));
end

%after exercise
Shanon_exe_PCG = zeros(1, length(PCG_exe_filtered));
for i = 1:length(PCG_exe_filtered)
    if PCG_exe_filtered(i) == 0
        Shanon_exe_PCG(i) = 0;
    else
        Shanon_exe_PCG(i) = ((abs(PCG_exe_filtered(i))^2) * log(abs(PCG_exe_filtered(i))^2));
    end
end


Shanon_exe_PCG_smooth = zeros(1, length(PCG_exe_filtered));
for i = 1:length(Shanon_exe_PCG)
    if i <= 100
        Shanon_exe_PCG_smooth(i) = (-(1/i))*sum(Shanon_exe_PCG(1:i));
    else
        Shanon_exe_PCG_smooth(i) = (-(1/100))*sum(Shanon_exe_PCG(i-100:i));
    end
end

Shanon_clean_exe_PCG = zeros(1, length(raw_exe_PCG));
for i = 1:length(Shanon_exe_PCG_smooth)
    Shanon_clean_exe_PCG(i) = (Shanon_exe_PCG_smooth(i) - mean(Shanon_exe_PCG_smooth)) / (std(Shanon_exe_PCG_smooth));
end


%% PCG signal Heart Rate Calculations

% rest
[peaks_values_PCG_rest,PCG_peaks_rest] = findpeaks(Shanon_clean_rest_PCG,'MinPeakDistance',0.75*Fs,'MinPeakHeight',0.7);

% HR calculation
ind_peaks_diff=diff(PCG_peaks_rest);
time_diff=ind_peaks_diff./Fs;
HR_rest_PCG = 60./time_diff;
HR_rest_PCG_mean = mean(HR_rest_PCG);
HR_rest_PCG_std = std(HR_rest_PCG);

% after exercise
[~,PCG_peaks_exe] = findpeaks(Shanon_clean_exe_PCG,'MinPeakDistance',0.5*Fs,'MinPeakHeight',0.004);%,'MinPeakDistance',ceil(1/max_HR));

% HR calculation
ind_peaks_diff=diff(PCG_peaks_exe);
time_diff=ind_peaks_diff./Fs;
HR_exe_PCG = 60./time_diff;
HR_exe_PCG_mean = mean(HR_exe_PCG);
HR_exe_PCG_std = std(HR_exe_PCG);

% plot peak detection
%unfiltered and filtered rest signals
figure;
subplot(2,1,1);
plot(t_rest,raw_rest_PCG,'b');
xlim([0 t_rest(end)]);
xlabel('time [sec]');
ylabel('ECG signal [mv]');
title('Unfiltered PCG- rest');
subplot(2,1,2);
plot(t_rest,PCG_rest_filtered,'b');
hold on
plot(t_rest,0.05*Shanon_clean_rest_PCG,'r','LineWidth',0.7);
hold on
scatter(t_rest(PCG_peaks_rest),0.05.*Shanon_clean_rest_PCG(PCG_peaks_rest),15,'filled','MarkerFaceColor','r');
legend('Filtered PCG signal','Shannon envelope');
title('PCG peak detection- rest');
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_rest(end)]);
ylim([min(PCG_rest_filtered)-0.05 max(0.05.*Shanon_clean_rest_PCG)+0.05])
hold off;

%plot unfiltered and filtered after exercise
figure;
subplot(2,1,1);
plot(t_exe,raw_exe_PCG,'b');
xlim([0 t_exe(end)]);
xlabel('time [sec]');
ylabel('ECG signal [mv]');
title('Unfiltered PCG- after exercise');
subplot(2,1,2);
plot(t_exe,PCG_exe_filtered,'b');
hold on
plot(t_exe,0.1.*Shanon_clean_exe_PCG,'r','LineWidth',0.7);
hold on
scatter(t_exe(PCG_peaks_exe),0.1.*Shanon_clean_exe_PCG(PCG_peaks_exe),15,'filled','MarkerFaceColor','r');
legend('Filtered PCG signal','Shannon envelope');
title('PCG peak detection- after exercise');
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_exe(end)]);
ylim([min(PCG_exe_filtered)-0.1 max(0.1.*Shanon_clean_exe_PCG)+0.05])
hold off;


%% 3. BCG signal 
%% Reading and filtering the signal

%taking the current chanel
BCG_rest = rest(:, channel_BCG); 
BCG_exe =exercise(:, channel_BCG); 

[b,a]=butter(5,2*[5 40]/Fs,'bandpass');         %creates butterworth BPF of n=5 order

%applying BPF on signal
BCG_rest_filtered=filter(b,a,BCG_rest); 
BCG_exe_filtered=filter(b,a,BCG_exe);   


%% Shannon energy envelope calculations
%rest
Shanon_rest_BCG = zeros(1, length(BCG_rest_filtered));
for i = 1:length(BCG_rest_filtered)
    if BCG_rest_filtered(i) == 0
        Shanon_rest_BCG(i) = 0;
    else
        Shanon_rest_BCG(i) = ((abs(BCG_rest_filtered(i))^2) * log(abs(BCG_rest_filtered(i))^2));
    end
end


Shanon_rest_BCG_smooth = zeros(1, length(BCG_rest_filtered));
for i = 1:length(Shanon_rest_BCG)
    if i <= 100
        Shanon_rest_BCG_smooth(i) = (-(1/i))*sum(Shanon_rest_BCG(1:i));
    else
        Shanon_rest_BCG_smooth(i) = (-(1/100))*sum(Shanon_rest_BCG(i-100:i));
    end
end

Shanon_clean_rest_BCG = zeros(1, length(Shanon_rest_BCG_smooth));
for i = 1:length(Shanon_rest_BCG_smooth)
    Shanon_clean_rest_BCG(i) = (Shanon_rest_BCG_smooth(i) - mean(Shanon_rest_BCG_smooth)) / (std(Shanon_rest_BCG_smooth));
end

%active
Shanon_exe_BCG = zeros(1, length(BCG_exe_filtered));
for i = 1:length(BCG_exe_filtered)
    if BCG_exe_filtered(i) == 0
        Shanon_exe_BCG(i) = 0;
    else
        Shanon_exe_BCG(i) = ((abs(BCG_exe_filtered(i))^2) * log(abs(BCG_exe_filtered(i))^2));
    end
end


Shanon_exe_BCG_smooth = zeros(1, length(BCG_exe_filtered));
for i = 1:length(Shanon_exe_BCG)
    if i <= 100
        Shanon_exe_BCG_smooth(i) = (-(1/i))*sum(Shanon_exe_BCG(1:i));
    else
        Shanon_exe_BCG_smooth(i) = (-(1/100))*sum(Shanon_exe_BCG(i-100:i));
    end
end

Shanon_clean_exe_BCG = zeros(1, length(Shanon_exe_BCG_smooth));
for i = 1:length(Shanon_exe_BCG_smooth)
    Shanon_clean_exe_BCG(i) = ( mean(Shanon_exe_BCG_smooth) - Shanon_exe_BCG_smooth(i)) / (std(Shanon_exe_BCG_smooth));
end


%% HR Calculations
% rest
beat_time=0.85; %[sec] for HR=70 beat/min
window=beat_time*Fs;
peak_rest_BCG=[];
for w=1:window:length(Shanon_clean_rest_BCG)-window
    cur_window=Shanon_clean_rest_BCG(w:w+window);
    [~,idx_max]=max(cur_window);
    peak_rest_BCG=[peak_rest_BCG,idx_max+w];
end
% Calculates HR
ind_peaks_diff=diff(peak_rest_BCG);
time_diff=ind_peaks_diff./Fs;
HR_rest_BCG = 60./time_diff;
HR_rest_BCG_mean = mean(HR_rest_BCG);
HR_rest_BCG_std = std(HR_rest_BCG);

[~,peak_exe_BCG] = findpeaks(Shanon_clean_exe_BCG,'MinPeakDistance',0.6*Fs);
% after exercise
beat_time2=0.6; %[sec] for HR=100 beat/min
window2=beat_time2*Fs;
peak_exe_BCG=[];
count=0;
for w=1:window2:length(Shanon_clean_exe_BCG)-window2
    cur_window=Shanon_clean_exe_BCG(w:w+window2);
    [~,idx_max]=max(cur_window);
    peak_exe_BCG=[peak_exe_BCG,idx_max+w];
end

% eliminate closoe peaks
diff_peaks=diff(peak_exe_BCG);
eliminate=find(diff_peaks<0.3*Fs);
peak_exe_BCG(eliminate+1)=[];

% Calculates HR
ind_peaks_diff=diff(peak_exe_BCG);
time_diff=ind_peaks_diff./Fs;
HR_exe_BCG = 60./time_diff;
HR_exe_BCG_mean = mean(HR_exe_BCG);
HR_exe_BCG_std = std(HR_exe_BCG);

%plot unfiltered and filtered after BCG peak detection
figure;
subplot(2,1,1);
plot(t_rest,BCG_rest,'b');
xlim([0 t_rest(end)]);
xlabel('time [sec]');
ylabel('ECG signal [mv]');
title('Unfiltered BCG- rest');
subplot(2,1,2);
plot(t_rest,BCG_rest_filtered,'b');
hold on
plot(t_rest,0.3.*Shanon_clean_rest_BCG,'r','LineWidth',0.7);
hold on
scatter(t_rest(peak_rest_BCG),0.3.*Shanon_clean_rest_BCG(peak_rest_BCG),15,'filled','MarkerFaceColor','r');
legend('Filtered BCG signal','Shannon envelope');
title('BCG peak detection- rest');
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_rest(end)]);
ylim([min(BCG_rest_filtered)-0.05 max(0.3.*Shanon_clean_rest_BCG)+0.05])
hold off;

figure;
subplot(2,1,1);
plot(t_exe,BCG_exe,'b');
xlim([0 t_exe(end)]);
xlabel('time [sec]');
ylabel('ECG signal [mv]');
title('Unfiltered BCG- after exercise');
subplot(2,1,2);
plot(t_exe,BCG_exe_filtered,'b');
hold on
plot(t_exe,0.5.*Shanon_clean_exe_BCG,'r','LineWidth',0.7);
hold on
scatter(t_exe(peak_exe_BCG),0.5.*Shanon_clean_exe_BCG(peak_exe_BCG),15,'filled','MarkerFaceColor','r');
legend('Filtered BCG signal','Shannon envelope');
title('BCG peak detection- after exercise');
xlabel('time [sec]');
ylabel('Voltage [mV]');
xlim([0 t_exe(end)]);
ylim([min(BCG_exe_filtered)-0.1 max(0.5.*Shanon_clean_exe_BCG)+0.05])
hold off;


%% 5. Statistics 

%% variance: 
V_rest_ECG=HR_rest_ECG_std^2;
V_exe_ECG=HR_exe_ECG_std^2;
V_rest_PCG=HR_rest_PCG_std^2;
V_exe_PCG=HR_exe_PCG_std^2;
V_rest_BCG=HR_rest_BCG_std^2;
V_exe_BCG=HR_exe_BCG_std^2;


%% Plot HR for all signals
figure();
plot(t_rest(t_peak_rest_ECG(2:end)),HR_rest_ECG,'-*m');
hold on;
plot(t_rest(PCG_peaks_rest(2:end)),HR_rest_PCG,'-*c');
hold on
plot(t_rest(peak_rest_BCG(2:end)),HR_rest_BCG,'-*g');
xlim([0 t_rest(end)]);
xlabel('time [sec]');
ylabel('HR [BPM]');
title('HR measurments- rest');
legend('ECG signal','PCG signal','BCG signal');
hold off;

figure;
plot(t_exe(t_peak_exe_ECG(2:end)),HR_exe_ECG,'-*m');
hold on;
plot(t_exe(PCG_peaks_exe(2:end)),HR_exe_PCG,'-*c');
hold on
plot(t_exe(peak_exe_BCG(2:end)),HR_exe_BCG,'-*g');
xlim([0 t_exe(end)]);
xlabel('time [sec]');
ylabel('HR [BPM]');
title('HR measurments- after exercise');
legend('ECG signal','PCG signal','BCG signal');
hold off;

%% ks-test
% rest
[h_rest(1),p_rest(1),D_rest(1)]=kstest2(HR_rest_ECG(1:end-4)',HR_rest_PCG(2:end-3)');
[h_rest(2),p_rest(2),D_rest(2)]=kstest2(HR_rest_ECG(1:end-4)',HR_rest_BCG(2:end-3)');

% active
[h_exe(1),p_exe(1),D_exe(1)]=kstest2(HR_exe_ECG',HR_exe_PCG(2:end)');
[h_exe(2),p_exe(2),D_exe(2)]=kstest2(HR_exe_ECG',HR_exe_BCG(3:end)');

%% relative error
% rest
e_rest(1)=abs((HR_rest_ECG_mean-HR_rest_PCG_mean)/HR_rest_ECG_mean);
e_rest(2)=abs((HR_rest_ECG_mean-HR_rest_BCG_mean)/HR_rest_ECG_mean);

% active
e_exe(1)=abs((HR_exe_ECG_mean-HR_exe_PCG_mean)/HR_exe_ECG_mean);
e_exe(2)=abs((HR_exe_ECG_mean-HR_exe_BCG_mean)/HR_exe_ECG_mean);