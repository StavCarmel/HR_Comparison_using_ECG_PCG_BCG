function signal_new = filter_the_signal(ECG_SIG,fs)

%Creating the time vector
time= 0:1/fs:((length(ECG_SIG))-1)*(1/fs);


%create a BPF with cutoff frequincies related to the ECG signal
my_filter=designfilt('bandpassfir','FilterOrder',1000,'CutoffFrequency1',4,'CutoffFrequency2',40,'SampleRate',fs);
% fvtool(my_filter); %plot the filter

%create a BSF in order to make sure the noise between 45-55 Hz will be
%cleanded
bsFilt = designfilt('bandstopfir','FilterOrder',400,'CutoffFrequency1',45,'CutoffFrequency2',55,'SampleRate',fs);
% fvtool(bsFilt)

%now we use this 2 filters we designed
signal_new = filtfilt(my_filter,ECG_SIG);
signal_new=filtfilt(bsFilt,signal_new);


end

