function r_detect=find_R_func(to_filter,ECG,fs)
if to_filter==1
    ECG=filter_the_signal(ECG,30,38,fs);
elseif to_filter==0
    ECG=smooth(ECG);%smooth the signal so the noise will not disturb the detection when checking the derivatives
    ECG=ECG';
    %check if the ECG signal is upside down
    if max(abs(ECG))> max(ECG)
    ECG = -(ECG);
    end
end

y0=abs([ECG(2:end,1)',0]-[0,ECG(1:end-1,1)']);%first derivative
y1=abs([ECG(3:end,1)',0,0]-2*ECG'+[0,0,ECG(1:end-2,1)']);%second derivative
y2=1.3*y0+1.1*y1;
counter=0;
r_detect=[];
time_diff_R=round(fs/3)-1;%requirement of the time difference between R peaks 
max_val_vec=[];
min_amp=max(ECG)*0.2;%minimal amplitude required for detection
for i=1:length(y2)
    if y2(i)>=90/fs
        counter=counter+1;
        if counter==round(fs/200)%minimal sequence reqiured 
            look_for_max=i-round(fs/30):i+round(fs/30);%look for the R peak as the max at this enviroment
            if i-time_diff_R>round(fs/3)%if we are at the end of the signal
                if i>length(y2)-round(fs/30)
                    look_for_max=i-round(fs/30):i+(length(y2)-i);
                end
                if i<round(fs/30)
                    look_for_max=1:i+round(fs/30);
                end
                r_vec=ECG(look_for_max);
                [max_val,i_r]=max(r_vec);
                if max_val>min_amp
                    r_detect=[r_detect,i-(round(fs/30)+1)+i_r];
                    time_diff_R=i;
                    counter=0;
                end
            end
        end
        continue
    end
    counter=0;
end