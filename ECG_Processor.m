%NAME: OWUSU ISAAC
%PROGRAM NAME: ASSIGNMENT_5
%UNIVERSITY NUMBER : N13863709
% DATE: 17/12/21

%Opening header_file
header = fopen('header_file.txt','r');
h = textscan(header,'%s');
fclose(header);
header = string(h{:});
%Initializing variables
sampling_freq = str2double(header(3));
num_Data = str2double(header(4));
num_Bits = str2double(header(6));
conversion_rate = str2double(header(7));
signal = load('a01.txt');
signal= (signal(:,2).*200);
f = sampling_freq;
%Taking 4 hours of total time to study
signal = double(signal(1:14400*sampling_freq));
time = 1/f:1/f:(length(signal)/f);
time=time.';
grid on
subplot(2,3,1);
plot(time(1:200),signal(1:200));
xlabel('time');
ylabel('Unsample signal');

new_f = 500;
%Resampling signals
resample_signal  = resample(signal, new_f,sampling_freq);
new_time = 1/new_f:1/new_f:(length(signal)/f);
new_time = new_time.';
subplot(2,3,2);
plot(new_time(1:1000),resample_signal (1:1000));
xlabel('time');
ylabel('1st sample signal');
%Passing low filter
d = designfilt('lowpassfir', 'Filterorder', 20, 'CutoffFrequency', 11, 'SampleRate', new_f);
filt_ECG1 = filter(d, resample_signal );
subplot(2,3,3);
%Plotting graphs
plot(new_time(1:1000),filt_ECG1(1:1000));
xlabel('time');
ylabel('signal');

hold on;
subplot(2,3,4);
d = designfilt('lowpassfir', 'Filterorder',5, 'CutoffFrequency', 11, 'SampleRate', new_f);
filt_ECG = filter(d, resample_signal );
%Slitting into 60 second domains
filt_ECG60= reshape(filt_ECG,[new_f*60, (length(filt_ECG)/30000)]).';
peak = [];
for i = 2:length(filt_ECG60(1,:))-1
    if(filt_ECG60(1, i)>100 && (filt_ECG60(1, i)>filt_ECG60(1, i-1))&&(filt_ECG60(1, i)> filt_ECG60(1, i+1)))
        peak = [peak, i];
    end
end
heart_rate = [];
length_peak = length(peak);
peak = peak/new_f;
for i=2:length_peak
    diff = peak(i)-peak(i-1);
    heart_rate = [heart_rate, 60/diff];
end
average_heart_rate = mean(heart_rate);
max_heart_rate = max (heart_rate);
min_heart_rate = min (heart_rate);
heart_rateData(1, 1) = average_heart_rate;
heart_rateData(1, 2) = max_heart_rate;
heart_rateData(1, 3) = min_heart_rate;
for j = 2:length(filt_ECG60(:,1))
peak = [];
%Computing peaks
for i = 2:length(filt_ECG60(j,:))-1
    if(filt_ECG60(j, i)>100 && (filt_ECG60(j, i)>filt_ECG60(j, i-1))&&(filt_ECG60(j, i)> filt_ECG60(j, i+1)))
        peak = [peak, i];
    end
end
heart_rate = [];
length_peak = length(peak);
peak = peak./new_f;
for i=2:length_peak
    diff = peak(i)-peak(i-1);
    heart_rate = [heart_rate, 60/diff];
end
peak = peak(2:end);
average_heart_rate = mean(heart_rate);
max_heart_rate = max (heart_rate);
min_heart_rate = min (heart_rate);
heart_rateData(j, 1) = average_heart_rate;
heart_rateData(j, 2) = max_heart_rate;
heart_rateData(j, 3) = min_heart_rate;
end
time_min = 1:1:240;

subplot(2,3,4);
plot(time_min,heart_rateData(:,1));
xlabel('time(minutes)');
ylabel('Heart rate');


PSDXArray = [];
filt_ECG_4(1,:) = resample(filt_ECG60(1,:),4,500);
Fs=4;
f=length(filt_ECG_4(1,:));
xdft = fft(filt_ECG_4(1,:)); %Do this to the ECG data, not the heart_rate
PSDX = (1/(Fs*f))* abs(xdft).^2;
PSDXArray = [PSDXArray, PSDX(1:f/2)];
freq = 1/f:Fs/f:Fs/2;
subplot(2,3,5);
plot(freq, PSDXArray(1,:));
xlabel('frequency');
ylabel('Power');

for i = 2:length(filt_ECG60(:,1))
filt_ECG_4(i,:) = resample(filt_ECG60(i,:),4,500);
Fs=4;
f=length(filt_ECG_4(i,:));
xdft = fft(filt_ECG_4(i,:)); %Do this to the ECG data, not the heart_rate
PSDX = (1/(Fs*f))* abs(xdft).^2;
PSDXArray = [PSDXArray; PSDX(1:f/2)];
end
subplot(2,3,6);
pcolor(PSDXArray.');
shading flat;
hold all;
%Graphing values
plot(heart_rateData(:,1),'r');
xlabel('time(minutes)');
ylabel('frequency');
n = fopen('a01.apn.txt','r');
m = textscan(n,'%c');
fclose(header);
z = char(m{:});
apnea1 = [];
for i=1:length(z);
if z(i)=='N'||z(i)=='A'
    apnea1 = [apnea1, z(i)];
end
end
apnea = [];
%Loop for comparison
for i=1:length(heart_rateData(:,1))
    if heart_rateData(i,1)>71
        apnea = [apnea, 'N'];
    else
        apnea = [apnea, 'A'];
    end
end
ctr=0;
for i=1:length(apnea)
    if apnea(i)==apnea1(i)
        ctr=ctr+1;
    end
end
%Computing accuracy
accuracy = (ctr/length(apnea))*100;
