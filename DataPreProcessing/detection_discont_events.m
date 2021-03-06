function [osc_start, osc_end,osc,rms,thr,minint,mindur,num_sd] = detection_discont_events(signal,fs,num_sd)
%detection using root-mean-square and gaussian noise fitting
%needs functions:
%calc_noice_thresh
%fit_noise
%ZeroPhaseFilter


%% osc_detection
parameters.osc_detection.min_interval=0.2;   %min interval between two events (s)
parameters.osc_detection.min_duration=1;     %min duration of event (s)
%parameters.osc_detection.thr_sd=3;         % 2,number of standard deviations


%% root-mean-square
lfp = ZeroPhaseFilter(signal,fs,[4,90]);
N = length(lfp);
rms = zeros(1,N);
win = 0.2*fs; %0.2
halfwin = round(win/2);
for i = halfwin : N-halfwin;
    rms(i) = norm(lfp(i-halfwin+1:i+halfwin))/sqrt(win);
end

%% calculate threshold and detect oscillations
%num_sd = parameters.osc_detection.thr_sd;   %number of standard deviations
thr = calc_noise_thresh(rms,num_sd);
osc = rms > thr;
osc([1,end]) = 0;

%% combine if short interval
osc_start = find(diff(osc)==1);
osc_end = find(diff(osc)==-1);
interval = (osc_start(2:end)-osc_end(1:end-1)) / fs;
for k=find(interval<parameters.osc_detection.min_interval)
    osc(osc_end(k):osc_start(k+1)) = 1;
end

%% delete if short event
osc_start = find(diff(osc)==1);
osc_end = find(diff(osc)==-1);
dur = (osc_end-osc_start) / fs;
for k = find(dur<parameters.osc_detection.min_duration)
    osc(osc_start(k):osc_end(k)) = 0;
end

%% start and end of events
osc_start = find(diff(osc)==1);
osc_end = find(diff(osc)==-1);

minint=parameters.osc_detection.min_interval;
mindur=parameters.osc_detection.min_duration;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%