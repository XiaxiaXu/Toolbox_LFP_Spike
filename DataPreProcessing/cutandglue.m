
function [cutandglued,Nwindows,win_points]=...
    cutandglue(param,lfp,timestamps_samp)

%last modified 25.10.2016 - to use samples as input
%written by S.Gretenkord 10.02.15
%to cut and glue a signal with REM or SWS or discont events,
%in this version, severel segments from each event are included(as many as can be fitted,
%starting from the middle, and chopping off the end).
%note that analysis like pwelch or coherence shouls be performed with the
%same window length as the cutting ang glueing.

%Inputs:
%lfp: unfiltered LFP (one channel)
%timestamps: time stamps pointing to start end of events - in samples!!
%parameters such as sampling rate fs and window length in second and min
%window number needed to include data for subsequent analysis
%param.win_sec=0.75;
%param.fs=fs;

%Outputs:
%cutandglued: structure with  chopped and glued signals,
%Nwindows: number of windows that could be cut and glue for the given
%signal
% 2016.11.15 Xiaxia
%get number and lengths of 'multi events'
[multievent_num, col] = size(timestamps_samp);

%number of samples per window
win_points = floor(param.win_sec.*param.fs); %window in samples

if mod(win_points,2) > 0  %window in samples needs to be even
    win_points = win_points + 1;
end

%calculate segments of events - as many segments as can be fitted -
%centralised, i.e. remainder cut off right and left of multi event -
%with bias to left part if remainder uneven

xn_seg=[];
xn_seg_detrend=[];
good_events = 0;
Time=[];

for ev = 1:multievent_num
    
    p_start = timestamps_samp(ev,1);
    p_end = timestamps_samp(ev,2);
    length_ev=p_end-p_start+1;
    
    %calculate the number of windows that can be fit into this event and
    %remainder
    N_windows=floor(length_ev/win_points);
    remain=mod(length_ev,win_points);
    
    if mod(remain,2) >0 %uneven
        remain=remain-1; %part at being of multievent that is discarded one sample shorter than the discarded part a the end
    end
    
    for ii=1:N_windows
        seg_start= p_start + remain/2 +  (ii-1)*win_points;
        if seg_start > 0
            good_events = good_events + 1;
            Time(good_events,:)=seg_start : seg_start+ win_points-1;
            %get segment of requested window length
            lfpSegX= lfp( seg_start : seg_start+ win_points-1);
            
            %create cut and glued trace
            xn_seg(good_events,:) = lfpSegX ;
            
            %detrended cut and glued trace
            %xn_seg_detrend(good_events,:) = lfpSegX - mean(lfpSegX);
            %xn_seg_detrend(good_events,:) = detrend(lfpSegX);% 2016.11.15 Xiaxia
            % subtract power line noise and its harmonics , 2016.11.15
            % Xiaxia , Chronux box
            %[xn_seg_no_harmonics(:,good_events),~]=rmlinesmovingwinc(xn_seg_detrend(good_events,:)',[0.1 0.05],100,param);
        end
    end %windowed one event
end %windowed many events

Nwindows=size(xn_seg,1);
cutandglued.timepoint=Time;
cutandglued.xn=xn_seg;
%cutandglued.xndetrend= xn_seg_detrend;
%cutandglued.xn_no_harmonics=xn_seg_no_harmonics';

end



