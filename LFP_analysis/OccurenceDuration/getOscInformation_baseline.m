function [OscStructure] = getOscInformation_baseline(time,signal, samplingrate,oscStart,oscEnd)
% filter trace
parameters=get_parameters;
signalLFP = ZeroPhaseFilter(signal, samplingrate,parameters.FrequencyBands.LFP);

oscDiff=oscEnd-oscStart;
[Ampl,~]=HilbertTransf(signalLFP);

durOsc=[];
amplOsc=[];
amplOsc_max=[];
all_prop_Osc=[];

for osc = 1:length(oscStart)
    durOsc(osc)=oscDiff(osc)/samplingrate;
    amplOsc(osc)=mean(Ampl(oscStart(osc):oscEnd(osc)));
    amplOsc_max(osc)=max(Ampl(oscStart(osc):oscEnd(osc)));
    all_prop_Osc(osc,:)=[durOsc(osc),amplOsc(osc),amplOsc_max(osc)];
end

% correct oscStartTime
    oscStart = oscStart./samplingrate*10^6+time(1);
    oscEnd = oscEnd./samplingrate*10^6+time(1);

meanDurOsc=[];
stdDurOsc=[];
semDurOsc=[];
meanAmplOsc=[];
stdAmplOsc=[];
semAmplOsc=[];
meanAmplOsc_max=[];
stdAmplOsc_max=[];
semAmplOsc_max=[];
durRecording=[];
occurrenceOsc=[];

if ~isempty(durOsc)
% mean duration
        meanDurOsc=mean(durOsc);
        stdDurOsc=std(durOsc);
        semDurOsc=std(durOsc)/sqrt(length(durOsc));
% mean amplitude

        meanAmplOsc=mean(amplOsc);
        stdAmplOsc=std(amplOsc);
        semAmplOsc=std(amplOsc)/sqrt(length(amplOsc));

% mean max amplitude
        meanAmplOsc_max=mean(amplOsc_max);
        stdAmplOsc_max=std(amplOsc_max);
        semAmplOsc_max=std(amplOsc_max)/sqrt(length(amplOsc_max));

% occurence
        durRecording=(length(time)/samplingrate)/60;    %minutes
        occurrenceOsc=length(oscStart)/durRecording;
else
    'No osc detected in this signal'
end
%% Spectrum
%         testvector_osc =sum(pxx_osc,2); %with testvector find those oscillations for whom mtspectrumc cannot be calculated (e.g. first or last oscillation)
%         a=find(testvector_osc)~=0;
%         mean_pxx_osc=mean(pxx_osc(a,:),1); %calculate power-mean for all osc that are not 0


OscStructure.example = ['mean-', 'std-', 'sem'];
OscStructure.meanDurOsc = [meanDurOsc,stdDurOsc,semDurOsc];
OscStructure.meanAmplOsc = [meanAmplOsc,stdAmplOsc,semAmplOsc];
OscStructure.meanAmplOsc_max = [meanAmplOsc_max,stdAmplOsc_max,semAmplOsc_max];
OscStructure.all_prop_Osc = all_prop_Osc;
OscStructure.occurrenceOsc = occurrenceOsc;
% OscStructure.mean_pxx_osc = mean_pxx_osc;
% OscStructure.f_osc = f_osc;
OscStructure.OscTimes = [oscStart;oscEnd];
end