function parameters = get_parameters

%% Basic information
param.stimFreq = [2 4 8 16 32 64];
param.powerEstFreq = 16;
param.constructs = {'CAG-ChR2ET/TC-2AtDimer2','CAG-2AtDimer2'};
param.CSCPFC= 20:32; % CSC to pick to calculate firingrate & probability for PFC
param.CSCHP= 4:16; % CSC to pick to calculate firingrate & probability for HP
param.Window_square = [-3 0 3 6]; % preStimStart, stimStart, stimEnd, postStimEnd
param.Window_sine = [-3 0 3 6];
param.Window_ramp = [-3 0 3 6];
% param.squareDuration = 3; % to be replade with param.Window_square
% param.sineDuration = 3; % to be replade with param.Window_sine
% param.rampDuration = 3; % to be replade with param.Window_ramp
param.stimulationTypes = {'square','ramp','sine'};
% param.Areas = {'PFCL2', 'intmCA1';'PFCL5', 'intmCA1'};
%% nlx_load
param.preprocessing.nlx_load_Opto.voltAdjust = 32.82; %adjust to microvolt
param.preprocessing.nlx_load_Opto.loadPeriod = [3 3]; % [preStim, postStim], in seconds
param.preprocessing.nlx_load_Opto.digital2binary.Threshold = 0.25;
param.preprocessing.nlx_load.digital2binary.Threshold= 0.25;
%% mainBaselineComparison
param.osc_detection.min_interval=0.2;   %min interval between two events (s)
param.osc_detection.min_duration=1;     %min duration of event (s)
param.osc_detection.thr_sd=2;           %number of standard deviations
param.osc_detection.sharpWaveThreshold = 10; %nr of std from baseline
%% Spike analysis
param.plotStimulationRasterPlot.binsize = 0.010; %sec
param.spikeanalysis.spikeDetection.threshold = 5; %nr of std from baseline
param.spikeanalysis.plotBaselineMUAPhase.freqBands = {'theta','beta','gamma'};
param.spikeanalysis.phaseFilters = {'LFP','theta', 'beta', 'gamma','HFO', 'MUA'};
param.spikeanalysis.ISI.binsize = 0.005; %sec
param.spikeanalysis.firingRate= 3.0; %sec
% param.spikeanalysis.firingRate= 0.1; %sec

% param.getStimulationWaveform.window = [0.005 0.020]; %sec, [preStim postStim]
param.getStimulationWaveform.window = [0.1 0.3]; %sec, [preStim postStim]

%getStimSpikeProbability
% param.getStimSpikeProbability.CSCprob = 20:32;
param.getStimSpikeProbability.spikeWindow = [0 0.015]; % time until spike has to occur

%mainMUAspiekResultsPowerEst
% param.mainMUAspikeResultsPowerEst.LaserOutputs = 0:5:100;

%getMeanFiringRate
% param.getMeanFiringRate.binsize = 1;
% param.getMeanFiringRate.MinNspikes = 100;
% param.getMeanFiringRate.spikeratio = 2;

%getStimulationFiringChannelsSquare
param.getStimulationFiringChannelsSquare.spikeProbThr=0.50;

%getStimulationFiringChannels
param.getStimulationFiringChannels.binsize = 1;
param.getStimulationFiringChannels.MinNspikes = 200; %150 for HP
param.getStimulationFiringChannels.spikeratio = 2;
%% nlx_load
param.preprocessing.nlx_load_Opto.voltAdjust = 32.82; %adjust to microvolt
param.preprocessing.nlx_load_Opto.loadPeriod = [3 3]; % [preStim, postStim], in seconds
param.preprocessing.nlx_load_Opto.digital2binary.Threshold = 0.25;


%% FrequencyBands
param.FrequencyBands.broadband       = [0.5 1500];
param.FrequencyBands.LFP             = [4 100];
param.FrequencyBands.theta           = [4 12];
param.FrequencyBands.beta            = [12 30];
param.FrequencyBands.gamma           = [30 100];
param.FrequencyBands.HFO             = [100 300];
param.FrequencyBands.MUA             = [500 5000];

%% oscViewer
param.mainOscViewer.preOsc=3; %sec
param.mainOscViewer.postOscStart=6; %sec

%% stimWavelet
param.stimWavelet.ymin = 4;
param.stimWavelet.ymax = 100;
param.wavelet_oscViewer.ymin = 4; %Hz
param.wavelet_oscViewer.ymax = 100; %Hz

%% cutandglue
param.powerSpectrum.windowSize = 1.5; %sec 
param.powerSpectrum.minwindowN = 40;
param.powerSpectrum.overlap = 0; %sec
param.powerSpectrum.nfft = 6400; % 0.5 Hz resolution
param.powerSpectrum.maxFreq = 400; % in Hz
%% coherence
param.coherence.minwindowN=40;
param.coherence.windowSize=1.5; % three oscillations inside (down to 1.5 Hz)
param.coherence.nfft=6400; % 0.5 Hz resolution
param.coherence.overlap=0; % sec
param.coherence.maxFreq=100; %in Hz
param.coherence.baseline.groups = {  'Cg', 'PL';...
                            'Cg', 'IL';...
                            'PL', 'IL';...
                            'Cg', 'HPreversal';...
                            'PL', 'HPreversal';...
                            'IL', 'HPreversal'};
% tapers chronux coherence
param.chronux=struct('tapers',[3 5],'pad',0,'Fs',3200,'fpass'...
                     ,[0 100],'err',[1 0.05],'segave',0);

                
                
parameters = param;

end