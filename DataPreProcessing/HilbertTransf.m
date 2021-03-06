%% Hilbert-Transformation

function [signal_ampl,signal_phase]=HilbertTransf(signal)
hilb_data=hilbert(signal);     % Hilbert-Transformation
clearvars signal
signal_ampl=abs(hilb_data);                %Amplitude
signal_phase=angle(hilb_data);        %Phase
end