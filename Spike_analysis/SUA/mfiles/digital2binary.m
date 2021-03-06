% By Joachim
function [recording_digital_binary] = digital2binary(recording_digital)
%% Converts the digital (TTL, Laser on/off) stimulus into a binary (1 (Laser On), 0 (Laser off)
% threshold is a parameter.
% threshold = 10 works normall
parameters = get_parameters;

stim_max                         = max(recording_digital);
recording_digital_binary = recording_digital > stim_max*parameters.preprocessing.nlx_load_Opto.digital2binary.Threshold;
end