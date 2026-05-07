clear dca
clc
close

% Step: Turn off firewall, allow inbound access


% Create connection to TI Radar board and DCA1000EVM Capture card
dca = dca1000("AWR6843ISK");

% Set recording duration to 100 seconds
dca.RecordDuration = 100;
% Define a variable to set the sampling rate in Hz for the
% phased.RangeResponse object. Because the dca1000 object provides the
% sampling rate in kHz, convert this rate to Hz.
fs = dca.ADCSampleRate*1e3;
% Define a variable to set the FMCW sweep slope in Hz/s for the
% phased.RangeResponse object. Because the dca1000 object provides the
% sweep slope in GHz/us, convert this sweep slope to Hz/s.
sweepSlope = dca.SweepSlope * 1e12;
% Define a variable to set the number of range samples
nr = dca.SamplesPerChirp;
% Create phased.RangeResponse System object that performs range filtering
% on fast-time (range) data, using an FFT-based algorithm
rangeresp = phased.RangeResponse(RangeMethod = 'FFT',...
                                 RangeFFTLengthSource = 'Property',...
                                 RangeFFTLength = nr, ...
                                 SampleRate = fs, ...
                                 SweepSlope = sweepSlope, ...
                                 ReferenceRangeCentered = false);
% The first call of the dca1000 object may take longer due to the
% configuration of the radar and the DCA1000EVM. To exclude the configuration
% time from the loop's duration, make the first call to the dca1000 object
% before entering the loop.
iqData = dca();

% Popup window to specify save filename
defaultName = "radar_capture.mat";
[file, path] = uiputfile("*.mat", "Save IQ data as", defaultName);

if isequal(file, 0)
    disp("Save canceled.");
else
    save(fullfile(path, file), "iqData", "fs", "sweepSlope", "nr");
    fprintf("Saved file: %s\n", fullfile(path, file));
end

clear dca