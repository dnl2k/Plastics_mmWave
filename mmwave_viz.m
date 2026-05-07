clear all

% Create connection to TI Radar board and DCA1000EVM Capture card
% dca = dca1000("AWR6843ISK");

% Define a variable to set the sampling rate in Hz for the
% phased.RangeResponse object. Because the dca1000 object provides the
% sampling rate in kHz, convert this rate to Hz.
% fs = dca.ADCSampleRate*1e3; % 6144000
fs = 6144000;

% Define a variable to set the FMCW sweep slope in Hz/s for the
% phased.RangeResponse object. Because the dca1000 object provides the
% sweep slope in GHz/us, convert this sweep slope to Hz/s.
% sweepSlope = dca.SweepSlope * 1e12; % 8.499999998956919e+13
sweepSlope = 8.499999998956919e+13;

% Define a variable to set the number of range samples
% nr = dca.SamplesPerChirp;   % 240
nr = 240;

clear dca

%%
close all
% Create phased.RangeResponse System object that performs range filtering
% on fast-time (range) data, using an FFT-based algorithm
rangeresp = phased.RangeResponse(RangeMethod = 'FFT',...
                                 RangeFFTLengthSource = 'Property',...
                                 RangeFFTLength = nr, ...
                                 SampleRate = fs, ...
                                 SweepSlope = sweepSlope, ...
                                 ReferenceRangeCentered = false);

% List your data files
datafiles = "D:\ClassCode\UB\CSE_711\Data\260501\"+[
    "bkg_4.mat"
    % "metal"
    % "background_5.mat"
    "HDPE_4.mat"
    "LDPE_4.mat"
    "PET_4.mat"
    "PP_4.mat"
    "PS_4.mat"
    "PVC_4.mat"
];

figure;

for i = 1:length(datafiles)

    % Load each file
    S = load(datafiles(i));   % assumes each file contains variable iqData
    iqData = S.iqData;

    % Get data from first receiver antenna
    iqData_rx1 = squeeze(iqData(:,1,:));
    
    iqData_rx1_avg = mean(iqData_rx1, 2);

    % Create subplot
    subplot(4, 2, i);

    % Plot range response
    plotResponse(rangeresp, iqData_rx1_avg);

    % Limit x-axis to 0–30 cm
    xlim([0 0.5]);

    % Optional formatting
    title(strrep(datafiles(i), "_", "\_"));
    xlabel("Range (m)");
end


%% Amplitude-Phase scatterplot at selected range bin
close all

% -----------------------------
% Settings
% -----------------------------
dataDir = "D:\ClassCode\UB\CSE_711\Data\260501\";

targetRange = 0.24;   % meters, e.g., 0.20 m = 20 cm
rx = 1;               % receiver antenna to use

materialTypes = ["HDPE", "LDPE", "PET", "PP", "PS", "PVC"];
% Add "bkg" if you also want background points:
% materialTypes = ["bkg", "HDPE", "LDPE", "PET", "PP", "PS", "PVC"];

% -----------------------------
% Radar range axis
% -----------------------------
c = physconst("LightSpeed");
Nrange = nr;

rangeAxis = (0:Nrange/2-1) * fs / Nrange * c / (2 * sweepSlope);

% Find closest range bin
[~, targetBin] = min(abs(rangeAxis - targetRange));

fprintf("Selected target range = %.3f m\n", targetRange);
fprintf("Closest range bin = %d\n", targetBin);
fprintf("Actual range bin center = %.3f m\n", rangeAxis(targetBin));

% -----------------------------
% Find all .mat files
% -----------------------------
files = dir(fullfile(dataDir, "*.mat"));

allAmp_dB = [];
allPhase_rad = [];
allLabel = strings(0);
allSampleNum = [];

% -----------------------------
% Load each file and extract amplitude/phase
% -----------------------------
for k = 1:length(files)

    filename = string(files(k).name);
    filepath = fullfile(files(k).folder, files(k).name);

    % Parse naming convention: type_sampleNumber.mat
    [~, nameNoExt, ~] = fileparts(filename);

    tokens = regexp(nameNoExt, "^(?<type>[A-Za-z]+)_(?<sample>\d+)$", "names");

    if isempty(tokens)
        fprintf("Skipping file with unmatched name: %s\n", filename);
        continue;
    end

    sampleType = string(tokens.type);
    sampleNum = str2double(tokens.sample);

    % Keep only selected material types
    if ~ismember(sampleType, materialTypes)
        continue;
    end

    % Load data
    S = load(filepath);   % assumes variable iqData exists
    iqData = S.iqData;

    % Get RX data: samples x chirps
    iqData_rx = squeeze(iqData(:, rx, :));

    % Average over chirps in complex domain
    iqData_rx_avg = mean(iqData_rx, 2);

    % Range FFT
    X = fft(iqData_rx_avg, Nrange);
    X = X(1:Nrange/2);

    % Complex value at selected range bin
    xBin = X(targetBin);

    % Amplitude and phase
    amp_dB = 20*log10(abs(xBin) + eps);
    phase_rad = angle(xBin);   % wrapped phase in [-pi, pi]

    % Store
    allAmp_dB(end+1, 1) = amp_dB;
    allPhase_rad(end+1, 1) = phase_rad;
    allLabel(end+1, 1) = sampleType;
    allSampleNum(end+1, 1) = sampleNum;

end

% -----------------------------
% Scatterplot: phase vs amplitude
% -----------------------------
figure;
gscatter(allPhase_rad, allAmp_dB, allLabel, [], "o", 8);

xlabel("Phase at selected range bin (rad)");
ylabel("Amplitude at selected range bin (dB)");
title(sprintf("Amplitude-Phase Scatter at Range %.1f cm", rangeAxis(targetBin)*100));
grid on;

xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({"-\pi", "-\pi/2", "0", "\pi/2", "\pi"});

legend("Location", "best");

%% Per-chirp amplitude-phase scatterplot at selected range bin
clear all
close all

% -----------------------------
% Radar parameters
% -----------------------------
fs = 6144000;                         % Hz
sweepSlope = 8.499999998956919e+13;   % Hz/s
nr = 240;                             % samples per chirp

% -----------------------------
% User settings
% -----------------------------
dataDir = "D:\ClassCode\UB\CSE_711\Data\260501";

targetRange = 0.22;   % meters, e.g., 0.20 m = 20 cm
rx = 1;               % receiver antenna index

materialTypes = ["HDPE", "LDPE", "PET", "PP", "PS", "PVC"];
% Use this instead if you also want background:
% materialTypes = ["bkg", "HDPE", "LDPE", "PET", "PP", "PS", "PVC"];

% -----------------------------
% Range axis
% -----------------------------
c = physconst("LightSpeed");
Nrange = nr;

rangeAxis = (0:Nrange/2 - 1) * fs / Nrange * c / (2 * sweepSlope);

% Find closest range bin
[~, targetBin] = min(abs(rangeAxis - targetRange));

fprintf("Requested range: %.3f m\n", targetRange);
fprintf("Closest range bin: %d\n", targetBin);
fprintf("Actual range bin center: %.3f m\n", rangeAxis(targetBin));

% -----------------------------
% Find all .mat files
% -----------------------------
files = dir(fullfile(dataDir, "*.mat"));

allAmp_dB = [];
allPhase_rad = [];
allLabels = strings(0,1);
allSampleNums = [];
allChirpNums = [];
allFileNames = strings(0,1);

% -----------------------------
% Load files and extract per-chirp amplitude/phase
% -----------------------------
for k = 1:length(files)

    filename = string(files(k).name);
    filepath = fullfile(files(k).folder, files(k).name);

    % Remove extension
    [~, nameNoExt, ~] = fileparts(filename);

    % Parse filename: type_sampleNumber.mat
    % Example: HDPE_0.mat, LDPE_4.mat, bkg_4.mat
    tokens = regexp(nameNoExt, "^(?<type>[A-Za-z]+)_(?<sample>\d+)$", "names");

    if isempty(tokens)
        fprintf("Skipping unmatched filename: %s\n", filename);
        continue;
    end

    sampleType = string(tokens.type);
    sampleNum = str2double(tokens.sample);

    % Keep only selected types
    if ~ismember(sampleType, materialTypes)
        continue;
    end

    % Load file
    S = load(filepath);

    if ~isfield(S, "iqData")
        fprintf("Skipping %s because it does not contain iqData\n", filename);
        continue;
    end

    iqData = S.iqData;   % expected size: samples x RX antennas x chirps

    % Check RX index
    if size(iqData, 2) < rx
        fprintf("Skipping %s because RX%d does not exist\n", filename, rx);
        continue;
    end

    % Extract selected RX antenna
    iqData_rx = squeeze(iqData(:, rx, :));   % samples x chirps

    % Range FFT for every chirp
    X = fft(iqData_rx, Nrange, 1);

    % Keep positive range bins only
    X = X(1:Nrange/2, :);   % range bins x chirps

    % Extract complex value at selected range bin for every chirp
    xBin = X(targetBin, :);   % 1 x numChirps

    % Per-chirp amplitude and phase
    amp_dB = 20*log10(abs(xBin) + eps);
    phase_rad = angle(xBin);   % wrapped phase, range [-pi, pi]

    numChirps = numel(xBin);

    % Store results
    allAmp_dB = [allAmp_dB; amp_dB(:)];
    allPhase_rad = [allPhase_rad; phase_rad(:)];
    allLabels = [allLabels; repmat(sampleType, numChirps, 1)];
    allSampleNums = [allSampleNums; repmat(sampleNum, numChirps, 1)];
    allChirpNums = [allChirpNums; (1:numChirps).'];
    allFileNames = [allFileNames; repmat(filename, numChirps, 1)];

end

% -----------------------------
% Plot scatter: phase vs amplitude
% -----------------------------
figure;

gscatter(allPhase_rad, allAmp_dB, categorical(allLabels), [], "o", 8);

xlabel("Phase at selected range bin (rad)");
ylabel("Amplitude at selected range bin (dB)");
title(sprintf("Per-Chirp Amplitude-Phase Scatter at %.1f cm", rangeAxis(targetBin)*100));

grid on;
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({"-\pi", "-\pi/2", "0", "\pi/2", "\pi"});

legend("Location", "best");

% -----------------------------
% Optional: create result table
% -----------------------------
resultTable = table( ...
    allFileNames, ...
    allLabels, ...
    allSampleNums, ...
    allChirpNums, ...
    allAmp_dB, ...
    allPhase_rad, ...
    'VariableNames', ["File", "Type", "SampleNumber", "ChirpNumber", "Amplitude_dB", "Phase_rad"]);

disp(resultTable);