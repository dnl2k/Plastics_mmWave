%% Per-chirp amplitude-phase scatterplot with average background subtraction
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

% -----------------------------
% Range axis
% -----------------------------
c = physconst("LightSpeed");
Nrange = nr;

rangeAxis = (0:Nrange/2 - 1) * fs / Nrange * c / (2 * sweepSlope);

[~, targetBin] = min(abs(rangeAxis - targetRange));

fprintf("Requested range: %.3f m\n", targetRange);
fprintf("Closest range bin: %d\n", targetBin);
fprintf("Actual range bin center: %.3f m\n", rangeAxis(targetBin));

% ---------------------------------------------------------
% 1. Compute average background from bkg_0.mat, bkg_1.mat, ...
% ---------------------------------------------------------
bkgFiles = dir(fullfile(dataDir, "bkg_*.mat"));

if isempty(bkgFiles)
    error("No background files found with pattern bkg_*.mat");
end

bkgSum = [];
numBkgChirpsTotal = 0;

for k = 1:length(bkgFiles)

    bkgPath = fullfile(bkgFiles(k).folder, bkgFiles(k).name);
    S = load(bkgPath);

    if ~isfield(S, "iqData")
        fprintf("Skipping %s because it does not contain iqData\n", bkgFiles(k).name);
        continue;
    end

    iqBkg = S.iqData;   % samples x RX antennas x chirps

    if isempty(bkgSum)
        bkgSum = zeros(size(iqBkg,1), size(iqBkg,2));
    end

    % Average over chirps in this background file, then accumulate
    % This keeps complex amplitude and phase information
    bkgSum = bkgSum + sum(iqBkg, 3);   % samples x RX
    numBkgChirpsTotal = numBkgChirpsTotal + size(iqBkg, 3);

end

avgBkg = bkgSum / numBkgChirpsTotal;   % samples x RX

fprintf("Computed average background from %d files and %d total chirps.\n", ...
    length(bkgFiles), numBkgChirpsTotal);

% ---------------------------------------------------------
% 2. Load plastic files, subtract background, extract per-chirp amp/phase
% ---------------------------------------------------------
files = dir(fullfile(dataDir, "*.mat"));

allAmp_dB = [];
allPhase_rad = [];
allLabels = strings(0,1);
allSampleNums = [];
allChirpNums = [];
allFileNames = strings(0,1);

for k = 1:length(files)

    filename = string(files(k).name);
    filepath = fullfile(files(k).folder, files(k).name);

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

    % Skip background files in the scatterplot
    if sampleType == "bkg"
        continue;
    end

    % Keep only selected material types
    if ~ismember(sampleType, materialTypes)
        continue;
    end

    S = load(filepath);

    if ~isfield(S, "iqData")
        fprintf("Skipping %s because it does not contain iqData\n", filename);
        continue;
    end

    iqData = S.iqData;   % samples x RX antennas x chirps

    % Check size compatibility
    if size(iqData,1) ~= size(avgBkg,1) || size(iqData,2) ~= size(avgBkg,2)
        fprintf("Skipping %s because its size does not match background size\n", filename);
        continue;
    end

    if size(iqData, 2) < rx
        fprintf("Skipping %s because RX%d does not exist\n", filename, rx);
        continue;
    end

    % -----------------------------------------------------
    % Background subtraction in complex IQ domain
    % -----------------------------------------------------
    avgBkgCube = repmat(avgBkg, 1, 1, size(iqData,3));
    iqData_sub = iqData - avgBkgCube;

    % Extract selected RX antenna
    iqData_rx = squeeze(iqData_sub(:, rx, :));   % samples x chirps

    % Range FFT for every chirp
    X = fft(iqData_rx, Nrange, 1);

    % Keep positive range bins only
    X = X(1:Nrange/2, :);   % range bins x chirps

    % Extract complex value at selected range bin for every chirp
    xBin = X(targetBin, :);   % 1 x numChirps

    % Per-chirp amplitude and phase
    amp_dB = 20*log10(abs(xBin) + eps);
    phase_rad = angle(xBin);   % wrapped phase in [-pi, pi]

    numChirps = numel(xBin);

    % Store results
    allAmp_dB = [allAmp_dB; amp_dB(:)];
    allPhase_rad = [allPhase_rad; phase_rad(:)];
    allLabels = [allLabels; repmat(sampleType, numChirps, 1)];
    allSampleNums = [allSampleNums; repmat(sampleNum, numChirps, 1)];
    allChirpNums = [allChirpNums; (1:numChirps).'];
    allFileNames = [allFileNames; repmat(filename, numChirps, 1)];

end

% ---------------------------------------------------------
% 3. Plot per-chirp amplitude-phase scatterplot
% ---------------------------------------------------------
figure;

gscatter(allPhase_rad, allAmp_dB, categorical(allLabels), [], "o", 8);

xlabel("Phase after background subtraction (rad)");
ylabel("Amplitude after background subtraction (dB)");
title(sprintf("Per-Chirp Amplitude-Phase Scatter at %.1f cm", rangeAxis(targetBin)*100));

grid on;
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({"-\pi", "-\pi/2", "0", "\pi/2", "\pi"});

legend("Location", "best");

% ---------------------------------------------------------
% 4. Optional: save result table
% ---------------------------------------------------------
resultTable = table( ...
    allFileNames, ...
    allLabels, ...
    allSampleNums, ...
    allChirpNums, ...
    allAmp_dB, ...
    allPhase_rad, ...
    'VariableNames', ["File", "Type", "SampleNumber", "ChirpNumber", "Amplitude_dB", "Phase_rad"]);

disp(resultTable);

%% ---------------------------------------------------------
% 5. Train/Validation/Test split: 60/20/20 by file
%    Neural network: 3 hidden layers, 20 neurons each
%    Input: amplitude + phase
% ---------------------------------------------------------

close all
rng(1);   % reproducibility

% -----------------------------
% Features and labels
% -----------------------------
X = [resultTable.Amplitude_dB, resultTable.Phase_rad];
Y = categorical(resultTable.Type);
fileNames = resultTable.File;

% Remove invalid rows
validRows = all(isfinite(X), 2);
X = X(validRows, :);
Y = Y(validRows);
fileNames = fileNames(validRows);

% -----------------------------
% Split by file, not by chirp
% This avoids putting chirps from the same file in train/valid/test.
% -----------------------------
classes = categories(Y);

trainFiles = strings(0,1);
validFiles = strings(0,1);
testFiles  = strings(0,1);

for c = 1:numel(classes)

    className = string(classes{c});

    classRows = Y == className;
    classFiles = unique(fileNames(classRows));

    % Shuffle files for this class
    classFiles = classFiles(randperm(numel(classFiles)));

    nFiles = numel(classFiles);

    if nFiles < 3
        warning("Class %s has fewer than 3 files. 60/20/20 split may be incomplete.", className);
    end

    % 60/20/20 split
    nTrain = round(0.60 * nFiles);
    nValid = round(0.20 * nFiles);
    nTest  = nFiles - nTrain - nValid;

    % Make sure each split gets at least one file when possible
    if nFiles >= 3
        nTrain = max(1, nTrain);
        nValid = max(1, nValid);
        nTest  = max(1, nTest);

        % Adjust if rounding made total too large
        while nTrain + nValid + nTest > nFiles
            nTrain = max(1, nTrain - 1);
        end
    end

    trainFiles = [trainFiles; classFiles(1:nTrain)];
    validFiles = [validFiles; classFiles(nTrain+1:nTrain+nValid)];
    testFiles  = [testFiles;  classFiles(nTrain+nValid+1:end)];
end

trainIdx = ismember(fileNames, trainFiles);
validIdx = ismember(fileNames, validFiles);
testIdx  = ismember(fileNames, testFiles);

XTrain = X(trainIdx, :);
YTrain = Y(trainIdx);

XValid = X(validIdx, :);
YValid = Y(validIdx);

XTest = X(testIdx, :);
YTest = Y(testIdx);

fprintf("Train chirps:      %d\n", size(XTrain,1));
fprintf("Validation chirps: %d\n", size(XValid,1));
fprintf("Test chirps:       %d\n", size(XTest,1));

if isempty(XValid)
    error("Validation set is empty. You need more files per class for a 60/20/20 file-level split.");
end

if isempty(XTest)
    error("Test set is empty. You need more files per class for a 60/20/20 file-level split.");
end

% -----------------------------
% Normalize using training set only
% -----------------------------
mu = mean(XTrain, 1);
sigma = std(XTrain, 0, 1);
sigma(sigma == 0) = 1;

XTrainNorm = (XTrain - mu) ./ sigma;
XValidNorm = (XValid - mu) ./ sigma;
XTestNorm  = (XTest  - mu) ./ sigma;

% -----------------------------
% Neural network architecture
% 3 hidden layers, 20 neurons each
% -----------------------------
numFeatures = size(XTrainNorm, 2);
numClasses = numel(categories(Y));

layers = [
    featureInputLayer(numFeatures, "Name", "input")

    fullyConnectedLayer(20, "Name", "fc1")
    reluLayer("Name", "relu1")

    fullyConnectedLayer(20, "Name", "fc2")
    reluLayer("Name", "relu2")

    fullyConnectedLayer(20, "Name", "fc3")
    reluLayer("Name", "relu3")

    fullyConnectedLayer(numClasses, "Name", "fc_out")
    softmaxLayer("Name", "softmax")
    classificationLayer("Name", "classification")
];

% -----------------------------
% Training options
% OutputNetwork = "best-validation-loss" returns the model
% with the best validation loss instead of the final epoch.
% -----------------------------
options = trainingOptions("adam", ...
    "MaxEpochs", 500, ...
    "MiniBatchSize", 32, ...
    "InitialLearnRate", 1e-3, ...
    "Shuffle", "every-epoch", ...
    "ValidationData", {XValidNorm, YValid}, ...
    "ValidationFrequency", 10, ...
    "ValidationPatience", 50, ...
    "OutputNetwork", "best-validation-loss", ...
    "Verbose", false, ...
    "Plots", "training-progress");

% -----------------------------
% Train network
% net is the best-validation-loss model
% -----------------------------
[net, trainInfo] = trainNetwork(XTrainNorm, YTrain, layers, options);

% ---------------------------------------------------------
% 6. Evaluate train, validation, and test sets
% ---------------------------------------------------------

YPredTrain = classify(net, XTrainNorm);
YPredValid = classify(net, XValidNorm);
YPredTest  = classify(net, XTestNorm);

trainAccuracy = mean(YPredTrain == YTrain);
validAccuracy = mean(YPredValid == YValid);
testAccuracy  = mean(YPredTest  == YTest);

fprintf("\nBest-validation model performance:\n");
fprintf("Train accuracy      = %.2f%%\n", trainAccuracy * 100);
fprintf("Validation accuracy = %.2f%%\n", validAccuracy * 100);
fprintf("Test accuracy       = %.2f%%\n", testAccuracy * 100);

% ---------------------------------------------------------
% 7. Confusion matrices
% ---------------------------------------------------------

figure;
cmTrain = confusionchart(YTrain, YPredTrain);
cmTrain.Title = sprintf("Train Confusion Matrix | Accuracy = %.2f%%", trainAccuracy * 100);
cmTrain.RowSummary = "row-normalized";
cmTrain.ColumnSummary = "column-normalized";

figure;
cmValid = confusionchart(YValid, YPredValid);
cmValid.Title = sprintf("Validation Confusion Matrix | Accuracy = %.2f%%", validAccuracy * 100);
cmValid.RowSummary = "row-normalized";
cmValid.ColumnSummary = "column-normalized";

figure;
cmTest = confusionchart(YTest, YPredTest);
cmTest.Title = sprintf("Test Confusion Matrix | Accuracy = %.2f%%", testAccuracy * 100);
cmTest.RowSummary = "row-normalized";
cmTest.ColumnSummary = "column-normalized";

% ---------------------------------------------------------
% 8. Scatterplots for train, validation, and test sets
% ---------------------------------------------------------

figure;

subplot(1,3,1);
gscatter(XTrain(:,2), XTrain(:,1), YTrain, [], "o", 8);
xlabel("Phase (rad)");
ylabel("Amplitude (dB)");
title("Train Set");
grid on;
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({"-\pi", "-\pi/2", "0", "\pi/2", "\pi"});

subplot(1,3,2);
gscatter(XValid(:,2), XValid(:,1), YValid, [], "o", 8);
xlabel("Phase (rad)");
ylabel("Amplitude (dB)");
title("Validation Set");
grid on;
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({"-\pi", "-\pi/2", "0", "\pi/2", "\pi"});

subplot(1,3,3);
gscatter(XTest(:,2), XTest(:,1), YTest, [], "o", 8);
xlabel("Phase (rad)");
ylabel("Amplitude (dB)");
title("Test Set");
grid on;
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({"-\pi", "-\pi/2", "0", "\pi/2", "\pi"});

% ---------------------------------------------------------
% 9. Decision boundary of best-validation model
% ---------------------------------------------------------

figure;

ampMin = min(X(:,1));
ampMax = max(X(:,1));

phaseMin = -pi;
phaseMax = pi;

[phaseGrid, ampGrid] = meshgrid( ...
    linspace(phaseMin, phaseMax, 300), ...
    linspace(ampMin, ampMax, 300));

XGrid = [ampGrid(:), phaseGrid(:)];
XGridNorm = (XGrid - mu) ./ sigma;

YGridPred = classify(net, XGridNorm);

% Decision regions
gscatter(phaseGrid(:), ampGrid(:), YGridPred, [], ".", 3);
hold on;

% Overlay actual data
gscatter(X(:,2), X(:,1), Y, [], "o", 8);

xlabel("Phase after background subtraction (rad)");
ylabel("Amplitude after background subtraction (dB)");
title(sprintf("Best-Validation NN Decision Boundary | Train %.1f%% | Valid %.1f%% | Test %.1f%%", ...
    trainAccuracy*100, validAccuracy*100, testAccuracy*100));

grid on;
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({"-\pi", "-\pi/2", "0", "\pi/2", "\pi"});

legend("Location", "bestoutside");

% ---------------------------------------------------------
% 10. Optional: plot validation accuracy/loss history
% ---------------------------------------------------------

figure;

subplot(2,1,1);
plot(trainInfo.TrainingLoss, "LineWidth", 1.2);
hold on;
if isfield(trainInfo, "ValidationLoss")
    plot(trainInfo.ValidationLoss, "LineWidth", 1.2);
    legend("Training Loss", "Validation Loss");
else
    legend("Training Loss");
end
xlabel("Iteration");
ylabel("Loss");
title("Training and Validation Loss");
grid on;

subplot(2,1,2);
if isfield(trainInfo, "ValidationAccuracy")
    plot(trainInfo.ValidationAccuracy, "LineWidth", 1.2);
    xlabel("Validation Check");
    ylabel("Validation Accuracy (%)");
    title("Validation Accuracy");
    grid on;
else
    text(0.1, 0.5, "ValidationAccuracy not available in trainInfo for this MATLAB version.");
    axis off;
end

% ---------------------------------------------------------
% 11. Save best-validation model
% ---------------------------------------------------------

modelFile = "amp_phase_nn_best_valid_model.mat";

save(modelFile, ...
    "net", ...
    "mu", ...
    "sigma", ...
    "trainAccuracy", ...
    "validAccuracy", ...
    "testAccuracy", ...
    "trainFiles", ...
    "validFiles", ...
    "testFiles", ...
    "trainInfo");

fprintf("Saved best-validation model to %s\n", modelFile);