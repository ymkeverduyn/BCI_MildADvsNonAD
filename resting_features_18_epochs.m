%% Load dataset 

% Load resting state eeg -- this needs to be changed for each subject manually
EEG = pop_loadset('filename', 'sub-012_task-Rest_eeg.set', 'filepath', 'C:\Users\jette\Documents\Uni\Tilburg\Semester2\BCI\Project\Dataset\ds003800-download\sub-012\eeg');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

%% Visual inspection 
% Plot the data to visually inspect it
%pop_eegplot(EEG, 1, 1, 1);

%% Feature extraction (split from below)
% 1) Discrete wavelet transform -- decompose into sub-bands
% 2) Power spectral density of frequency bands 

% Ensure the Wavelet Toolbox is installed
if ~license('test', 'Wavelet_Toolbox')
    error('Wavelet Toolbox is not installed.');
end

% Ensure the Signal Processing Toolbox is installed
if ~license('test', 'Signal_Toolbox')
    error('Signal Processing Toolbox is not installed.');
end

% Parameters
channels = 1:EEG.nbchan; % Channels to analyze
fs = EEG.srate;          % Sampling frequency
windowLength = 3;        % Window length in seconds
numWindows = 18;         % Number of windows
samplesPerWindow = windowLength * fs; % Number of samples per window
bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
freqRanges = [0.5 4; 4 8; 8 12; 12 30; 30 100]; % Frequency ranges

% Choose Daubechies wavelet and level of decomposition
wavelet = 'db4';
level = 6; % Level of decomposition

% Initialize cell arrays to hold sub-bands and results
subBands = cell(length(channels), length(bands), numWindows);
psdResults = cell(length(channels), length(bands), numWindows);
varResults = zeros(length(channels), length(bands), numWindows);
ampSumResults = zeros(length(channels), length(bands), numWindows);
arParams = cell(length(channels), length(bands), numWindows);
arOrder = zeros(length(channels), length(bands), numWindows);

% Channel pairs for interhemispheric coherence
channelPairs = [1, 9; 8, 16; 2, 10; 7, 15; 3, 11; 6, 14; 4, 12; 5, 13];
coherenceResults = zeros(size(channelPairs, 1), length(bands), numWindows);

%% Create Directory for Saving Results
outputDir = 'New_EEG_Results_sub012';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Initialize results table
resultsTable = table();

%% DWT & PSD feature extraction for all channels and windows

% Loop through each window
for w = 1:numWindows
    % Define the start and end indices for the current window
    startIdx = (w - 1) * samplesPerWindow + 1;
    endIdx = startIdx + samplesPerWindow - 1;
    
    % Loop through each channel
    for ch = 1:length(channels)
        % Extract data for the current channel and window
        eegData = EEG.data(channels(ch), startIdx:endIdx);
        
        % Perform discrete wavelet transform (DWT)
        [c, l] = wavedec(eegData, level, wavelet);
        
        % Loop through each frequency band
        for b = 1:length(bands)
            % Find the levels corresponding to the frequency range
            freqRange = freqRanges(b, :);
            approxLevel = find((fs ./ (2.^(1:level))) < freqRange(2), 1);
            detailLevel = find((fs ./ (2.^(1:level))) < freqRange(1), 1);
            
            % Extract approximation and detail coefficients
            if ~isempty(approxLevel)
                A = wrcoef('a', c, l, wavelet, approxLevel);
            else
                A = zeros(size(eegData));
            end
            if ~isempty(detailLevel)
                D = wrcoef('d', c, l, wavelet, detailLevel);
            else
                D = zeros(size(eegData));
            end
            
            % Sum the coefficients within the frequency band
            subBands{ch, b, w} = A + D;
            
            % Convert sub-band to double precision
            subBandDouble = double(subBands{ch, b, w});
            
            % Calculate PSD using Welch's method
            [pxx, f] = pwelch(subBandDouble, [], [], [], fs);
            psdResults{ch, b, w} = [f, pxx];
            
            % Calculate variance
            varResults(ch, b, w) = var(subBandDouble);
            
            % Calculate amplitude summation
            ampSumResults(ch, b, w) = sum(abs(subBandDouble));
            
            % Determine optimum AR model order using AIC
            numOrders = 10; % Maximum order to consider
            aicValues = zeros(1, numOrders);
            for order = 1:numOrders
                model = ar(subBandDouble, order, 'burg');
                aicValues(order) = aic(model);
            end
            [~, bestOrder] = min(aicValues);
            
            % Estimate AR model parameters using the optimal order
            arModel = ar(subBandDouble, bestOrder, 'burg');
            arParams{ch, b, w} = arModel.A;
            arOrder(ch, b, w) = bestOrder;
        end
    end
    
    %% Calculate Interhemispheric Coherence for Specified Channel Pairs

    % Loop through each channel pair
    for pairIdx = 1:size(channelPairs, 1)
        ch1 = channelPairs(pairIdx, 1);
        ch2 = channelPairs(pairIdx, 2);
        
        % Extract data for the channel pair
        eegData1 = EEG.data(ch1, startIdx:endIdx);
        eegData2 = EEG.data(ch2, startIdx:endIdx);
        
        % Loop through each frequency band
        for b = 1:length(bands)
            % Convert sub-band to double precision
            subBand1 = double(subBands{ch1, b, w});
            subBand2 = double(subBands{ch2, b, w});
            
            % Calculate CPSD using Welch's method
            [pxx1, f] = cpsd(subBand1, subBand2, [], [], [], fs);
            
            % Calculate PSD for each channel
            [psd1, ~] = pwelch(subBand1, [], [], [], fs);
            [psd2, ~] = pwelch(subBand2, [], [], [], fs);
            
            % Calculate coherence as described
            coherence = abs(pxx1) ./ sqrt(psd1 .* psd2);
            
            % Amplitude summation of coherence values
            coherenceResults(pairIdx, b, w) = sum(coherence);
        end
    end
end

%% Combine results into a single table

% Combine amplitude summation and variance results
for w = 1:numWindows
    windowLabel = sprintf('Window%d', w);
    for ch = 1:length(channels)
        for b = 1:length(bands)
            resultsTable.(sprintf('%s_Ch%d_%s_Var', windowLabel, ch, bands{b})) = varResults(ch, b, w);
            resultsTable.(sprintf('%s_Ch%d_%s_AmpSum', windowLabel, ch, bands{b})) = ampSumResults(ch, b, w);
        end
    end
    % Combine coherence results
    for pairIdx = 1:size(channelPairs, 1)
        for b = 1:length(bands)
            resultsTable.(sprintf('%s_Ch%d_Ch%d_%s_Coh', windowLabel, channelPairs(pairIdx, 1), channelPairs(pairIdx, 2), bands{b})) = coherenceResults(pairIdx, b, w);
        end
    end
end

% Save the results table to a CSV file
csvFilename = fullfile(outputDir, 'EEG_Features_Windows_sub012.csv');
writetable(resultsTable, csvFilename);

%%
plot(time, subBands{ch, b, w}); 