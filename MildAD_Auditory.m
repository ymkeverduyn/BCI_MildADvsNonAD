eeglab;

% Define subjects to process
subjects = {'002', '004', '005', '007', '011'};
base_path = '/Users/ymkeverduyn/Documents/MATLAB/ds003800-download/';


% Define channels to measure for each ERP component
channels = {'Fz', 'Cz', 'Pz'};

% Define the time windows for each ERP component
windows = struct('P300', [250, 500], 'N200', [180, 250], 'P200', [150, 220], 'N100', [80, 150]);

% Define the pre-stimulus and post-stimulus time in milliseconds
pre_stimulus_time = 0;  % 0 ms before the event
post_stimulus_time = 800;  % 800 ms after the event

for s = 1:length(subjects)
    subj = subjects{s};
    subj_path = fullfile(base_path, ['sub-' subj], 'eeg');
    
    % Load dataset for the current subject
    EEG = pop_loadset('filename', ['sub-' subj '_task-AuditoryGammaEntrainment_eeg.set'], 'filepath', subj_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    
    % Import the events data
    EEG = pop_importevent(EEG, 'event', fullfile(subj_path, ['sub-' subj '_task-AuditoryGammaEntrainment_events.tsv']), ...
                          'fields', {'seconds' 'duration' 'sample' 'trial_type' 'response_time' 'stim_file' 'value'}, ...
                          'skipline', 1, 'timeunit', 1, 'optimalign', 'off');
    
    % Clean up the event types
    for i = 1:length(EEG.event)
        if isempty(EEG.event(i).type) || ~ischar(EEG.event(i).type)
            % Assuming alternate pattern, assign types manually based on index
            if mod(i, 2) == 0  % Assuming even indices are 'rest' events
                EEG.event(i).type = '1';
            else  % Odd indices are 'auditory' events
                EEG.event(i).type = '2';
            end
        else
            % Ensure the event type is a character string
            EEG.event(i).type = char(EEG.event(i).type);
        end
    end
    
    % Debug: Print the imported events
    fprintf('Subject %s: Imported events:\n', subj);
    for i = 1:length(EEG.event)
        fprintf('Event %d: Type=%s, Latency=%.2f\n', i, EEG.event(i).type, EEG.event(i).latency);
    end

    % Check if event types are as expected
    event_types = {EEG.event.type}; % Ensure all event types are character arrays
    unique_event_types = unique(event_types(~cellfun('isempty', event_types))); % Filter out empty elements
    fprintf('Unique event types: %s\n', strjoin(unique_event_types, ', '));

    % Read the channel names from TSV file and assign to EEG
    chanlocs = readtable(fullfile(subj_path, ['sub-' subj '_task-AuditoryGammaEntrainment_channels.tsv']), 'FileType', 'text');
    for i = 1:height(chanlocs)
        EEG = pop_chanedit(EEG, 'changefield', {i 'labels' chanlocs.name{i}});
    end
    
    % Save the dataset with updated channel information
    EEG = pop_saveset(EEG, 'filename', ['sub-' subj '_eeg_with_channels.set'], 'filepath', subj_path);
    
    % Load the updated EEG dataset
    EEG = pop_loadset('filename', ['sub-' subj '_eeg_with_channels.set'], 'filepath', subj_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    
    % Plot the data to visually inspect it (optional)
    pop_eegplot(EEG, 1, 1, 1);
    
    % Define the event types for the start of auditory windows and rest windows
    auditory_event_type = '2';  % Event type indicating the start of the auditory window
    rest_event_type = '1';      % Event type indicating start of rest
    
    % Find the indices of the auditory start events
    auditory_start_indices = find(strcmp({EEG.event.type}, auditory_event_type));
    
    % Ensure no auditory event is too close to the start of the recording
    % Adjust min_latency based on your specific needs
    min_latency = -pre_stimulus_time * EEG.srate / 1000; % Minimum latency to allow for pre-stimulus interval
    if pre_stimulus_time > 0
        min_latency = 0; % Adjust for no pre-stimulus period
    end
    valid_auditory_start_indices = auditory_start_indices([EEG.event(auditory_start_indices).latency] > min_latency);
    
    % Debug: Print the number of auditory events found and their latencies
    fprintf('Found %d valid auditory start events\n', length(valid_auditory_start_indices));
    for idx = 1:length(valid_auditory_start_indices)
        fprintf('Event %d at latency %.2f seconds\n', idx, EEG.event(valid_auditory_start_indices(idx)).latency / EEG.srate);
    end
    
    % Find indices for each channel
    chan_indices = struct();
    for c = 1:length(channels)
        chan_indices.(channels{c}) = find(strcmp({EEG.chanlocs.labels}, channels{c}));
    end
    
    % Create epochs around the auditory events
    EEG = pop_epoch(EEG, {auditory_event_type}, [pre_stimulus_time/1000 post_stimulus_time/1000]);
    
    % Perform baseline correction (optional)
    % EEG = pop_rmbase(EEG, [pre_stimulus_time 0]);
    
    % Extract the ERP components for each epoch separately
    num_epochs = size(EEG.data, 3);
    peak_amplitudes = struct();
    peak_latencies = struct();
    
    % Initialize structures for peak amplitudes and latencies
    for c = 1:length(channels)
        chan = channels{c};
        peak_amplitudes.(chan) = struct('P300', zeros(num_epochs, 1), 'N200', zeros(num_epochs, 1), 'P200', zeros(num_epochs, 1), 'N100', zeros(num_epochs, 1));
        peak_latencies.(chan) = struct('P300', zeros(num_epochs, 1), 'N200', zeros(num_epochs, 1), 'P200', zeros(num_epochs, 1), 'N100', zeros(num_epochs, 1));
    end
    
    % Convert time windows to samples
    window_samples = structfun(@(x) round((x - pre_stimulus_time) / (post_stimulus_time - pre_stimulus_time) * size(EEG.data, 2)), windows, 'UniformOutput', false);
    
    % Create time vector
    time_vector = linspace(pre_stimulus_time, post_stimulus_time, size(EEG.data, 2));
    
    % Extract peak amplitudes and latencies for each ERP component
    for i = 1:num_epochs
        epoch_data = EEG.data(:, :, i);
        
        for c = 1:length(channels)
            chan = channels{c};
            chan_idx = chan_indices.(chan);
            
            % P300
            [peak_amplitudes.(chan).P300(i), peak_idx_p300] = max(epoch_data(chan_idx, window_samples.P300(1):window_samples.P300(2)));
            peak_latencies.(chan).P300(i) = time_vector(window_samples.P300(1) + peak_idx_p300 - 1);
            
            % N200
            [peak_amplitudes.(chan).N200(i), peak_idx_n200] = min(epoch_data(chan_idx, window_samples.N200(1):window_samples.N200(2)));
            peak_latencies.(chan).N200(i) = time_vector(window_samples.N200(1) + peak_idx_n200 - 1);
            
            % P200
            [peak_amplitudes.(chan).P200(i), peak_idx_p200] = max(epoch_data(chan_idx, window_samples.P200(1):window_samples.P200(2)));
            peak_latencies.(chan).P200(i) = time_vector(window_samples.P200(1) + peak_idx_p200 - 1);
            
            % N100
            [peak_amplitudes.(chan).N100(i), peak_idx_n100] = min(epoch_data(chan_idx, window_samples.N100(1):window_samples.N100(2)));
            peak_latencies.(chan).N100(i) = time_vector(window_samples.N100(1) + peak_idx_n100 - 1);
        end
    end
    
    % Plot the ERP waveforms for each epoch and each component
    components = fieldnames(windows);
    
    for c = 1:length(channels)
        chan = channels{c};
        
        for comp = 1:length(components)
            component = components{comp};
            figure;
            hold on;
            for i = 1:num_epochs
                plot(time_vector, EEG.data(chan_indices.(chan), :, i));
            end
            xlabel('Time (ms)');
            ylabel('Amplitude (µV)');
            title([component ' Components at ' chan ' for Each Epoch']);
            legend(arrayfun(@(x) ['Epoch ' num2str(x)], 1:num_epochs, 'UniformOutput', false));
            hold off;
            
            % Adjust the figure's PaperPosition property and save the figure
            fig = gcf;
            if isvalid(fig)  % Check if the figure handle is valid
                fig.PaperPositionMode = 'auto';
                fig_pos = fig.PaperPosition;
                fig.PaperSize = [fig_pos(3) fig_pos(4)];
                % Save the figure
                saveas(fig, fullfile(subj_path, [component '_' chan '_epochs.png']));
            else
                warning('Figure for component %s at channel %s could not be saved because it is not valid.', component, chan);
            end
        end
    end
    
    % Save the peak amplitudes and latencies to a CSV file
    for c = 1:length(channels)
        chan = channels{c};
        results_table = table((1:num_epochs)', peak_amplitudes.(chan).P300, peak_latencies.(chan).P300, ...
            peak_amplitudes.(chan).N200, peak_latencies.(chan).N200, peak_amplitudes.(chan).P200, peak_latencies.(chan).P200, ...
            peak_amplitudes.(chan).N100, peak_latencies.(chan).N100, ...
            'VariableNames', {'Epoch', ['PeakAmplitude_P300_' chan], ['PeakLatency_P300_' chan], ...
            ['PeakAmplitude_N200_' chan], ['PeakLatency_N200_' chan], ['PeakAmplitude_P200_' chan], ['PeakLatency_P200_' chan], ...
            ['PeakAmplitude_N100_' chan], ['PeakLatency_N100_' chan]});
        writetable(results_table, fullfile(subj_path, ['ERP_peak_features_' chan '.csv']));
    end
end
