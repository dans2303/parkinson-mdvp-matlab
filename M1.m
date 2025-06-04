clc;
close all;
clear all;

% Directory containing the audio files
audio_dir = 'E:\Kuliah NCU\PhD\Research\MDVP Matlab\SAMPLES';  % Replace with the path to your files
audio_files = dir(fullfile(audio_dir, '*.wav'));  % Get a list of all .wav files

% Initialize an empty table to store the results
results_table = table;

% Loop through each file
for i = 1:length(audio_files)
    % Get the full file path
    file_name = audio_files(i).name;
    full_file_path = fullfile(audio_dir, file_name);
    
    % Load the audio file
    [y, fs] = audioread(full_file_path);

    % Convert to mono if stereo
    if size(y, 2) == 2
        y = mean(y, 2);
    end

    % Calculate the duration of the recording or TSAM
    duration = length(y) / fs;

    % Pre-Emphasis Filter
    alpha = 0.97; % Pre-emphasis coefficient
    y_preemph = filter([1 -alpha], 1, y);

    % Normalize the pre-emphasized signal
    y_norm = y_preemph / max(abs(y_preemph));

    % High-Pass Filter Design (IIR for efficiency)
    high_cutoff_freq = 70 / (fs/2); % Normalized cutoff frequency (e.g., 70 Hz)
    hp_order = 4;  % Filter order (IIR)
    hp_filter = designfilt('highpassiir', 'FilterOrder', hp_order, ...
                           'HalfPowerFrequency', high_cutoff_freq, ...
                           'DesignMethod', 'butter');
    y_hp = filtfilt(hp_filter, y_norm);  % Apply the IIR filter with zero-phase filtering


    % Low-Pass Filter Design (IIR for efficiency)
    low_cutoff_freq = 5000 / (fs/2); % Normalized cutoff frequency (e.g., 5000 Hz)
    lp_order = 4;  % Lowered filter order for efficiency
    lp_filter = designfilt('lowpassiir', 'FilterOrder', lp_order, ...
                           'HalfPowerFrequency', low_cutoff_freq, ...
                           'DesignMethod', 'butter');
    y_filtered = filtfilt(lp_filter, y_hp);  % Zero-phase filtering

    % Parameters for pitch detection
    window_size = round(0.03 * fs); % 30 ms window
    overlap = round(0.015 * fs); % 50% overlap

    % Buffer the signal
    windowed_signal = buffer(y_filtered, window_size, overlap, 'nodelay');

    % Apply Hamming window to each segment
    hamming_window = hamming(window_size);
    windowed_signal = windowed_signal .* repmat(hamming_window, 1, size(windowed_signal, 2));

    % Number of segments (SEG)
    SEG = size(windowed_signal, 2);

    % Define pitch range (Hz)
    min_F0 = 75;  % Minimum expected F0 (Hz)
    max_F0 = 300; % Maximum expected F0 (Hz)
    min_lag = round(fs / max_F0);  % Maximum expected pitch period
    max_lag = round(fs / min_F0);  % Minimum expected pitch period

    % Pitch detection using optimized autocorrelation
    pitch_periods = zeros(1, SEG);
    N = 0;
    NUV = 0; % Number of Unvoiced Segments

    for j = 1:SEG
        segment = windowed_signal(:, j);
        [r, lags] = xcorr(segment, 'coeff');
        r = r(lags >= 0); % Consider only non-negative lags
        lags = lags(lags >= 0);

        % Focus on the lag range that corresponds to plausible pitch periods
        r = r(lags >= min_lag & lags <= max_lag);
        lags = lags(lags >= min_lag & lags <= max_lag);
    
        % Smooth the autocorrelation function
        r_smooth = smooth(r, 5);  % Moving average smoothing
    
        % Peak detection with validation
        min_peak_height = max(r_smooth) * 0.2;
        [~, locs] = findpeaks(r_smooth, 'MinPeakHeight', min_peak_height, 'MinPeakProminence', 0.1);
    
        if ~isempty(locs)
            pitch_period = lags(locs(1)); % Period corresponding to the first valid peak
        
            % Check if the detected peak might be a harmonic
            if length(locs) > 1 && (lags(locs(2)) > 1.5 * lags(locs(1)) && lags(locs(2)) < 2.5 * lags(locs(1)))
               pitch_period = lags(locs(2));  % Use the harmonic
            end
        
            pitch_periods(j) = pitch_period; % Store pitch period
            N = N + 1; % Increment the count of detected pitch periods
        else
            pitch_periods(j) = NaN; % If no valid peak is found, mark as NaN
            NUV = NUV + 1; % Increment the count of unvoiced segments
        end
    end

    
    % Remove NaN values
    pitch_periods = pitch_periods(~isnan(pitch_periods));

    % Convert pitch periods to F0 values
    F0_values = fs ./ pitch_periods;

    % Fundamental Frequency Information Measurements
    % Calculate F0 (the average F0)
    average_F0 = mean(F0_values);
    % Calculate the average pitch period (T0)
    T0 = mean(pitch_periods);
    % Calculate the highest F0
    highest_F0 = max(F0_values);
    % Calculate the lowest F0
    lowest_F0 = min(F0_values);
    % STD of F0
    F0_std = std(F0_values);
    % Calculate the Phonatory F0 Range (PFR)
    lowest_semitone = 12 * log2(lowest_F0 / 55);
    highest_semitone = 12 * log2(highest_F0 / 55);
    PFR = highest_semitone - lowest_semitone;

    % Short and Long-term Frequency Perturbation Measurements
    % Calculate Jitter Absolute (Jita)
    pitch_periods = pitch_periods(:)';
    abs_diffs = abs(diff(pitch_periods));
    Jita = sum(abs_diffs) / (N - 1);
    % Calculate Jitter in Percent (Jitt)
    Jitt = (Jita / T0) * 100;
    % Calculate Relative Average Perturbation (RAP)
    rap_sum = 0;
    for j = 2:(length(pitch_periods)-1)
        rap_sum = rap_sum + abs((pitch_periods(j-1) + pitch_periods(j) + pitch_periods(j+1)) / 3 - pitch_periods(j));
    end
    RAP = (1 / (N-2)) * (rap_sum / T0);
    % Calculate Pitch Period Perturbation Quotient (PPQ)
    ppq_sum = 0;
    for j = 1:(length(pitch_periods)-4)
        sum_5 = 0;
        for r = 0:4
            sum_5 = sum_5 + abs(pitch_periods(j+r) - pitch_periods(j+2));
        end
        ppq_sum = ppq_sum + (sum_5 / 5);
    end
    PPQ = (1 / (N-4)) * (ppq_sum / T0);
    % Calculate Smoothed Pitch Period Perturbation Quotient (sPPQ)
    smoothing_factor = 55;
    sppq_sum = 0;
    for j = 1:(length(pitch_periods) - smoothing_factor + 1)
        sum_sf = 0;
        for r = 0:(smoothing_factor - 1)
            sum_sf = sum_sf + abs(pitch_periods(j+r) - pitch_periods(j + floor(smoothing_factor / 2)));
        end
        sppq_sum = sppq_sum + (sum_sf / smoothing_factor);
    end
    sPPQ = (1 / (N - smoothing_factor + 1)) * (sppq_sum / T0);
    % Calculate Coefficient of Fundamental Frequency Variation (vF0)
    vF0 = (F0_std / average_F0) * 100;

    % Short and Long-term Amplitude Perturbation Measurements
    amp_peaks = max(windowed_signal) - min(windowed_signal);
    % Calculate Shimmer in dB (ShdB)
    shdB_sum = 0;
    for j = 1:(length(amp_peaks)-1)
        shdB_sum = shdB_sum + abs(20 * log10(amp_peaks(j+1) / amp_peaks(j)));
    end
    ShdB = (1 / (N-1)) * shdB_sum;
    % Calculate Shimmer in Percent (Shim)
    shim_sum = 0;
    for j = 1:(length(amp_peaks)-1)
        shim_sum = shim_sum + abs(amp_peaks(j+1) - amp_peaks(j));
    end
    average_amp = mean(amp_peaks);
    Shim = (shim_sum / (N-1)) / average_amp * 100;
    % Calculate Amplitude Perturbation Quotient (APQ)
    apq_sum = 0;
    for j = 1:(length(amp_peaks)-4)
        sum_5 = 0;
        for r = 0:4
            sum_5 = sum_5 + abs(amp_peaks(j+r) - amp_peaks(j+2));
        end
        apq_sum = apq_sum + (sum_5 / 5);
    end
    APQ = (1 / (N-4)) * (apq_sum / average_amp);
    % Calculate Smoothed Amplitude Perturbation Quotient (sAPQ)
    sapq_sum = 0;
    for j = 1:(length(amp_peaks) - smoothing_factor + 1)
        sum_sf = 0;
        for r = 0:(smoothing_factor - 1)
            sum_sf = sum_sf + abs(amp_peaks(j+r) - amp_peaks(j + floor(smoothing_factor / 2)));
        end
        sapq_sum = sapq_sum + (sum_sf / smoothing_factor);
    end
    sAPQ = (1 / (N - smoothing_factor + 1)) * (sapq_sum / average_amp);
    % Calculate Coefficient of Amplitude Variation (vAm)
    mean_amp = mean(amp_peaks);
    std_amp = std(amp_peaks);
    vAm = (std_amp / mean_amp) * 100;

    % Voice Break Related Measurements
    % Calculate Degree of Voice Breaks (DVB)
    threshold = 0.02;  % Adjust this value as needed
    min_break_duration = 0.05;  % Minimum duration for a voice break in seconds
    breaks = y < threshold;
    break_durations = [];
    current_break_length = 0;
    for j = 1:length(breaks)
        if breaks(j)
            current_break_length = current_break_length + 1;
        else
            if current_break_length / fs >= min_break_duration
                break_durations = [break_durations, current_break_length / fs];
            end
            current_break_length = 0;
        end
    end
    total_break_length = sum(break_durations);
    DVB = total_break_length / duration;
    % Calculate Number of Voice Breaks (NVB)
    NVB = length(break_durations);

    % Sub-Harmonic Components Related Measurements
    subharmonic_segments = 0;
    NSH = 0; % Number of Subharmonic Segments

    % Iterate only over the length of pitch_periods
    for j = 1:length(pitch_periods)
        if ~isnan(pitch_periods(j))
            if (pitch_periods(j) > T0*1.5 && pitch_periods(j) < T0*2.5) || ...
               (pitch_periods(j) > T0*2.5 && pitch_periods(j) < T0*3.5)
                subharmonic_segments = subharmonic_segments + 1;
                NSH = NSH + 1; % Increment the count of subharmonic segments
            end
        end
    end

    % Calculate Degree of Subharmonics (DSH)
    DSH = (subharmonic_segments / length(pitch_periods)) * 100;

    % Voice Irregularity Related Measurements
    % Degree of Voiceless (DUV)
    DUV = (NUV / SEG) * 100;

    % Noise Related Measurements
    % Noise-to-Harmonic Ratio (NHR)
    window_length = round(0.08192 * fs);
    windowed_signal_nhr = buffer(y_filtered, window_length, round(window_length / 2), 'nodelay');
    lp_order_nhr = 22;
    lp_cutoff_freq_nhr = 6000 / (fs / 2);
    lp_filter_nhr = fir1(lp_order_nhr, lp_cutoff_freq_nhr, 'low', hamming(lp_order_nhr + 1));
    y_lp_nhr = filter(lp_filter_nhr, 1, y_filtered);
    fs_downsampled = 12500;
    y_downsampled = resample(y_lp_nhr, fs_downsampled, fs);
    y_analytic = hilbert(y_downsampled);
    nfft = 2048;
    harmonic_energy = 0;
    inharmonic_energy = 0;
    for j = 1:size(windowed_signal_nhr, 2)
        segment = windowed_signal_nhr(:, j);
        Y = fft(segment, nfft);
        P = abs(Y).^2;
        freqs = linspace(0, fs_downsampled / 2, nfft / 2 + 1);
        harmonic_band = (freqs >= 70 & freqs <= 4500);
        inharmonic_band = (freqs >= 1500 & freqs <= 4500);
        harmonic_energy = harmonic_energy + sum(P(harmonic_band));
        inharmonic_energy = inharmonic_energy + sum(P(inharmonic_band));
    end
    NHR = (inharmonic_energy / harmonic_energy);

    % Voice Turbulence Index (VTI)
    windowed_signal_vti = buffer(y_filtered, window_length, round(window_length / 2), 'nodelay');
    harmonic_energy_vti = 0;
    inharmonic_energy_vti = 0;
    for j = 1:size(windowed_signal_vti, 2)
        segment = windowed_signal_vti(:, j);
        Y = fft(segment, nfft);
        P = abs(Y).^2;
        freqs = linspace(0, fs_downsampled / 2, nfft / 2 + 1);
        harmonic_band_vti = (freqs >= 70 & freqs <= 4500);
        inharmonic_band_vti = (freqs >= 2800 & freqs <= 5800);
        harmonic_energy_vti = harmonic_energy_vti + sum(P(harmonic_band_vti));
        inharmonic_energy_vti = inharmonic_energy_vti + sum(P(inharmonic_band_vti));
    end
    VTI = (inharmonic_energy_vti / harmonic_energy_vti);


    % Store the results in the table
    results_table.NO(i) = i;  % Sequence number
    results_table.FILE_NAME{i} = file_name;  % File name
    results_table.AVG_F0(i) = average_F0;
    results_table.AVG_T0(i) = T0;
    %results_table.HF0(i) = highest_F0;
    %results_table.LF0(i) = lowest_F0;
    results_table.STD_F0(i) = F0_std;
    results_table.PFR(i) = PFR;
    results_table.JITA(i) = Jita;
    results_table.JITT(i) = Jitt;
    results_table.RAP(i) = RAP;
    results_table.PPQ(i) = PPQ;
    results_table.SPPQ(i) = sPPQ;
    results_table.VF0(i) = vF0;
    results_table.SHDB(i) = ShdB;
    results_table.SHIM(i) = Shim;
    results_table.APQ(i) = APQ;
    results_table.SAPQ(i) = sAPQ;
    results_table.VAM(i) = vAm;
    results_table.DVB(i) = DVB;
    results_table.NVB(i) = NVB;
    results_table.DSH(i) = DSH;
    results_table.NSH(i) = NSH;
    %results_table.NUV(i) = NUV;
    %results_table.DUV(i) = DUV;
    results_table.NHR(i) = NHR;
    results_table.VTI(i) = VTI;
end

% Define the file name for the CSV output
output_csv = 'mdvp_results.csv';

% Write the table to a CSV file
writetable(results_table, output_csv);

% Display the final table to check the results
disp(results_table);
