function [peaks, relphases, relphases_unwrap] = compute_relphase(mov, subject, task, task_number, stim_plat_indices, stim, samp_freq)
% COMPUTE_RELPHASE Computes relative phase values between movement and stimulus.
%
% INPUTS:
%   mov - Array, size N-by-2, first column for the N movement values,
%         second column for corresponding time.
%   subject - String, the subject identifier.
%   task - String, the task identifier ('sync' or 'syncslow').
%   task_number - Integer, the task number.
%   stim_plat_indices - Array, size P-by-2, each row contains the index of 
%                       the stimulus at the beginning and at the end of the
%                       P-th plateau.
%   stim - Array, size S-by-1, time of the onset of the S stimuli.
%   samp_freq - Sampling frequency in Hz.
%
% OUTPUTS:
%   peaks - Array, size R-by-2, first column for the peak value, second
%           column for the time of occurence.
%   relphases - Array, size R-by-2, first column for the wrapped relative
%               phase values (bounded between -pi and +pi), second column
%               for the time of occurence. Each peak has a corresponding
%               relative phase value.
%   relphases_unwrap - Array, size R-by-2, first column for the unbounded
%                      relative phase values, second column for the time of
%                      occurence. 
%
% The function identifies flexion peaks in the movement data and computes 
% the relative phase values between these peaks and the closest stimulus.
% The relative phase is calculated as delta_T / P * (2*pi), where delta_T
% is the time difference between peak onset and stimulus onset and P is the
% time period between the two stimuli surrounding the peak.
%
% Version 1.0, aug. 2024 by XXX, tested with Matlab R2021b 
% (9.11.0.1769968) running on a maci64 computer
% Contact: XXX

peaks = [];
threshold = 65;
min_sum = 60;

switch subject
    case "IN01"
        if task == "syncslow" && task_number ~= 3, min_sum = 20; end
    case "FR04"
        if task == "sync" && task_number == 3, min_sum = 40; end
    case "IN05"
        if task == "syncslow"
            if task_number == 3, threshold = 80;
            elseif task_number == 4, threshold = 70; end
        end
    case "FR06"
        if task == "sync" && task_number == 2, min_sum = 40; end
    case "FR08"
        if task == "sync" && task_number == 1, min_sum = 40; end
    case "IN09"
        if task == "sync" && task_number == 2, min_sum = 40; end
    case "IN11", threshold = 79; min_sum = 30;
    case "IN12", threshold = 60;
    case "IN13"
        if task == "sync" && task_number == 1, min_sum = 30; end
end

peaks(:,2) = peak_identification(mov(:,1), stim_plat_indices, round(stim*500), threshold, min_sum);
peaks(:,1) = mov(peaks(:,2),1);
peaks(:,2) = peaks(:,2) / samp_freq;

if subject == "IN05" && task == "syncslow" && task_number == 3
    peaks([112, 115, 121, 135, 153], :) = [];
elseif subject == "IN05" && task == "syncslow" && task_number == 4
    peaks(116, :) = [];
end

reference = [];
reference(:,2) = stim;
reference(:,1) = 100;

% Function relphase calculates relative phases by finding the distance
% between each peak, then divide it by the distance between the two 
% surrounding stims (reference)
% It deals with cases where there are more stim than peaks and conversely
relphases_unwrap = relphase(reference, peaks);
relphases_unwrap(:,1) = - relphases_unwrap(:,1);  % Convention in the litterature

% Born between -pi and +pi
relphases = relphases_unwrap;
relphases(:,1) = (mod(relphases(:,1) + pi, 2*pi) - pi);
end


%% SUBFUNCTION PEAK_IDENTIFICATION
function peaks = peak_identification(mov, stim_plat_indices, stim, threshold, min_sum)
% PEAK_IDENTIFICATION Identifies peak flexion moments with exclusion of erroneously detected peaks.
%
% Author: XXX
%        
%
% INPUTS:
%   mov - Array, size N-by-2, first column for the N movement values,
%         second column for corresponding time.
%   stim_plat_indices - Array, size P-by-2, each row contains the index of 
%                       the stimulus at the beginning and at the end of the
%                       P-th plateau.
%   stim - Array, size S-by-1, contains the time of the onset of the S
%          stimuli.
%   threshold - Double, the threshold to detect true peaks.
%   min_sum - Double, minimum sum of absolute value of velocity between 
%            peaks to consider a peak valid.
%
% OUTPUTS:
%   peaks - Array, size P-by-1, time of occurence or the peaks.
%
% This function identifies the peaks in the movement signal, then only keep
% those who corresponding to peak flexion. Depending on the shape of the
% movement, findpeaks identifies peaks that do not corresponding to peak 
% flexion, hence a sorting is made based on a movement threshold (the 
% movement is normalized between 0 and 100 for each plateau) and the total 
% quantity of displacement between two consecutive peaks. If the identified
% peak truly is a flexion peak, then it has completed a full cycle of
% flexion-extension-flexion, hence the sum of the absolute values of
% velocity between two peaks should have at least a minimum value.

peaks = [];

% Loop through plateaus: peak detection is not the same
for plat_idx = 1:length(stim_plat_indices)

    end_plat = stim(stim_plat_indices(plat_idx, 2)) + 1;
    if plat_idx == 1
        beg_plat = stim(1);
    else
        beg_plat = stim(stim_plat_indices(plat_idx - 1, 2)) - 1;
    end

    % Normalize movement of the plateau between 0 and 100
    mov_plat = mov(beg_plat:end_plat);
    mov_plat = (mov_plat - min(mov_plat)) / (max(mov_plat - min(mov_plat))) * 100;


    % Detected peaks not necessarily peak flexion at low frequency
    %   - Use algorithm until 8th plateau
    %   - Use findpeaks afterward
    if plat_idx < 8
        % Don't take into account peaks under 65 u.a. (signal noise, false
        % start)
        [~, peaks_plateau] = findpeaks(mov_plat, "MinPeakHeight", threshold);

        rowsDel = [];
        values = zeros(numel(peaks_plateau)-1, 1);

        % DETERMINE WHETHER THE FIRST PEAK OF THE PLATEAU IS KEPT
        if plat_idx > 1
            if plat_idx == 2
                beg_prec_plat = stim(1);
            else
                beg_prec_plat = stim(stim_plat_indices(plat_idx - 2, 2)) + 1;
            end

            % Calculate index of the total variation of movement
            % If it does not exceeds a threhsold, the peak is not
            % considered a peak flexion
            sum_abs_dx = sum(abs(diff(mov(last_peak_plat + beg_prec_plat : peaks_plateau(1) + beg_plat))));
            if sum_abs_dx < 80
                rowsDel = [rowsDel; 1];
            end
        end

        % LOOP THROUGH DETECTED PEAKS TO FILTER THEM
        for peak_idx = 2 : numel(peaks_plateau)
            prec_peak = find(~isnan(peaks_plateau(1:peak_idx-1)), 1, "last");

            % Calculate index of the total variation of movement
            % If it does not exceeds a threhsold, the peak is not
            % considered a peak flexion
            dx = diff(mov_plat(peaks_plateau(prec_peak) : peaks_plateau(peak_idx)));
            sum_abs_dx = sum(abs(dx));
            values(peak_idx) = sum_abs_dx;
            if sum_abs_dx < min_sum
                peaks_plateau(peak_idx) = nan;
            end
        end

        % Suppress all the rows containing nan
        rowsDel =[rowsDel; find(isnan(peaks_plateau))];
        peaks_plateau(rowsDel) = [];

    else  % FOR PLATEAU >= 8
        [~, peaks_plateau] = findpeaks(mov_plat);
    end

    last_peak_plat = peaks_plateau(end);  % Save for comparing against first peak of the following plateau
    peaks = [peaks; peaks_plateau + beg_plat - 1];  % Add right peak indexes to the array
end

% Suppress the first peak if it is noth truly a peak
mov_before_first = mov(peaks(1)-1000 : peaks(1));
if sum(diff(abs(mov_before_first))) < 10
    peaks(1) = [];
end

% When detected twice, keep the peak only once
double_peaks = find(diff(peaks) == 0);
peaks(double_peaks) = [];

end

