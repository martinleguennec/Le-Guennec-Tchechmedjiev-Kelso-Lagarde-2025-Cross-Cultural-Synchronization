function [stim, mov] = read_file(path_file, samp_freq)
% READ_FILE Reads and processes movement data from a file.
%
% INPUTS:
%   path_file - String, the path to the data file.
%   samp_freq - Sampling frequency in Hz.
%
% OUTPUTS:
%   stim - Array, size S-by-1, voltage of the signal sent to the speakers.
%   mov - Array, size M-by-2, first column for the movement signal, second
%         column for corresponding time.
%
% The function reads movement (voltage of a goniometer) and stimulus data 
% (voltage sent to speakers acquired with acquisition card) from a text
% file, sampled at 5 kHz. Stimulus signal does not undergo any treatment.
% Movement signal is downsampled to 500 Hz, low-pass filter and normalized.
%
% Version 1.0, aug. 2024 by XXX, tested with Matlab R2021b 
% (9.11.0.1769968) running on a maci64 computer
% Contact: XXX

file_id = fopen(path_file);
if file_id == -1
    error('Unable to open file: %s', path_file);
end

data = textscan(file_id, '%f %f');
gonio = data{1};
stim = data{2};
fclose(file_id);

gonio = gonio(1 : 10 : end);
t = (0 : 1 : size(gonio, 1) - 1) ./ samp_freq;
t = t';
mov = ButtLowPass(samp_freq, 5, gonio);
mov = (mov - min(mov)) / (max(mov) - min(mov)) * 100;
mov(:,2) = t;
