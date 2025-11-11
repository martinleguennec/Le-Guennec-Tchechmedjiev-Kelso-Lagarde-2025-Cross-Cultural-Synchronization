function DT = dwell_time(phase_osc, threshold)
% DWELL_TIME Calculates the total time spent near attractors in phase space.
%
% INPUTS:
%   phase_osc - Array, unwrapped relative phase values.
%   threshold - Double, threshold value to determine proximity to attractors.
%
% OUTPUTS:
%   DT - Double, the total dwell time as a percentage of the total time.
%
% This function calculates the dwell time, which is the time spent near
% an attractor, indicating the coupling (all other things being equal,
% dwell time can diminushes because of decreasing concentration for e.g.).
% When there is synchronization, the relative phase is constant and hence,
% its derivative is null (Pikovsky et al., 2001). By counting the number of
% points where the derivative is below a threshold, we count the time where
% the participant was synchonized. If the frequency of the metronome
% varies, a high dwell time means that the participant was synchronized up
% to high frequencies, and hence that the coupling strength between the
% participant and the metronome was high (Arnold tongues; Pikovski et al.,
% 2001).
% _________________________________________________________________________
%  REFERENCES:
%  - Pikovsky, A., Rosenblum, M., & Kurths, J. (2003). Synchronization : A 
%    Universal Concept in Nonlinear Sciences (1st edition). Cambridge 
%    University Press.
%  - Zelic, G., Mottet, D., & Lagarde, J. (2012). Behavioral impact of 
%    unisensory and multisensory audio-tactile events: pros and cons for 
%    interlimb coordination in juggling. PLoS One, 7(2), e32308.
%       This code was initially created for this article
%
% Original version by Julien Lagarde, documented and tested aug. 2024 by

% Compute the difference of the phase oscillation
diff_phase_rel = diff(phase_osc);

% Apply median filter twice to remove peaks in the derivative
diff_phase_rel = medfilt1(medfilt1(diff_phase_rel, 3), 3);

% Filter parameters for moving average smoothing
a = 1;
b = [1/4 1/4 1/4 1/4];
aver_diff = filter(b, a, diff_phase_rel);

% Identify the points where the smoothed derivative is below the threshold
P = find(abs(aver_diff(3:end)) < threshold);

% Calculate dwell time as the percentage of points near attractors
DT = length(P) / (length(phase_osc) - 3) * 100;
