function signalFiltered = ButtLowPass(sampFreq, cutoffFreq, signal)
    % buttLowPass : low pass filter of a signal with dual pass Butterworh
    %
    % Calling Sequence
    %   filteredSignal = buttLowPass(samplingFreq, cutoffFreq, signal)
    %
    % Parameters
    %   order           : number,  the order of the filter
    %   sampFreq        : number,  sampling frequency in Hz
    %   cutoffFreq      : number,  cutoff frequency in Hz
    %   signal          : vector,  the signal to filter
    %
    % Output
    %   filteredSignal  : vector,  the signal filtered, the same size as
    %                              signal (input)
    %
    % Description
    %   buttLowPass : low pass filter of a signal with dual pass Butterworh. 
    %   The filteredSignal is a vector (same size as signal)

    
    % Authors
    %   XXX
    


    
    % The signal must be a row column
    if size(signal, 1) < size(signal, 2)
        signal = signal';
    end
    

    % Prepare the parameters then create the filter coefficients b and a
    nyquistWs = sampFreq/2;                                  % We apply the Nyquist frequency
    order = 2;                                               % The order is set to 2 because of the correction
    cutoffCorrect = cutoffCorrection(cutoffFreq, sampFreq);  % Correction of the cutoffFreq
    
    [b, a] = butter(order, (cutoffCorrect ./ nyquistWs), 'Low');
    
    
    % Dual-pass filtering in the time domain 
    signalFiltered = fltsflts(signal, b, a) ; 
end

%% Subfunction cutoffCorrection
function cutoffFreqCorrected = cutoffCorrection(cutoffFreq, sampFreq)
    % cutoffCorrection : corrects the cutoff frequency for a second order
    %                    Butterworth filter
    %
    % Calling Sequence
    %   cutoffFreqCorrected = cutoffCorrection(cutoffFreq, sampFreq)
    %
    % Parameters
    %   cutoffFreq           : number,  the wanted cutoff frequency
    %   sampFreq             : number,  the sampling frequency of the
    %                                   signal to filter
    %
    % Output
    %   cutoffFreqCorrected  : number,  the cutoff frequency corrected                     
    %
    % Description
    %   cutoffFreqCorrected corrects the cutoff frequency to apply to cut
    %   the signal's frequencies at the wanted frequency. This function is
    %   used before using a 2nd order Butterworth filter
    
    % Authors
    %   Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %   Version 1.0.0 -- M. Le Guennec -- 2023-05-04
    %     First version 
    %
    % Adapted from
    %   https://www.codeproject.com/Articles/1267916/Multi-pass-Filter-Cutoff-Correction
    %
    % Source 
    %   Biomechanics and Motor Control of Human Movement, 4th-edition (page 69)


    filterPasses = 2;                            % filtfilt use 2 passes
    C = (((2^(1 / filterPasses)) - 1)^(1 / 4));  % David A. Winter butterworth correction factor
    Wc = 2*pi*cutoffFreq;                        % angular cutoff frequency
    Uc = tan(Wc/(2*sampFreq));                   % adjusted angular cutoff frequency
    Un = Uc / C;                                 % David A. Winter correction
    cutoffFreqCorrected = atan(Un)*sampFreq/pi; 
end

%% Subfunction fltsflts
function signalFiltered = fltsflts(signal, b, a)
    % fltsflts : filters the signal S using the filter coefficients b and a
    %            with dual pass
    %
    % Calling Sequence
    %   signalFiltered = fltsflts(signal, b, a)
    %
    % Parameters
    %   signal          : vector,  the input signal
    %   b               : vector,  the transfer function coefficients of
    %                              the numerator
    %   a               : vector,  the transfer function coefficients of
    %                              the denominator
    %
    % Output
    %   signalFiltered  : vector,  the signal filtered, same size as 
    %                              the  input                      
    %
    % Description
    %   fltsflts filters the signal (input) using the transfer function 
    %   coefficients a and b to filter the signal and give signalFiltered.
    %   The signal S is filterd twice, in chronological and anti-
    %   chronological directions. 
    %   This method cancels the phase lag of the causal filter and doubles 
    %   the order of the filter. 
    %   Signal is duplicated at the beginning and at the end to delete the
    %   ending and starting effect
    
    % Authors
    %   Denis Mottet - Univ. Montpellier - France
    %   Julien Lagarde - Univ. Montpellier - France
    %   Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %   Version 1.0.0 -- M. Le Guennec -- 2023-05-04
    %     Adaptation of fltsflts function from D. Mottet for Matlab
    %     https://github.com/DenisMot/ScilabDataAnalysisTemplate/blob/master/PRG/LIB_Signal/fltsflts.sci
    %     Changes in the method of duplication thanks to butterworth filter from J. Lagarde
    

    
    % Duplicate the signal before and after it to delete the starting and
    % ending effects of the filter
    duplique = [rot90(rot90(signal)) ; signal ; rot90(rot90(signal))] ;

    % Now we can apply the filter
    dupliqueFilter1 = filter(b, a, duplique);                % First filter
    dupliqueFilter1Reverse = dupliqueFilter1(end : -1 : 1);  % Reverse order
    dupliqueFilter2 = filter(b, a, dupliqueFilter1Reverse);  % Filter again (anti-chronological)
    dupliqueFilter2Reverse = dupliqueFilter2(end : -1 : 1);  % Go back to chronological order

    % Cut the added values before and after signal
    coupe = dupliqueFilter2Reverse(numel(signal)+1:2*numel(signal)) ;


    signalFiltered = coupe;
end