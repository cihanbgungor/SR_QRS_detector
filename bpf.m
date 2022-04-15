% Band-pass filter that is used to filter ECG signal + noise

% Written by Cihan Berk Gungor 
% Please cite: C. B. Güngör, P. P. Mercier, and H. Töreyin, "A Stochastic Resonance Electrocardiogram Enhancement Algorithm for Robust QRS Detection" in IEEE Journal of Body Healtcare Informatics (IEEE JBHI).

% bpf function is used in Main_SR_ECG_for_MITBIH_Arrhythmia.m file which is
% the main file to run SR-based pre-emphasis method on ECG recordings from 
% MIT-BIH Arrhythmia database 

% It is generated with built-in Fiter Designer tool of MATLAB

function b = bpf

% All frequency values are in Hz.
Fs = 360;                 % Sampling Frequency - This should be set accordingly when QT or European ST-T database is processed.
Fstop1 = 0.01;            % First Stopband Frequency
Fpass1 = 0.05;            % First Passband Frequency
Fpass2 = 100;             % Second Passband Frequency
Fstop2 = 110;             % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation = -60 dB
Dpass  = 0.057501127785;  % Passband Ripple = -24.8 dB
Dstop2 = 0.0001;          % Second Stopband Attenuation = -80 dB
flag   = 'scale';         % Sampling Flag

% Calculate the order from the parameters using KAISERORD.
[N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
    1 0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
Hd = dsp.FIRFilter( ...
    'Numerator', b);