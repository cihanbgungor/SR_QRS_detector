% Band-pass filter that is used to filter SR output to remove bistable
% behavior

% Written by Cihan Berk Gungor 
% Please cite: C. B. Güngör, P. P. Mercier, and H. Töreyin, "A Stochastic Resonance Electrocardiogram Enhancement Algorithm for Robust QRS Detection" in IEEE Journal of Body Healtcare Informatics (IEEE JBHI).

% hpf function is used in Main_SR_ECG_for_MITBIH_Arrhythmia.m file which is
% the main file to run SR-based pre-emphasis method on ECG recordings from 
% MIT-BIH Arrhythmia database 

% It is generated with built-in Fiter Designer tool of MATLAB


function b = hpf

% All frequency values are in Hz.
Fs = 360;                % Sampling Frequency - This should be set accordingly when QT or European ST-T database is processed.
  
Fstop = 15;              % Stopband Frequency
Fpass = 20;              % Passband Frequency
Dstop = 0.0001;          % Stopband Attenuation = -60 dB
Dpass = 0.057501127785;  % Passband Ripple = -24.8 dB
flag  = 'scale';         % Sampling Flag

% Calculate the order from the parameters using KAISERORD.
[N,Wn,BETA,TYPE] = kaiserord([Fstop Fpass]/(Fs/2), [0 1], [Dpass Dstop]);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
Hd = dsp.FIRFilter( ...
    'Numerator', b);
