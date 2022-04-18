%% The main matlab file that recordings loaded from the database and QRS enhancement done by Stochastic Resonance (SR) 

% Written by Cihan Berk Gungor 

% The main file to run SR QRS enhancement algorithm on MIT-BIH Arrhythmia
% database (Other Physionet databases may also run easily with small changes)

% Following functions are required to run this main code. 
% ud_mono_SR_for_MITBIH_Arrhythmia.m, bpf.m, and hpf.m

%% Recordings and annotations of MIT-BIH Arrhythmia database are loaded from Physionet databases

% WFDB toolbox should be installed before running this code. Instructions 
% to install it are available in https://archive.physionet.org/physiotools/matlab/wfdb-app-matlab/

% The Physionet databases other than MIT-BIH Arrhytmia database can also 
% be easily loaded by WFDB toolbox


clear all; 
clc;

db_name = 'mitdb';                                                          % db_name is set for MIT-BIH Arrhythmia Database. It should be changed for other databases
recording = '100';                                                          % Recording number is set. For each recording, this should be updated

dataset_and_rec = strcat(db_name, '/', recording);                          % String that includes database name and recording number is constructed to use rdsamp function of WFDB toolbox
[sig, Fs, tm] = rdsamp(dataset_and_rec, 1);                                 % rdsamp() is used to extract signal, time axis and Fs sampling frequency for the indicated recording
[ann,anntype,subtype,chan,num,comments]=rdann(dataset_and_rec, 'atr');      % rdann() is used to extract all information about annotations that the annotation type 'atr' is given as the second input to the function
ann_orig = ann;

% Removal of the portions of annatotions which includes not-QRS annotations % 
ann_edit = 0;
count = 1;
count2 = 1;
anntype_not_QRS = '';
for j=1:length(ann)
    if(anntype(j) ~= 'N')                                                                               
        if(anntype(j) ~= 'F')                                                                           
            if(anntype(j) ~= 'V')                                                                       
                if(anntype(j) ~= 'S')                                                                   
                    if(anntype(j) ~= 'J')                                                               
                        if(anntype(j) ~= 'a')                                                           
                            if(anntype(j) ~= 'Q')                                                       
                                if(anntype(j) ~= 'A')
                                    if(anntype(j) ~= '/')
                                        if(anntype(j) ~= 'f')
                                            if(anntype(j) ~= 'e')
                                                if(anntype(j) ~= 'R')
                                                    if(anntype(j) ~= 'j')
                                                        anntype_not_QRS(count) = anntype(j);
                                                        count = count + 1;
                                                    else
                                                        ann_edit(count2) = ann(j);
                                                        count2 = count2 + 1;
                                                    end
                                                else
                                                    ann_edit(count2) = ann(j);
                                                    count2 = count2 + 1;
                                                end
                                            else
                                                ann_edit(count2) = ann(j);
                                                count2 = count2 + 1;
                                            end
                                        else
                                            ann_edit(count2) = ann(j);
                                            count2 = count2 + 1;
                                        end
                                    else
                                        ann_edit(count2) = ann(j);
                                        count2 = count2 + 1;
                                    end
                                else
                                    ann_edit(count2) = ann(j);
                                    count2 = count2 + 1;
                                end
                            else
                                ann_edit(count2) = ann(j);
                                count2 = count2 + 1;
                            end
                        else
                            ann_edit(count2) = ann(j);
                            count2 = count2 + 1;
                        end
                    else
                        ann_edit(count2) = ann(j);
                        count2 = count2 + 1;
                    end
                else
                    ann_edit(count2) = ann(j);
                    count2 = count2 + 1;
                end
            else
                ann_edit(count2) = ann(j);
                count2 = count2 + 1;
            end
        else
            ann_edit(count2) = ann(j);
            count2 = count2 + 1;
        end   
    else
        ann_edit(count2) = ann(j);
        count2 = count2 + 1;
    end
end

ann = 0;            
ann = ann_edit;                                                             % Finally, only QRS annatations are assigned to ann array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Time = tm;                                                                  % Obtained time and signal for the recording is assigned 
u_ref = sig;                                                                % Obtained time and signal for the recording is assigned   

u_ref_nm = ((u_ref-mean(u_ref))/(max(u_ref)-min(u_ref)));                   % Normalized verison of the raw ECG signal is stored 

% ===================== Figure to show raw ECG signal ===================== %
figure;
plot(Time,u_ref)                                                            % Raw ECG signal is plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

u = u_ref;                                                                  % The raw ECG signal assigned to u variable. At this step one can add noise on the raw signal

% ================== Band-pass filtering the raw ECG signal ============== %
u=filter(bpf(),1,u);                                                        % Band-pass filter is generated
u=u(floor(length(bpf())/2):end);                                            % Band-pass filtered signal with reduced length due to bpf operation
ext = ones((length(u_ref)-length(u)), 1); ext = ext*u(length(u));           % Length of the band-pass filtered signal is corrected to be competable with time     
u = [u' ext'];                                                              % Length of the band-pass filtered signal is corrected to be competable with time    
u_bpf_out = u;                                                              % bpf output signal is stored in new variable for potential usage of it for visualization
u=u';                                                                       % Length of the band-pass filtered signal is corrected to be competable with time   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

u_nm = ((u-mean(u))/(max(u)-min(u)));                                       % Normalized verison of the band-pass filtered signal is stored

Se_srout = 0; Pp_srout = 0; F1_srout = 0;                                   % Detection performance reuslts are cleared


% =========== Parameter setup and SR solution with 4th order RK =========== %
% Parameters of the monostable well, a and b, is assigned in the following
% part.
% Step size of the solution (h) of the Langevin equation with 4th order RK
% method is assigned.
% Initial value of SR output, x_srout1 is set to zero. 
% ud_mono_SR_for_MITBIH_Arrhythmia(a,b,dr,h,u,xin1,yin1); solves Langevin 
% equation with 4th order RK method a and b are monostable well parameters, 
% dr is the damping coefficient, h step size of the solution, u is 
% signal+noise, and x1 is the initial value

% For each recording, h parameter should be set accordingly

h = 1/37.14;                                                                % For recordings 100, 103, 106, 107, 109, 115, 117, 122, 205, 208, 213, 214, 215, 219, 220, 221, 231, and 234 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/47.14; For recording 101 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/37.24; For recording 200 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/38.14; For recordings 102, 116, 118, 203, 210, 217, and 233 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/39.14; For recordings 104, 108, 121, and 123 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/50.14; For recordings 105, 207, 209, 222, 223, and 228 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/44.14; For recording 111 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/45.14; For recordings 112, 113, 114, 119, 201, 202, 212, and 230 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/40.14; For recordings 124, and 232 of the MIT-BIH Arrhthmia database
                                                                            % h = 1/50.14; For recording 207 of the MIT-BIH Arrhthmia database

a = -638;                                                                   % For all recordings of MIT-BIH Arrhythmia database

b = 0.001;                                                                  % For all recordings of MIT-BIH Arrhythmia database     

dr = 12.62;                                                                 % For all recordings of MIT-BIH Arrhythmia database 

xsrout1=0;
ysrout1=0;

x_srout_ud_m=0;                                                             % Output of SR system is initally set to zero
      
[x_srout_ud_m]=ud_mono_SR_for_MITBIH_Arrhythmia(a,b,dr,h,u,xsrout1,ysrout1);% Langevin equation is solved with 4th order Runge-Kutta method and result assigned to x_srout_ud_m

x_srout_ud_m_nm = ((x_srout_ud_m-mean(x_srout_ud_m))/(max(x_srout_ud_m)-min(x_srout_ud_m))); % SR output normalized to limit y values of it in order to get better visulization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ============= High-pass filtering applied on the SR output ============== %
x_srout_ud_m_hpf_fil=filter(hpf(),1,x_srout_ud_m);                          % High-pass filter is generated
x_srout_ud_m_hpf_fil=x_srout_ud_m_hpf_fil(floor(length(hpf())/2):end);      % High-pass filtered signal with reduced length due to hpf operation     
                                                                             
ext = ones((length(u_ref)-length(x_srout_ud_m_hpf_fil)), 1);                % Length of the band-pass filtered signal is corrected to be competable with time  
ext = ext*x_srout_ud_m_hpf_fil(length(x_srout_ud_m_hpf_fil));               % Length of the band-pass filtered signal is corrected to be competable with time   
x_srout_ud_m_hpf_fil = [x_srout_ud_m_hpf_fil ext'];                         % Length of the band-pass filtered signal is corrected to be competable with time  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_srout_ud_m_hpf_fil_nm = ((x_srout_ud_m_hpf_fil-mean(x_srout_ud_m_hpf_fil))/(max(x_srout_ud_m_hpf_fil)-min(x_srout_ud_m_hpf_fil))); % Normalized verison of the high-pass filtered SR output signal is stored

% Constant threshold is applied, performance metrics and detection performance are obtained %

in_QRS_loc = tm(ann);                                                       % Ground truth R-peak locations are assigned as a timestamps
THRESHOLD = 0.1;                                                            % Constant threshold value is set to 0.1

threshold_out_temp = 0;                                                     % Temporary array to store '0' and '1' values after threshold is applied
for i=1:length(x_srout_ud_m_hpf_fil_nm)             
    if(x_srout_ud_m_hpf_fil_nm(i) > THRESHOLD)
        threshold_out_temp(i) = 1;                                          % If a sample the high-pass filtered SR output higher than the threshold, that sample is marked with 1
    else 
        threshold_out_temp(i) = 0;                                          % If a sample the high-pass filtered SR output lower than the threshold, that sample is marked with 0
    end
end
diff_threshold_out_temp = diff(threshold_out_temp);                         % '0' to '1' and '1' to '0' transitions are identified
diff_threshold_out_temp = [diff_threshold_out_temp 0];

% The following for loop and if statements obtain the points of the middle
% point of '0' to '1' and '1' to '0' transitions. By doing this, the
% deceted QRS-wave locations are obtained
srout_QRS_counter = 1;
for i=1:length(diff_threshold_out_temp)
    if(diff_threshold_out_temp(i) == 1)
        for j=i:(i+6)
            if(diff_threshold_out_temp(j) == -1)
                srout_QRS(srout_QRS_counter) = round((i+j)/2);
                srout_QRS_counter = srout_QRS_counter + 1;
            end
        end
    end
end    

srout_QRS_loc = Time(srout_QRS);                                            % Obtained detected QRS locations are stored as timestamps

threshold_srout_line = THRESHOLD*ones(1, length(x_srout_ud_m_hpf_fil_nm));  % Constant threshold line is generated to use it in the plots for visualization

% TP and FN cases are obtained by following for loop.    
% Each value of chechformatch array indicates whether output spike
% higher than threshold is actually spike or there should be spike
% eventhough SR output won't be able to pass threshold at that
% location.
checkformatch = 0;                                                          % checkformatch is initially set to zero. 
for i=1:length(in_QRS_loc)
    checkformatch(i)=min(abs(in_QRS_loc(i)-srout_QRS_loc));
end

% FP cases are obtained by following loop.
% Each value of chechformatch2 array indicates the cases where the SR
% output spike higher than the threshold value but there is no spike in
% the ground truth.  
checkformatch2 = 0;                                                         % checkformatch2 is initially set to zero. 
for i=1:length(srout_QRS_loc)
    checkformatch2(i)=min(abs(srout_QRS_loc(i)-in_QRS_loc));
end

% Performance metrics and detection performance results of SR output are initially set to zero. 
TP_srout = 0; FN_srout = 0; FP_srout = 0; TN_srout = 0; 
Se_srout = 0; Pp_srout = 0; DER_srout = 0; Acc_srout = 0;

% Total values of TP, FN, and FP values are calculated with length()
% function. 
% +- 150 ms is the limit of the time difference between ground truth
% QRS-wave time to SR output detected QRS-wave time. If the value of the 
% checkformatch is higher than 150 ms, there sould be a detected QRS-wave 
% at the SR output but it is not. If the value of the checkformatch is 
% lower than the time difference limit of 150 ms, a detected QRS-wave at 
% the SR output is correspond to QRS-wave at the ground truth.
% If the value of the chechformatch2 is higher than the time difference
% limit of 150 ms, a detected QRS-wave at the SR output is not 
% corresponds to a QRS-wave at ground truth and it is counted as FP.

TP_srout = length(find(checkformatch<0.15));                                % Number of True Positives
FN_srout = length(find(checkformatch>0.15));                                % Number of False Negatives
FP_srout = length(find(checkformatch2>0.15));                               % Number of False Positives


Se_srout = TP_srout/(TP_srout+FN_srout)*100;                                % Sensitivity
Pp_srout = TP_srout/(TP_srout+FP_srout)*100;                                % Positive predictive value

DER_srout = 100*(FP_srout+FN_srout)/TP_srout;                               % Detection error rate
Acc_srout = 100*TP_srout/(TP_srout+FP_srout+FN_srout);                      % Accuracy 
F1_srout = 2*((Se_srout*Pp_srout)/(Se_srout+Pp_srout));                     % F1-Score

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Figure to show waveforms from different steps of the algorithm along with the true and detected QRS-waves %

figure('WindowState', 'maximized');
hold on; 
plot(Time, u_ref_nm, 'k');                                                  % Raw ECG signal (nm)          
plot(Time, u_nm, 'b');                                                      % Band-pass filtered ECG signal (nm)  
plot(Time, x_srout_ud_m_hpf_fil_nm, 'r');                                   % High-pass filtered SR output (nm)  
plot(Time, threshold_srout_line,'g');                                       % Constant threshold line                 

plot(Time(ann), 0.6, 'k.', 'MarkerSize',12);                                % Ground truth QRS-wave locations are marked 
plot(Time(srout_QRS), 0.5, 'r.', 'MarkerSize',12);                          % Detected QRS-wave locations are marked 
legend('Raw ECG (nm)', 'BPF ECG (nm)','HPF SR output (nm)', 'Threshold');
set(gca,'FontSize',15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

