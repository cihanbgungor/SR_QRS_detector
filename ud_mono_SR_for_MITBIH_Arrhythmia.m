% Function that includes impelemented 4th order Runke-Kutta method to solve
% Langevin equation with underdamped monostable well with damping
% coefficient adjustment

% Written by Cihan Berk Gungor 
% Please cite: C. B. Güngör, P. P. Mercier, and H. Töreyin, "A Stochastic Resonance Electrocardiogram Enhancement Algorithm for Robust QRS Detection" in IEEE Journal of Body Healtcare Informatics (IEEE JBHI).

% ud_mono_SR_for_MITBIH_Arrhythmia(a,b,dr_in, h,u,x1) function is used in 
% Main_SR_ECG_for_MITBIH_Arrhythmia.m file which is the main file to run 
% SR-based pre-emphasis method on ECG recordings from MIT-BIH Arrhythmia 
% database 

% Discrete solution of movement of the underdamped Brownian particle in
% monostable well. 
% a and b are monostable well parameters.
% dr_in is the initial damping coefficient that is adjusted depending on
% the input signal amplitude
% h is the step size of the 4th order Runge Kutta (RK) method solution
% u is the input signal s_noise(t)
% x1 is the initial value of the solution
% y1 is the initial value of dU(x)/dx
% x is the solution returned by function. 

% Please refer to the paper for detailed description.

function [ x ] = ud_mono_SR_for_MITBIH_Arrhythmia(a,b,dr_in,h,u,x1,y1)
    x=zeros(1,length(u));
    y=zeros(1,length(u));
    x(1)=x1;  % initial value of solution
    y(1)=y1; % initial value of dU(x)/dx
    N=length(u);
    u_pp_amp = max(u)-min(u);

    % Following for loop calculates discrete samples of x, x(i), by 4th
    % order RK method.

    % Equation that is solved here is 
    % d^2x(t)/dt^2+gamma*(dx(t)/dt)=(-dU0(x,t))/dx)+s_noise(t).  

    % dU0(x,t)/dx = -a*x(i)-b*(x(i)^3 and
    % u(t) = s_noise(t).

    for i=1:N-1
        
        % In the following if-else statements, the damping coeficient is
        % dynamically adjusted depending on the input amplitude.

        % For each recording of the MIT-BIH Arrhythmia database
        % corresponding lines should be uncommented

        % DR adjustment for reocrdings 100, 101, 102, 103, 105, 106, 107, 
        % 111, 112, 113, 115, 116, 117, 118, 119, 121, 122, 123, 124, 200, 
        % 201, 202, 203, 205, 208, 209, 212, 213, 214, 217, 219, 220, 221, 
        % 222, 223, 228, 230, 231, 232, 233, 234 of the MIT-BIH Arrhythmia 
        % database
        mul = 10;
        if u(i) < (u_pp_amp/10)
            dr = dr_in*mul;
        elseif u(i) < (u_pp_amp/5)
            dr = dr_in/100;
        else
            dr = dr_in/200;
        end

%         % DR adjustment for recording 104 of the MIT-BIH Arrhythmia database 
%         mul = 10;
%         if u(i) < (u_pp_amp/5)
%             dr = dr_in*mul;
%         else
%             dr = dr_in*8;
%         end

%         % DR adjustment for recording 207 of the MIT-BIH Arrhythmia database 
%         mul = 10;
%         if u(i) < (u_pp_amp/10)
%             dr = dr_in*mul;
%         elseif u(i) < (u_pp_amp/2)
%             dr = dr_in/100;
%         else
%             dr = dr_in/200;
%         end

%         % DR adjustment for recording 108 of the MIT-BIH Arrhythmia database
%         mul = 10;
%         if u(i) < (u_pp_amp/15)
%             dr = dr_in*mul;
%         elseif u(i) < (u_pp_amp/12)
%             dr = dr_in/mul;
%         elseif u(i) < (u_pp_amp/10)
%             dr = dr_in/10;
%         elseif u(i) < (u_pp_amp/5)
%             dr = dr_in/100;
%         else
%             dr = dr_in/200;
%         end
        
%         % DR adjustment for recording 109 of the MIT-BIH Arrhythmia database 
%         mul = 10;
%         if u(i) < (u_pp_amp/15)
%             dr = dr_in*mul;
%         elseif u(i) < (u_pp_amp/10)
%             dr = dr_in/200;
%         elseif u(i) < (u_pp_amp/5)
%             dr = dr_in/200;
%         else
%             dr = dr_in/200;
%         end

%         % DR adjustment for recording 114 of the MIT-BIH Arrhythmia database
%         mul = 10;
%         if u(i) < (u_pp_amp/10)
%             dr = dr_in*mul;
%         elseif u(i) < (u_pp_amp/5)
%             dr = dr_in/100;
%         else
%             dr = dr_in/200;
%         end   
        
%         % DR adjustment for recording 210 of the MIT-BIH Arrhythmia database
%         mul = 10;
%         if u(i) < (u_pp_amp/15)
%             dr = dr_in*mul;
%         elseif u(i) < (u_pp_amp/10)
%             dr = dr_in/mul;
%         elseif u(i) < (u_pp_amp/5)
%             dr = dr_in/100;
%         else
%             dr = dr_in/200;
%         end           
        
%         % DR adjustment for recording 215 of the MIT-BIH Arrhythmia database
%         mul = 10;
%         if u(i) < (u_pp_amp/15)
%             dr = dr_in*mul;
%         elseif u(i) < (u_pp_amp/10)
%             dr = dr_in/mul;
%         elseif u(i) < (u_pp_amp/5)
%             dr = dr_in/100;
%         else
%             dr = dr_in/200;
%         end        

        y1=y(i);
        k1=(a*x(i)-b*x(i)^3)-dr*y1+u(i);
        y2=y(i)+k1*h/2;
        k2=(a*(x(i)+y1*h/2)-b*(x(i)+y1*h/2)^3)-dr*y2+u(i);
        y3=y(i)+k2*h/2;
        k3=(a*(x(i)+y2*h/2)-b*(x(i)+y2*h/2)^3)-dr*y3+u(i+1);
        y4=y(i)+k3*h;
        k4=(a*(x(i)+y3*h)-b*(x(i)+y3*h)^3)-dr*y4+u(i+1);
        x(i+1)=x(i)+(y1+2*y2+2*y3+y4)*h/6;
        y(i+1)=y(i)+(k1+2*k2+2*k3+k4)*h/6; 
    end
end

