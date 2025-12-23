%% Clear Workspace and Display Header
clearvars; 
disp('--------------Source------------------'); 
choice=1;
%% ----------------- Transmission Line Length -----------------
tr2 = input('Put The Length of TL(km) --> ');  % TL length   
if tr2<80
   choice=input('Do you want: (1) Transmission Line Performance Analysis or (2) Best Capacitor Compensation Configuration? Enter 1 or 2: ');
 while choice ~= 1 && choice ~= 2
    choice = input('Invalid choice. Please enter 1 or 2: ');
end
elseif tr2>240
    choice=3;
end

if choice==1
    %% ----------------- User Input for Source -----------------
f  = input('Operating frequency in (Hz) = ');   % Operating frequency

%% ----------------- User Input for Load -----------------
disp('--------------Load--------------------'); 
Vr = input('phase - phase voltage (rms) in (V) = ');  % Load voltage (phase-to-phase)
V_L=Vr;
P  = input('Active power in (watt) = ');              % Active power
pf = input('Power factor = ');                        % Power factor
while pf>1
    disp('please make sure that the power factor is less than 1')
pf = input('Power factor = '); 
end
tr = input('1-lead or 2-lag --> ');                  % Type of reactive power (lead or lag)
while tr~=1&&tr~=2
    disp('enter 1 or 2 only')
tr = input('1-lead or 2-lag --> ');   
end
%% ----------------- Reactive Power Calculation -----------------
if tr == 1
   
    Q = P * (tand(acosd(pf)));        % Capacitive reactive power for lead
    S = P - Q * j;                    % Apparent power
elseif tr == 2 
    Q = P * (tand(acosd(pf)));        % Inductive reactive power for lag
   
    S = P + Q * j;                    % Apparent power
end
V_r_phase = Vr / sqrt(3);
I_r = conj(S) / (sqrt(3) * Vr);
  % Material properties
    Eo    = 8.85*10^-12;                   % Permittivity of air
   disp('Enter the material properties:');
Ro_25 = input('Enter the resistivity of the conductor at 25°C (Ohm*m) = ');  % Resistivity at 25°C ( resistivity of Al at 25 *C =1.72 * 10^-8 ) most common

    % User Inputs
    Length  = tr2;  
    Tamb    = input('Ambiant temperature in (°C) = '); % Ambient temperature
    Sp      = input('Spiralling factor = ');           % Spiralling factor
    SKeff   = input('Skin effect factor = ');          % Skin effect factor
    T_const = input('conductor Temperature constant= '); 
    % Resistivity at desired temperature
    Ro_Tamb = Ro_25 * ((T_const + Tamb)/(T_const + 25)) * (1 + Sp) * (1 + SKeff);  

    m = input('1-phase or 3-phases --> ');            % Number of phases
while m~=1&&m~=3
    disp('enter 1 or 3 only')
 m = input('1-phase or 3-phases --> ');  
end
    %% ----------------- Geometry Calculation -----------------
    if m == 1 
        D = input('Distance between the two lines in (m) = ');  
        
        GMDeq = D; 
    elseif m == 3
        Dab = input('Distance between phase A & B in (m) = ');
        Dbc = input('Distance between phase B & C in (m) = ');
        Dca = input('Distance between phase C & A in (m) = ');
        GMDeq = nthroot(Dab * Dbc * Dca , 3);  % Geometric Mean Distance for 3-phase
    end

    %% ----------------- Bundle Conductor Parameters -----------------
    disp('From 1 to 4')
    Nb = input('Number of bundles per phase = ');
    while Nb > 4
    disp('Please Input The number of bundles between 1 and 4')
Nb = input('Number of bundles per phase = ');
end
    rb = input('Radius of one bundles in (cm) = '); 
    s  = input('Distance between bundles in (cm) = ');  
    Ns = input('Number of strands per bundles = '); 
    rs = input('Radius of one strand in (cm) = ') * 10^-2; % Convert to meters

    % GMR Calculation based on number of conductors
    if Nb == 1 
        GMRL = rb * exp(-0.25) * 10^-2; 
        GMRC = rb * 10^-2;
    elseif Nb == 2        
        GMRL = sqrt(rb * exp(-0.25) * s) * 10^-2;
        GMRC = sqrt(rb * s) * 10^-2; 
    elseif Nb == 3  
        GMRL = nthroot(rb * exp(-0.25) * s^2 , 3) * 10^-2;
        GMRC = nthroot(rb * s^2 , 3) * 10^-2;
    elseif Nb == 4   
        GMRL = nthroot(rb * exp(-0.25* s^3 * sqrt(2))^4, 4) * 10^-2;
        GMRC = nthroot(rb * s^3 * sqrt(2),4) * 10^-2;
    end

    %% ----------------- Transmission Line Parameters -----------------
    if m == 1
        L = 4e-07 * log(GMDeq/GMRL) * 10^3;   % Inductance per km
        disp(['L = ', num2str(L), ' H/km']); 
        R = (Ro_Tamb)/(pi * rs^2 * Ns * Nb) * 10^3;  % Resistance per km
        disp(['R = ', num2str(R), ' ohm/km']); 
        C = (4 * pi * Eo)/log(GMDeq/GMRC) * 10^3;    % Capacitance per km
        disp(['C = ', num2str(C), ' F/km']);
    elseif m == 3
        L = 2e-07 * log(GMDeq/GMRL) * 10^3;  
        disp(['L = ', num2str(L), ' H/km']);
        R = (Ro_Tamb)/(pi * rs^2 * Ns * Nb) * 10^3; 
        disp(['R = ', num2str(R), ' ohm/km']); 
        C = (2 * pi * Eo)/log(GMDeq/GMRC) * 10^3;
        disp(['C = ', num2str(C), ' F/km']);
    end

    %% ----------------- ABCD Constants Calculation -----------------
    Ir = conj(S/(sqrt(3) * Vr));    % Receiving current
    Z  = (R + 2 * pi * f * L * j) * Length; 
    y  = 2 * pi * f * C * j * Length;
end
%% ----------------- long Transmission Line -----------------
if tr2 > 240
 open_system('simulink_part');   
sim('simulink_part');           
end
%% ----------------- Medium Transmission Line -----------------
if 80 < tr2 && tr2 <= 240
    disp('------------Medium Transmission line--------------'); 
    
  
    tr3 = input('1-For PI Model TL or 2-For T Model TL: ');
    while tr3~=1&&tr3~=2
    disp('enter 1 or 2 only')
 tr3 = input('1-For PI Model TL or 2-For T Model TL: ');
  
end
    if tr3 == 1
        disp('------------Pi Model Transmission line--------------'); 
        A = 1 + (Z * y)/2; 
        B = Z; 
        C = y * (1 + Z*y/4);
        D = (1 + Z*y/2);
        Vs = A * Vr/sqrt(3) + B * Ir; 
        Vs_mag = abs(Vs) * sqrt(3);
        Vs_angle = angle(Vs); 
        A_s = abs(A);
        disp('------------ABCD Constants (Generalised Circuit Constants of TL)--------------'); 
        disp(['A = ', num2str(A), ' dimensionless']); 
        disp(['B = ', num2str(B), ' ohm']); 
        disp(['C = ', num2str(C), ' Siemens']); 
        disp(['D = ', num2str(D), ' dimensionless']); 
    elseif tr3 == 2
        disp('------------T Model Transmission line--------------'); 
        A = 1 + (Z * y)/2; 
        B = Z * (1 + Z*y/4); 
        C = y;
        D = 1 + (Z*y/2);
        Vs = A * Vr/sqrt(3) + B * Ir; 
        Vs_mag = abs(Vs) * sqrt(3);
        Vs_angle = angle(Vs); 
        A_s = abs(A);
        disp('------------ABCD Constants (Generalised Circuit Constants of TL)--------------'); 
        disp(['A = ', num2str(A), ' dimensionless']); 
        disp(['B = ', num2str(B), ' ohm']); 
        disp(['C = ', num2str(C), ' Siemens']); 
        disp(['D = ', num2str(D), ' dimensionless']); 
    end
end

%% ----------------- Short Transmission Line -----------------
if tr2 <= 80 &&choice==1
    disp('------------Short Transmission line--------------'); 
    
    % ABCD constants for short line (no capacitance effect)
    A = 1; 
    B = Z;  % Series impedance of the line
    C = 0; 
    D = 1;

    % Sending voltage calculation
    Vs = A * Vr/sqrt(3) + B * Ir; 
    Vs_mag = abs(Vs) * sqrt(3);       % Magnitude of sending voltage
    Vs_angle = angle(Vs);             % Angle of sending voltage

    % Display ABCD constants
    disp('------------ABCD Constants (Generalised Circuit Constants of TL)--------------'); 
    disp(['A = ', num2str(A), ' dimensionless']); 
    disp(['B = ', num2str(B), ' ohm']); 
    disp(['C = ', num2str(C), ' Siemens']); 
    disp(['D = ', num2str(D), ' dimensionless']); 
 
elseif tr2 <= 80 &&choice==2
  E = input('The maximum Vs in KV: ');  % Voltage in KV
I = input('The Maximum I in conductor in A: ');  % Current in A
Pmax = sqrt(3) * E * I / 1000;  % Maximum Power Calculation (in MW)
fprintf('The Maximum Power = %.2f MW\n', Pmax);  % Display Power in MW

vs = input('The Voltage at Sending End in KV: ');  % Sending Voltage
vr = input('The Voltage at Receiving End in KV: ');  % Receiving Voltage

PA = input('The Power Angle in degrees: ');  % Power Angle
stages = input('The Number of Stages: ');  % Number of Stages
xcs = input('The Capacitance Reactance in OHMS: ');  % Capacitance Reactance for stages
xls = input('The  Inductive Reactance per stage in OHMS: ');  % Inductive Reactance for stages
xl=input('The series Inductive Reactance  in OHMS: ');

% Effective Reactance Calculation
Xeff = (vr * vs * sind(PA)) / Pmax;
fprintf('The Effective Reactance is: %.2f OHMS\n', Xeff);

% Intermediate Reactance Calculation
xq = (xls * xcs) / (xcs - xls);

% Loop Initialization
i_closed = 0;  % Counter for open stages
j_open = stages;  % Counter for closed stages

% Variables to store the optimal Xcal and corresponding stage configuration
optimal_Xcal = inf;
optimal_i_closed = 0;

% While loop for stages
while i_closed <= stages
    Xcal = xl + (stages - i_closed) * xq - i_closed * xcs;  % Update Xcal
    
    if Xcal >= Xeff && Xcal < optimal_Xcal
        % Update optimal values if the current Xcal is valid and lower
        optimal_Xcal = Xcal;
        optimal_i_closed = i_closed;
    end
    
    i_closed = i_closed + 1;  % Update number of open stages
end

% Calculate corresponding closed and open stages for optimal Xcal
j_open = stages - optimal_i_closed;

% Final Power Calculation for optimal configuration
pmax_cal = (vr * vs * sind(PA)) / optimal_Xcal;

% Display Results
fprintf('The Optimal Number of Closed Stages = %d, The Optimal Number of Open Stages = %d\n', j_open, optimal_i_closed);
fprintf('The Minimum Xcal (still > Xeff) = %.2f OHMS, The Maximum Power calc = %.2f MW\n', optimal_Xcal, pmax_cal);


end
if choice==1
%%
%  PERFORMANCE ANALYSIS
%% ====================================================

fprintf('\n--- MODULE 3: PERFORMANCE ANALYSIS ---\n');

V_s = Vs;
I_s = C * V_r_phase + D * I_r;

Vr_nl = abs(V_s / A);
Vr_fl = abs(V_r_phase);

voltage_regulation = (Vr_nl - Vr_fl) / Vr_fl * 100;
P_r = 3 * real(V_r_phase * conj(I_r));
P_s = 3 * real(V_s * conj(I_s));
P_loss = P_s - P_r;

efficiency = (P_r / P_s) * 100;

fprintf('Voltage Regulation = %.2f %%\n', voltage_regulation);
fprintf('Efficiency          = %.2f %%\n', efficiency);
end
%% ====================================================


fprintf('\n=========== END OF PROGRAM ===========\n');