function Cin = InputGenerator(t)

% Input function parameters
Cin_pars = [2.5E-2,1.36E0,6.26E-3,1.23E-1,2.75E-3,3.77E-3];

A1 = Cin_pars(1);
k1 = Cin_pars(2);
A2 = Cin_pars(3);
k2 = Cin_pars(4);
A3 = Cin_pars(5);
k3 = Cin_pars(6);

% Parameters for determining max concentration and time of occurrence
CO = 85.05;           %mL/min
TBV = 20.70;          %mL
Weight = 0.277;       %kg
Dose = 2;             %mg/kg

% Triple exponential with linear ramp
t_peak = TBV/CO; 
C_peak = Dose*Weight/TBV;

if t <= 0
    Cin = 0;
elseif t > 0 && t <= t_peak
    Slope = C_peak/t_peak;
    Cin = Slope*t; %Linear ramp
else
    Cin = A1*exp(-k1*t)+A2*exp(-k2*t)+A3*exp(-k3*t);
end

end