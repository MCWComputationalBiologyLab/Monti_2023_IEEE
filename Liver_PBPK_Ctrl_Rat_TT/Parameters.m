function p = Parameters(mpars)

% Anatomical parameters are from Hall etal (2011) and Brown etal (1997)
p.Vtot = 11.17;         %volume of total liver of rat, units: ml
p.fV1 = 0.18;           %fractional volume of liver that is capillary
p.fV2 = 0.82;           %fractional volume of liver that is hepatocytes
p.fV3 = 0.004;          %fractional volume of liver that is bile duct
p.fVrbc = 0.40;         %fractional volume of rbc in capillary
p.fVpl = 0.60;          %fractional volume of plasma in capillary
p.V1 = p.Vtot*p.fV1;    %volume of capilary; units: ml
p.V2 = p.Vtot*p.fV2;    %volume of Hepatocyte compartment; units: ml
p.V3 = p.Vtot*p.fV3;    %volume of Bile compartment; units: ml
p.Vrbc = p.fVrbc*p.V1;  %volume of red blood cells; units: ml
p.Vpl = p.fVpl*p.V1;    %volume of plasma compartment; units: ml
p.F = 11.80;            %flow in and out of capilary; units: ml/min

% Parameters defining partition coefficient for dye between rbc vs plasma
p.lamda_max = 9;        %Maximum partition coefficient 
p.lamda_Km = 5e-3;      %Km for partition coefficient
p.lamda_nH = 3;         %nH for partition coefficient

% Estimated parameter set
p.Vmax1 = mpars(1);              %Max Velocity of OATPB1 transporter for transfer between capillary and hepatocytes
p.Km1 = 9.3E-3;                  %MM Constant of SF for OATPB1 transporter for transfer between capillary and hepatocytes
p.Vmax2 = mpars(2);              %Max Velocity for SF conjugation within hepatocytes
p.Km2 = 5.8E-2;                  %MM Constant for SF conjugation within hepatocytes
p.Vmax3 = mpars(3);              %Max Velocity of MRP2 transporter for transfer between hepatocytes and bile duct
p.Km3 = 5.8E-4;                  %MM Constant for SF for MRP2 transporter for transfer between hepatocytes and bile duct
%p.ke = 1.0/2;                   %Original clearance rate of SF and SF conjugate from the bile duct
p.ke = (8.65*2.7)/(1000*p.V3);   % From Charles River poster for a 270 g rat in mL/min/100 g. 8.65 mL/min/100g is average of AM and PM sampling rates 
%Flow rate adjusted from Dr. Yang's email (works out to 5.23E-1/min). Bile flow rate is in uL/min and p.V3 is in mL 
end
