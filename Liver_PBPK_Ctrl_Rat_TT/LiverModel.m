function dy = LiverModel(t,y,p)
%% Explanation 
% SF enters via empirical input function (Cin) and either leaves via 
% outflow (C1) or is transported into hepatocytes (C3). Once in the 
% hepatocytes SF can either be converted to SFG (C4) or is transported into
% the bile (C5). SFG can then exit the hepatocyte back into the blood (C2)
% or be transported into the bile (C6)

%% Assumptions
% In this version of the model, we assume that
% 1. Once SFG is in the blood, it cannot reenter the hepatocyte
% 2. The transport of SF and SFG by transport proteins has the same Vmax
% and Km regardless of the substrate.
% 3. The OAT1B1/MRP3 coupled transporter acts as a single, unidirectional
% transporter.

% NOTE: Since SFG cannot re-enter the hepatocyte from the blood and it uses
% the same transporter as SF, the transport of SFG from the hepatocyte into
% the blood is Vmax1*C4/(Km1+C4), where C4 is the concentration of SFG in 
% the hepatocyte, NOT Vmax1*(C2-C4)/(Km1+C2+C4). 

%% State variables
C1 = y(1);     %Concenration of SF in blood
C2 = y(2);     %Concentration of conjugated SF in blood
C3 = y(3);     %Concentration of SF in hepatocytess
C4 = y(4);     %Concentration of conjugated SF in hepatocytes
C5 = y(5);     %Concentration of SF in bile duct
C6 = y(6);     %Concentration of conjugated SF in bile duct

%% Model parameters
V1 = p.V1;          %Volume of capillary blood
V2 = p.V2;          %Volume of Hepatocyte compartment
V3 = p.V3;          %Volume of Bile compartment
F = p.F;            %Flow in and out of capilary
Vmax1 = p.Vmax1;    %Max Velocity of OATPB1 transporter for transfer between capillary and hepatocytes
Km1 = p.Km1;        %MM Constant of SF for OATPB1 transporter for transfer between capillary and hepatocytes
Vmax2 = p.Vmax2;    %Max Velocity for SF conjugation within hepatocytes
Km2 = p.Km2;        %MM Constant for SF conjugation within hepatocytes
Vmax3 = p.Vmax3;    %Max Velocity of MRP2 transporter for transfer between hepatocytes and bile duct
Km3 = p.Km3;        %MM Constant for SF for MRP2 transporter for transfer between hepatocytes and bile duct
ke = p.ke;          %Clearance rate of drug from bile

% Parameters for calculating plasma concentrations
Vrbc = p.Vrbc;
Vpl = p.Vpl;
lam_max = p.lamda_max;
lam_Km = p.lamda_Km;
lam_nH = p.lamda_nH;

%% Arterial fluorescent concentration (mg/mL)
Cin = InputGenerator(t);
%Calculate plasma concentration of SF
lamda = 0.5 + lam_max*Cin.^lam_nH./(lam_Km^lam_nH+Cin.^lam_nH); % partition coefficient
Cpl_In = ((Vrbc+Vpl)./(Vrbc./lamda+Vpl)).*Cin; %concentration in plasma based on concentration in blood
%Calculate plasma concentration of SFG
One_miunus_alpha = InputGenerator_SFC(t); % ratio of SF_F to SF_C
denom = 1./(One_miunus_alpha)-1; %for reading clarity
Cin_C_pl = Cpl_In./denom;
%Calculate blood concentration of SFG
Cin_C = (Vrbc./lamda+Vpl).*Cin_C_pl/(Vrbc+Vpl);

%% Differential equations
dy = zeros(6,1);
dy(1) = (F*(Cin-C1) - Vmax1*(C1-C3)/(Km1+C1+C3))/V1; % SF dynamics in blood
dy(2) = (F*(Cin_C-C2) - Vmax1*(C2-C4)/(Km1+C2+C4))/V1; %SF conjugate dynamics io blood. Only has unidirectional flow out of hepatocyte via MRP3
dy(3) = (Vmax1*(C1-C3)/(Km1+C1+C3) - Vmax2*C3/(Km2+C3) - Vmax3*C3/(Km3+C3))/V2; % SF dynamics in hepatocytes
dy(4) = (Vmax1*(C2-C4)/(Km1+C2+C4) + Vmax2*C3/(Km2+C3) - Vmax3*C4/(Km3+C4))/V2;  % SF conjugate dynamics in hepatocytes. Only has unidirectional flow out of hepatocyte via MRP3
dy(5) = (Vmax3*C3/(Km3+C3) - ke*V3*C5)/V3;  % SF dynamics in bile duct
dy(6) = (Vmax3*C4/(Km3+C4) - ke*V3*C6)/V3;  % SF conjugate dynamics in bile duct
