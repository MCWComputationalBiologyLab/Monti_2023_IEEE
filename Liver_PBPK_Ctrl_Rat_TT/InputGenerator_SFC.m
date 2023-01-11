function SF_C_R = InputGenerator_SFC(t)

% Input function parameters
%pars = [1.4569e+00,6.5277e+01,1.5797e+00]; %original from 0-60 min fitting
pars = [9.4526e-01,1.0033e+01,1.8667e+00]; %scaled 4 hr fitting

Vmax = pars(1);
Km = pars(2);
n = pars(3);

SF_C_R = Vmax*t.^n./(Km^n+t.^n);

end