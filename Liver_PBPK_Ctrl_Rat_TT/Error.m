function Err = Error(mpar,tData,BileData)

% Simulate model
p = Parameters(mpar); %load parameters into model
y0 = [0 0 0 0 0 0]; % Initial concentrations in mg/ml 
              %(should all be zero since no drug in system initially)

options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:3), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Calculate model solution
ODE_FH = @(t,y)LiverModel(t,y,p);
sols = ode15s(ODE_FH,tData,y0,options);
y = deval(tData,sols);

% Assigning models
BileModel = y(5,:)'; %mg/ml

% Calculate SSE
Err_Bile = sum(((BileModel-BileData)/max(BileData)).^2,'omitnan');
Err = Err_Bile;
end
