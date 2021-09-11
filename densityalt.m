function h = densityalt(rho,varargin)
% DENSITYALT  Returns altitude corresponding to the given array of air densities
%   in the standard or non-standard atmosphere.
%
%   H = DENSITYALT(RHO) returns altitude, h, as a function of air density, rho.
% 
%   The input RHO can be followed by parameter/value pairs for further control
%   of DENSITYALT. Possible parameters are:
%     inputUnits     - String for units of input RHO, either kg/m³ or slug/ft³.
%                      [{SI}|kg/m3|kg/m^3  |  US|slug/ft3|slug/ft^3]
%     outputUnits    - String for units of output H, either meters or feet.
%                      [{SI}|m|meters  |  US|ft|feet]
% 
%   If the input RHO is a DimVar, inputUnits and outputUnits will be ignored and
%   the output will be a DimVar.
%
%   Example: Compute density altitude 10,000 ft. for +15C (27F) hot day.
%     h = densityalt( atmos(10000, 'tOffset', 27, 'units', 'US'), ...
%         'inputUnits', 'US', 'outputUnits', 'US')
% 
%   See also
%     ATMOS     - http://www.mathworks.com/matlabcentral/fileexchange/28135,
%     U, UNITS  - http://www.mathworks.com/matlabcentral/fileexchange/38977.
%
%   H = DENSITYALT(RHO,Param1,Val1,Param2,Val2,...)

%   Copyright 2012-2013, 2015 Sky Sartorius
%   Author contact: mathworks.com/matlabcentral/fileexchange/authors/101715


%% Parse inputs:
p = inputParser;
fName = 'densityalt';
p.FunctionName = fName;

p.addRequired('rho',@(x) validateattributes(x,{'numeric','DimVar'},...
    {'positive'},'','density, rho,'));

% Ignored if input is DimVar:
p.addParameter('inputUnits','SI');
p.addParameter('outputUnits','SI');

parse(p,rho,varargin{:});
i = p.Results;
rho = i.rho;

% Validate inputs strings:
i.inputUnits = validatestring(i.inputUnits,...
    {'SI','kg/m3','kg/m^3','US','slug/ft3','slug/ft^3'},fName,'inputUnits');
i.outputUnits = validatestring(i.outputUnits,...
    {'SI','m','meters','US','ft','feet'},fName,'outputUnits');

%% Process input density and return density in units of kg/m^3:
dimVarOut = false; % Flag to convert output to a DimVar.
if isa(rho,'DimVar')
    rho = rho/(u.kg/u.m^3);
    dimVarOut = true; 
elseif any(strcmpi(i.inputUnits,{'US','slug/ft3','slug/ft^3'}))
    % Convert from imperial units.
    rho = rho * 515.3788183931961;
end
% Otherwise, input rho is already in kg/m^3.

%% Find the altitude:

h = da_analytical_metric(rho);

%% Process output altitude from meters into desired units:
if dimVarOut
    h = h*u.m;
elseif any(strcmpi(i.outputUnits,{'US','ft','feet'})) 
    % Convert to imperial units if necessary.
    h = h / 0.3048;
end
% Otherwise, output h is already in meters.

end

function h = da_analytical_metric( rho )

%  Lapse rate Base Temp       Base Geop. Alt    Base Pressure
%   Ki (K/m)  Ti (K)          Hi (m)            P (Pa)
D =[-0.0065   288.15          0                 101325            % Troposphere
    0         216.65          11000             22632.04059693474 % Tropopause
    0.001     216.65          20000             5474.877660660026 % Stratosph. 1
    0.0028    228.65          32000             868.0158377493657 % Stratosph. 2
    0         270.65          47000             110.9057845539146 % Stratopause
    -0.0028   270.65          51000             66.938535373039073% Mesosphere 1
    -0.002    214.65          71000             3.956392754582863 % Mesosphere 2
    0         186.94590831019 84852.04584490575 .373377242877530];% Mesopause

% Constants:
rho0 = 1.225; % Sea level density, kg/m^3
g0 = 9.80665; % m/sec^2

K = D(:, 1); % K/m
T = D(:, 2); % K
H = D(:, 3); % m
P = D(:, 4); % Pa

R = P(1) / T(1) / rho0; % N-m/kg-K

% Density at base of each atmosphere layer.
rho_tab = P ./ (T * R);

% Find index of base layer

dim = size(rho);

% Initialize to layer one to handle below sea level cases.
ibase = ones(dim);
for ilayer = 1:length(rho_tab)
    msk = rho <= rho_tab( ilayer );
    ibase( msk ) = ilayer;
end

% Indexing like foo(ibase) results in a column vector when
% ibase is a vector.  The proliferation of reshape in the following code
% fights this to maintain a consistent shape whether the input is scalar,
% column vector, row vector, or matrix.

% Density ratio to base of layer
rhoonrhoi = rho ./ reshape(rho_tab(ibase), dim);

% It is easiest to calculate both the gradient and isothermal results for
% all points and then index in the correct result after the computation.

% Gradient layer
TonTi = rhoonrhoi .^ (-1.0 ./ (1.0 + g0 ./ (reshape(K(ibase), dim) * R)));
Temp = TonTi .* reshape(T(ibase), dim);
h = reshape(H(ibase), dim) + (Temp - reshape(T(ibase), dim)) ./ reshape(K(ibase), dim);

% Isothermal layer
hiso = reshape(H(ibase), dim) - log(rhoonrhoi) .* R .* reshape(T(ibase), dim) ./ g0;

% Mask to identify isothermal layers
isomask = reshape(K(ibase) == 0, dim);

% Place isothermal results in solution
h(isomask) = hiso(isomask);

end
