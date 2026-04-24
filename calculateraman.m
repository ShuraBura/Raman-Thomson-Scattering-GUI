function [spectrum] = calculateraman(wavelength_axis_nm, laser_nm, species, temperature_K, pressure_Pa, profile, linewidth_nm)
% CALCULATERAMAN_FIXED    Calculates a theoretical rotational Raman spectrum.
%   SPECTRUM = CALCULATERAMAN_FIXED(WAVELENGTH_AXIS_NM, LASER_NM, SPECIES, ...
%   TEMPERATURE_K, PRESSURE_PA, PROFILE, LINEWIDTH_NM) calculates the Raman 
%   spectrum for a given SPECIES ('N2', 'O2', or 'CO2').
%
%   Inputs:
%       wavelength_axis_nm - Wavelength axis for output spectrum (nm)
%       laser_nm           - Laser wavelength (nm)
%       species            - 'N2', 'O2', or 'CO2'
%       temperature_K      - Temperature (K)
%       pressure_Pa        - Pressure (Pa)
%       profile            - 'gaussian' or 'lorentzian'
%       linewidth_nm       - Instrumental linewidth (nm)
%
%   This function includes corrections for:
%       1. Accurate polarizability anisotropy values
%       2. Physically correct differential cross-section formula
%       3. Accurate partition function via summation
%       4. Proper handling of O-branch (J->J-2) and S-branch (J->J+2)
%       5. Correct Placzek-Teller coefficients for each branch

spectrum = zeros(size(wavelength_axis_nm));

% --- Physical Constants ---
c = 299792458;              % Speed of light (m/s)
h = 6.62607015e-34;         % Planck constant (J*s)
kB_J = 1.380649e-23;        % Boltzmann constant (J/K)
e0 = 8.854187817e-12;       % Permittivity of vacuum (F/m)
laser_m = laser_nm * 1e-9;  % Laser wavelength (m)

% --- Species-Specific Parameters ---
switch species
    case 'N2' % Nitrogen
        y2 = 6.31e-81;                      % Polarizability anisotropy (F^2*m^4)
        B_J = 1.98957 * c * 1e2 * h;        % Rotational constant (J) from 1.98957 cm^-1
        I = 1;                               % Nuclear spin quantum number
        g_func = @(J) 6 - (3 * mod(J, 2));  % Statistical weight (6 even, 3 odd)

    case 'O2' % Oxygen
        y2 = 1.47e-80;                      % Polarizability anisotropy (F^2*m^4)
        B_J = 1.43764 * c * 1e2 * h;        % Rotational constant (J) from 1.43764 cm^-1
        I = 0;                               % Nuclear spin quantum number
        g_func = @(J) mod(J, 2);            % Statistical weight (0 even, 1 odd)
        
    case 'CO2' % Carbon Dioxide
        y2 = 5.57e-80;                      % Polarizability anisotropy (F^2*m^4)
        B_J = 0.39021 * c * 1e2 * h;        % Rotational constant (J) from 0.39021 cm^-1
        I = 0;                               % Nuclear spin quantum number
        g_func = @(J) 1 - mod(J, 2);        % Statistical weight (1 even, 0 odd)
    
    otherwise
        error('Species must be ''N2'', ''O2'', or ''CO2''.');
end

% --- Calculate Partition Function ---
% Sum over sufficient states for convergence
Q_J_states = 0:150;
Q = sum(g_func(Q_J_states) .* (2*Q_J_states + 1) .* ...
        exp(-B_J * Q_J_states .* (Q_J_states + 1) / (kB_J * temperature_K)));

% --- Calculate Total Number Density ---
density_total = pressure_Pa / (kB_J * temperature_K);

% --- Define Line Shape Profile Function ---
if ischar(profile) && strcmpi(profile, 'gaussian')
    profilefun = @(x) exp(-x.^2);
elseif ischar(profile) && strcmpi(profile, 'lorentzian')
    profilefun = @(x) 1./(1 + x.^2);
elseif isnumeric(profile)
    % Arbitrary instrumental profile
    yy = (profile - min(profile)) / (max(profile) - min(profile));
    range = (max(wavelength_axis_nm) - min(wavelength_axis_nm)) / 2;
    xx = linspace(-range, range, length(profile));
    profilefun = @(x) interp1(xx, yy, x, 'linear', 0);
else
    error('Profile must be ''gaussian'', ''lorentzian'', or a numeric array.');
end

% --- Determine Dynamic J_max based on Boltzmann cutoff (1e-6) ---
J_test = 0:100; % Test a reasonable range
EJ_test = B_J * J_test .* (J_test + 1);
boltzmann_factor = exp(-EJ_test / (kB_J * temperature_K));
J_max_dynamic = max(J_test(boltzmann_factor > 1e-6));
if isempty(J_max_dynamic), J_max_dynamic = 0; end % Fallback if no states meet cutoff
J_max_O = max(2, J_max_dynamic); % Ensure O-branch starts at J=2
J_max_S = J_max_dynamic; % S-branch starts at J=0

% --- O-Branch: J -> J-2 (Stokes lines) ---
for J = 2:J_max_O
    K = J - 2;
    
    % Calculate transition properties
    EJ = B_J * J * (J + 1);
    deltaE = B_J * ((K^2 + K) - (J^2 + J));
    lambda_scattered_m = 1 / (1/laser_m - deltaE/(h*c));
    
    % Check if within wavelength range
    lambda_scattered_nm = lambda_scattered_m * 1e9;
    if lambda_scattered_nm < min(wavelength_axis_nm) || lambda_scattered_nm > max(wavelength_axis_nm)
        continue;
    end
    
    % Placzek-Teller coefficient for O-branch
    b_O = (3/2) * J * (J - 1) / ((2*J + 1) * (2*J - 1));
    
    % Differential cross section
    diffcross = (pi^2 / (45 * e0^2)) * b_O * y2 / (lambda_scattered_m^4);
    
    % State population
    g_J = g_func(J);
    statedensity = (density_total / Q) * g_J * (2*J + 1) * ...
                   exp(-EJ / (kB_J * temperature_K));
    
    % Add peak to spectrum
    peak_intensity = diffcross * statedensity;
    spectrum = spectrum + peak_intensity * ...
               profilefun((wavelength_axis_nm - lambda_scattered_nm) / linewidth_nm);
end

% --- S-Branch: J -> J+2 (anti-Stokes lines) ---
for J = 0:J_max_S
    K = J + 2;
    
    % Calculate transition properties
    EJ = B_J * J * (J + 1);
    deltaE = B_J * ((K^2 + K) - (J^2 + J));
    lambda_scattered_m = 1 / (1/laser_m - deltaE/(h*c));
    
    % Check if within wavelength range
    lambda_scattered_nm = lambda_scattered_m * 1e9;
    if lambda_scattered_nm < min(wavelength_axis_nm) || lambda_scattered_nm > max(wavelength_axis_nm)
        continue;
    end
    
    % Placzek-Teller coefficient for S-branch
    b_S = (3/2) * (J + 1) * (J + 2) / ((2*J + 1) * (2*J + 3));
    
    % Differential cross section
    diffcross = (pi^2 / (45 * e0^2)) * b_S * y2 / (lambda_scattered_m^4);
    
    % State population
    g_J = g_func(J);
    statedensity = (density_total / Q) * g_J * (2*J + 1) * ...
                   exp(-EJ / (kB_J * temperature_K));
    
    % Add peak to spectrum
    peak_intensity = diffcross * statedensity;
    spectrum = spectrum + peak_intensity * ...
               profilefun((wavelength_axis_nm - lambda_scattered_nm) / linewidth_nm);
end

end