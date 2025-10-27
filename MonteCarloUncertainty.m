% Classical Lamination Theory with Thermal Effects - Corrected Version
% Author: Leander Tenbarge 

clear all

function [Output] = CompLT(Parameters,PlyAngles, LoadsVector, tau, deltaT)
    % Classical Lamination Theory Analysis with Thermal Effects
    %
    % INPUTS:
    %   Parameters: The Monte-Carlo parameters for each interation -> Smax,LTmax,ET1,TTmax,ET2,LCmax,TCmax,nuF,nuM,nu12,Vf,Ef,Em,thermalVector(1),thermalVector(2)
    %   elasticityMatricies: Nx4 matrix [E1, E2, V12, G12] for each ply (Pa)
    %   thermalExpansion: Nx2 matrix [alpha1, alpha2] for each ply (1/°C)
    %   PlyAngles: Nx1 vector of ply angles (degrees)
    %   LoadsVector: 6x1 vector [Nx, Ny, Nxy, Mx, My, Mxy] (N/m, N·m/m)
    %   tau: scalar, thickness of each ply (m)
    %   deltaT: scalar, temperature change from stress-free state (°C)
    %   TWP: Tsai-Wu parameters [F1; F2; F11; F22; F66; F12]
    %
    % OUTPUTS:
    %   Output structure with fields:
    %       .StrainVector: 6x1 midplane strains and curvatures
    %       .PlyStresses: structure with stress data
    %       .ABD: 6x6 stiffness matrix
    %       .abd: 6x6 compliance matrix
    %       .Q: 3x3xN stiffness matrices (local)
    %       .Qbar: 3x3xN stiffness matrices (global)
  

    % Reimporting the parameters:
    
    Smax = Parameters(1);  % Shear strength (Pa)
    LTmax = Parameters(2); % Longitudinal tensile strength (Pa)
    ET1 = Parameters(3);   % Longitudinal modulus (Pa)
    TTmax = Parameters(4); % Transverse tensile strength (Pa)
    ET2 = Parameters(5);   % Transverse modulus (Pa)
    LCmax = Parameters(6); % Longitudinal compressive strength (Pa)
    TCmax = Parameters(7); % Transverse compressive strength (Pa)
    nuF = Parameters(8);  % Fiber Poisson's ratio
    nuM = Parameters(9);  % Matrix Poisson's ratio
    nu12 = Parameters(10); % Major Poisson's ratio
    xi = 1;                 % Halpin-Tsai parameter
    Vf = Parameters(11);   % Fiber volume fraction
    Ef = Parameters(12);   % Fiber modulus (Pa)
    Em = Parameters(13);   % Matrix modulus (Pa)
    thermalVector = [Parameters(14),Parameters(15)];  % Thermal expansion coefficients (1/°C)

    
    % Calculate G12 using Halpin-Tsai
    Gm = Em / (2 * (1 + nuM));
    Gf = Ef / (2 * (1 + nuF));
    eta = ((Gf/Gm) - 1) / ((Gf/Gm) + xi);
    G12 = Gm * ((1 + xi*eta*Vf) / (1 - eta*Vf));
    
    % Setup material matrices
    propertyVector = [ET1, ET2, nu12, G12];
    elasticityMatricies = repmat(propertyVector, size(PlyAngles,1), 1);
    thermalExpansion = repmat(thermalVector,size(PlyAngles,1), 1);
    
    % Tsai-Wu parameters:
    F1 = (1/LTmax) - (1/LCmax);
    F2 = (1/TTmax) - (1/TCmax);
    F11 = 1 / (LTmax * LCmax);
    F22 = 1 / (TTmax * TCmax);
    F66 = 1 / Smax^2;
    F12 = -0.5 * sqrt(F11 * F22);
    TWP = [F1; F2; F11; F22; F66; F12];


    % Initialize
    numPlies = size(elasticityMatricies, 1);
    Output.Q = zeros(3, 3, numPlies);
    Output.Qbar = zeros(3, 3, numPlies);
    alpha_local = zeros(3, numPlies);
    alpha_global = zeros(3, numPlies);
    z = zeros(numPlies + 1, 1);
    A = zeros(3, 3);
    B = zeros(3, 3);
    D = zeros(3, 3);
    Output.ABD = zeros(6, 6);
    N_thermal = zeros(3, 1);
    M_thermal = zeros(3, 1);
    PlyAngles_rad = deg2rad(PlyAngles);
    
    
    % Calculate Q matrices (stiffness in principal directions)
    for i = 1:numPlies
        E1 = elasticityMatricies(i, 1);
        E2 = elasticityMatricies(i, 2);
        V12 = elasticityMatricies(i, 3);
        G12 = elasticityMatricies(i, 4);
        V21 = E2 * V12 / E1;
        
        Q11 = E1 / (1 - V12 * V21);
        Q12 = (V21 * E1) / (1 - V12 * V21);
        Q22 = E2 / (1 - V12 * V21);
        Q66 = G12;
        
        Output.Q(:, :, i) = [Q11, Q12, 0;
                             Q12, Q22, 0;
                             0, 0, Q66];
        
        alpha_local(:, i) = [thermalExpansion(i, 1); 
                             thermalExpansion(i, 2); 
                             0];
    end
    
    % Calculate Qbar matrices and transformed thermal expansion
    for i = 1:numPlies
        theta = PlyAngles_rad(i);
        c = cos(theta);
        s = sin(theta);
        
        Q11 = Output.Q(1, 1, i);
        Q12 = Output.Q(1, 2, i);
        Q22 = Output.Q(2, 2, i);
        Q66 = Output.Q(3, 3, i);
        
        % Transformed Q components
        Qbar11 = Q11*c^4 + 2*(Q12 + 2*Q66)*c^2*s^2 + Q22*s^4;
        Qbar12 = (Q11 + Q22 - 4*Q66)*c^2*s^2 + Q12*(c^4 + s^4);
        Qbar22 = Q11*s^4 + 2*(Q12 + 2*Q66)*c^2*s^2 + Q22*c^4;
        Qbar16 = (Q11 - Q12 - 2*Q66)*c^3*s - (Q22 - Q12 - 2*Q66)*c*s^3;
        Qbar26 = (Q11 - Q12 - 2*Q66)*c*s^3 - (Q22 - Q12 - 2*Q66)*c^3*s;
        Qbar66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*c^2*s^2 + Q66*(c^4 + s^4);
        
        Output.Qbar(:, :, i) = [Qbar11, Qbar12, Qbar16;
                                Qbar12, Qbar22, Qbar26;
                                Qbar16, Qbar26, Qbar66];
        
        % Transform thermal expansion coefficients
        T_epsilon = [c^2,      s^2,      2*s*c;
                     s^2,      c^2,     -2*s*c;
                     -s*c,     s*c,      c^2-s^2];
        
        alpha_global(:, i) = T_epsilon * alpha_local(:, i);
    end
    
    % Determine z-coordinates
    totalThickness = tau * numPlies;
    z(1) = -totalThickness / 2;
    
    for k = 2:numPlies + 1
        z(k) = z(k-1) + tau;
    end
    
    % Calculate A, B, D matrices and thermal resultants
    for i = 1:numPlies
        z_k = z(i + 1);
        z_k_minus_1 = z(i);
        
        A = A + Output.Qbar(:, :, i) * (z_k - z_k_minus_1);
        B = B + 0.5 * Output.Qbar(:, :, i) * (z_k^2 - z_k_minus_1^2);
        D = D + (1/3) * Output.Qbar(:, :, i) * (z_k^3 - z_k_minus_1^3);
        
        N_thermal = N_thermal + Output.Qbar(:, :, i) * alpha_global(:, i) * deltaT * (z_k - z_k_minus_1);
        M_thermal = M_thermal + 0.5 * Output.Qbar(:, :, i) * alpha_global(:, i) * deltaT * (z_k^2 - z_k_minus_1^2);
    end
    
    % Assemble ABD matrix
    Output.ABD(1:3, 1:3) = A;
    Output.ABD(4:6, 1:3) = B;
    Output.ABD(1:3, 4:6) = B;
    Output.ABD(4:6, 4:6) = D;
    
    % Calculate compliance matrix
    Output.abd = inv(Output.ABD);
    
    % Adjust loads for thermal effects
    LoadsVector_total = LoadsVector - [N_thermal; M_thermal];
    
    % Calculate midplane strains and curvatures
    Output.StrainVector = Output.abd * LoadsVector_total;
    epsilon0 = Output.StrainVector(1:3);
    kappa = Output.StrainVector(4:6);
    
    % Initialize stress storage
    Output.PlyStresses.global = zeros(numPlies, 3, 3);
    Output.PlyStresses.local = zeros(numPlies, 3, 3);
    Output.PlyStresses.thermal_local = zeros(numPlies, 3, 3);
    Output.PlyStresses.mechanical_local = zeros(numPlies, 3, 3);
    Output.PlyStresses.z_coords = zeros(numPlies, 3);
    
    % Calculate stresses in each ply
    for i = 1:numPlies
        theta = PlyAngles_rad(i);
        c = cos(theta);
        s = sin(theta);
        
        % Transformation matrix for stresses
        T_stress = [c^2,   s^2,    2*s*c;
                    s^2,   c^2,   -2*s*c;
                    -s*c,  s*c,    c^2-s^2];
        
        % Calculate at three locations: bottom, middle, top
        z_locations = [z(i), (z(i) + z(i+1))/2, z(i+1)];
        Output.PlyStresses.z_coords(i, :) = z_locations;
        
        for j = 1:3
            z_loc = z_locations(j);
            
            % Mechanical strains
            epsilon_mech = epsilon0 + z_loc * kappa;
            
            % Thermal strains
            epsilon_thermal = alpha_global(:, i) * deltaT;
            
            % Mechanical stresses
            sigma_mech_global = Output.Qbar(:, :, i) * epsilon_mech;
            sigma_mech_local = T_stress * sigma_mech_global;
            Output.PlyStresses.mechanical_local(i, j, :) = sigma_mech_local;
            
            % Thermal stresses
            sigma_thermal_global = Output.Qbar(:, :, i) * epsilon_thermal;
            sigma_thermal_local = T_stress * sigma_thermal_global;
            Output.PlyStresses.thermal_local(i, j, :) = sigma_thermal_local;
            
            % Total stresses
            sigma_global = sigma_mech_global + sigma_thermal_global;
            Output.PlyStresses.global(i, j, :) = sigma_global;
            sigma_local = T_stress * sigma_global;
            Output.PlyStresses.local(i, j, :) = sigma_local;
        end
    end
    
    % Calculate Tsai-Wu Failure Criterion and Factor of Safety
    for i = 1:numPlies
        for j = 1:3
            s1 = Output.PlyStresses.local(i, j, 1);
            s2 = Output.PlyStresses.local(i, j, 2);
            t12 = Output.PlyStresses.local(i, j, 3);
            
            % Tsai-Wu criterion value
            TW_value = s1*TWP(1) + s2*TWP(2) + TWP(3)*s1^2 + ...
                       TWP(4)*s2^2 + TWP(5)*t12^2 + 2*TWP(6)*s1*s2;
            
            % Store the Tsai-Wu value for reference
            Output.TW_criterion(i, j) = TW_value;
            
            
        end
    end
end



%% Structural Load Calculation Function
function [CltforceMatrix, controlPoints] = Structural(coordinates, Cp, CD, CL, conditions, tau,Tadditional)
% Overview:
% Conducts Force and Moment analysis on a thin airfoil skin for classical laminate theory derived from compressible flow relations and mechanics of materials;

% Solution Processes:
% 1.) Calculating the Chordlength
% 2.) Aerodynamic Calcs. 
% 3.) Control Point Calculations.
% 4.) Centroid Calcs.
% 5.) Mechanics of Materials.
% 6.) Force calcs.
% 7.) Displaying the results. 

% Inputs:
% coordinates: Nx2 matrix of [x, y] airfoil coordinates (m)
% Cp: Nx1 pressure coefficient distribution (dimensionless)
% CD: Drag coefficient (dimensionless)
% CL: Lift coefficient (dimensionless)
% conditions: [T_o; P_o; R; gamma; Mu; Mach]
%   T_o: Total temperature (K)
%   P_o: Total pressure (Pa)
%   R: Specific gas constant (J/(kg·K))
%   gamma: Ratio of specific heats (dimensionless)
%   Mu: Dynamic viscosity (Pa·s)
%   Mach: Mach number (dimensionless)
% tau: skin thickness (m)

% Outputs:
% CltforceMatrix: Nx6 load matrix [Nx, Ny, Nxy, Mx, My, Mxy]
%   Nx: In-plane force per unit span (N/m)
%   Ny: In-plane force per unit span (N/m)
%   Nxy: Shear force per unit span (N/m)
%   Mx: Bending moment per unit span (N/m) - NOTE: Currently returns force/span, not moment/span
%   My: Bending moment per unit span (N·m/m)
%   Mxy: Twisting moment per unit span (N·m/m)
%   controlPoints: Nx2 control point coordinates (m)


    % Initialize arrays
    nPoints = size(coordinates, 1) - 1;
    controlPoints = zeros(nPoints, 2);
    controlLength = zeros(nPoints, 1);
    controlAngles = zeros(nPoints, 1);
    controlArea = zeros(nPoints, 1);
    controlNormals = zeros(nPoints, 1);
    Xsum = zeros(nPoints, 1);
    Ysum = zeros(nPoints, 1);
    areavector = zeros(nPoints, 1);
    CltforceMatrix = zeros(nPoints, 6);
    

    %% 1.) Calculating the Chordlength
    % Find midpoint index
    midIdx = ceil(size(coordinates, 1) / 2);
    
    % Get leading edge and trailing edge points
    leadingEdge = coordinates(midIdx, :);  % LE typically at midpoint
    trailingEdge = coordinates(1, :);       % TE at start (or end)
    
    % Calculate chord as distance between LE and TE
    chord = sqrt((trailingEdge(1) - leadingEdge(1))^2 + ...
                 (trailingEdge(2) - leadingEdge(2))^2);

    %% 2.) Aerodynamic Calcs. 
    Pinf = conditions(2) / (1 + ((conditions(4) - 1)/2) * conditions(6)^2)^(conditions(4) / (conditions(4) - 1));
    Tinf = conditions(1) / (1 + ((conditions(4) - 1)/2) * conditions(6)^2);
    rhoinf = Pinf / (conditions(3) * Tinf);
    Vinf = conditions(6) * sqrt(conditions(4) * conditions(3) * Tinf);
    Q = 0.5 * rhoinf * (Vinf^2);
    Re = rhoinf*Vinf*chord/conditions(5);
    
    %% 3.) Control Point Calculations.
    for i = 1:nPoints
        controlPoints(i, 1) = (coordinates(i, 1) + coordinates(i+1, 1)) / 2;
        controlPoints(i, 2) = (coordinates(i, 2) + coordinates(i+1, 2)) / 2;
        controlLength(i) = sqrt((coordinates(i+1, 1) - coordinates(i, 1))^2 + (coordinates(i+1, 2) - coordinates(i, 2))^2);
        controlAngles(i) = atan2(coordinates(i+1, 2) - coordinates(i, 2), coordinates(i+1, 1) - coordinates(i, 1));
        controlNormals(i) = controlAngles(i) + pi/2;
        controlArea(i) = controlLength(i) * tau;
        Xsum(i) = controlArea(i) * controlPoints(i, 1);
        Ysum(i) = controlArea(i) * controlPoints(i, 2);
    end
    
    %% 4.) Centroid Calcs.
    centroid = [sum(Xsum) / sum(controlArea), sum(Ysum) / sum(controlArea)];
    
    %% 5.) Mechanics of Materials.
    Ixx = 0;
    Iyy = 0;
    for i = 1:nPoints
        % Area using shoelace formula
        areavector(i) = coordinates(i, 1) * coordinates(i+1, 2) - coordinates(i+1, 1) * coordinates(i, 2);
        
        % Second moments of area
        Ixx_local = (1/12) * tau * controlLength(i)^3 * cos(controlAngles(i))^2 + controlArea(i) * (controlPoints(i, 2) - centroid(2))^2;
        Iyy_local = (1/12) * tau * controlLength(i)^3 * sin(controlAngles(i))^2 + controlArea(i) * (controlPoints(i, 1) - centroid(1))^2;
        
        % Accumulate
        Ixx = Ixx + Ixx_local;
        Iyy = Iyy + Iyy_local;
    end
    
    % Calculate moment coefficient about structural centroid from Cp distribution
    Cm_centroid = 0;
    for i = 1:nPoints
        % Average Cp at this panel
        Cp_avg = (Cp(i) + Cp(i+1)) / 2;
        
        % Non-dimensional moment arm components
        x_norm = (controlPoints(i, 1) - centroid(1)) / chord;
        y_norm = (controlPoints(i, 2) - centroid(2)) / chord;
        dL_norm = controlLength(i) / chord;
        
        % Moment coefficient increment (perpendicular moment arm)
        dCm = Cp_avg * dL_norm * (x_norm * sin(controlNormals(i)) - y_norm * cos(controlNormals(i)));
        
        Cm_centroid = Cm_centroid + dCm;
    end
    
    % Calculate torsional moment and shear flow
    M_torsion = Cm_centroid * Q * chord^2;  % N·m/m (torsional moment per unit span)
    areaMean = 0.5 * abs(sum(areavector));  % Enclosed area
    shearflow = (M_torsion + Tadditional) / (2 * areaMean); % N/m (shear flow)
    
    
    
    %% 6.) Force calcs.
    for i = 1:nPoints

        % Local bending moment per unit span from pressure
        Cp_avg = (Cp(i) + Cp(i+1)) / 2;
        MxLocal = -(Cp_avg * Q * controlLength(i) );
        
        % Force resultants (Classical Laminate Theory format)
        CltforceMatrix(i, 1) = 0;                         % Nx (N/m)
        CltforceMatrix(i, 2) = 0;                         % Ny (N/m)
        CltforceMatrix(i, 3) = shearflow;                 % Nxy (N/m)
        CltforceMatrix(i, 4) = MxLocal;                   % Mx (N·m/m)
        CltforceMatrix(i, 5) = 0;                         % My (N·m/m)
        CltforceMatrix(i, 6) = 0;                         % Mxy (N·m/m)
    end

    %% 7.) Displaying the results

    % Header:
    fprintf('\n\nStructural Analysis of a  Composite Wing Section for Classical Lamination Theory\n')
    fprintf('Author: Leander Tenbarge\n\n')

    % Input Conditions:
    fprintf('\n=== Input Conditions ===================\n');
    fprintf('T₀     = %.4f K (Stagnation Temperature)\n', conditions(1));
    fprintf('P₀     = %.4f Pa (Stagnation Pressure)\n', conditions(2));
    fprintf('R      = %.4f J/(kg·K) (Gas Constant)\n', conditions(3));
    fprintf('γ      = %.4f (Specific Heat Ratio)\n', conditions(4));
    fprintf('μ      = %.4e Pa·s (Dynamic Viscosity)\n', conditions(5));
    fprintf('L      = %.4e m (Chordlength)\n', chord);
    fprintf('Mach   = %.4f (Mach Number)\n', conditions(6));
    fprintf('=======================================\n\n');
    
    % Freestream conditions:
    fprintf('\n=== Freestream Conditions =============\n');
    fprintf('P∞     = %.4f Pa\n', Pinf);
    fprintf('T∞     = %.4f K\n', Tinf);
    fprintf('ρ∞     = %.4f kg/m³\n', rhoinf);
    fprintf('V∞     = %.4f m/s\n', Vinf);
    fprintf('Q∞     = %.4f Pa\n', Q);
    fprintf('Re     = %.4e (Reynolds Number)\n', Re);
    fprintf('=======================================\n\n');
    
    % Structural properties:
    fprintf('\n=== Structural Properties =============\n');
    fprintf('Cx     = %.4f m (Centroid X Position)\n', centroid(1));
    fprintf('Cy     = %.4f m (Centroid Y Position)\n', centroid(2));
    fprintf('Ixx    = %.4e m^4 (Moment of Inertia about x-axis)\n', Ixx);
    fprintf('Iyy    = %.4e m^4 (Moment of Inertia about y-axis)\n', Iyy);
    fprintf('A_enc  = %.4e m^2 (Enclosed Area)\n', areaMean);
    fprintf('=======================================\n\n');
    
    % Aerodynamic Loads
    fprintf('\n=== Aerodynamic Loads ================\n');
    fprintf('CL     = %.4f (Lift Coefficient)\n', CL);
    fprintf('CD     = %.4f (Drag Coefficient)\n', CD);
    fprintf('Cm     = %.4f (Moment Coefficient about Centroid)\n', Cm_centroid);
    fprintf('L      = %.4f N/m (Lift per unit span)\n', CL * Q * chord);
    fprintf('D      = %.4f N/m (Drag per unit span)\n', CD * Q * chord);
    fprintf('M_tor  = %.4f N·m/m (Torsional Moment)\n', M_torsion);
    fprintf('q      = %.4f N/m (Shear Flow)\n', shearflow);
    fprintf('======================================\n\n');
end


















%% Structural Analysis:

% Importing the and processing the structural data:
structuralInput = readtable("StructuralData.xlsx");

% Extract structural inputs:
Xcoords = table2array(structuralInput(:, 1));
Ycoords = table2array(structuralInput(:, 2));
Cpdist = table2array(structuralInput(:, 3));
Cl = table2array(structuralInput(1, 4));
Cd = table2array(structuralInput(1, 5));
Conditions = table2array(structuralInput(1:6, 6));
Tau = table2array(structuralInput(1, 7));
NumPlies = table2array(structuralInput(2, 7));
Coordinates = [Xcoords, Ycoords];
Tadditional = 0;

% Perform structural analysis:
[CltforceMatrix, controlPoints] = Structural(Coordinates, Cpdist, Cd, Cl, Conditions, Tau,Tadditional);




%% Monte-Carlo Sensitivity Analysis setup:

% Inputs
N = 1000; % The Number of Iterations
PlyAngles = [0, 45, 90, 45, 0]; % Gives us the number of plys and the layup
deltaT = 0; % Change in temperature


% Importing and processing the Material data (Means for Monte-Carlo):
materialInput = readtable('MaterialData.xlsx');

% Extracting the Material Properties and defining the means:
Smax =table2array(materialInput(1,1));  % Shear strength (Pa)
LTmax = table2array(materialInput(1,2)); % Longitudinal tensile strength (Pa)
ET1 = table2array(materialInput(1,3));   % Longitudinal modulus (Pa)
TTmax = table2array(materialInput(1,4)); % Transverse tensile strength (Pa)
ET2 = table2array(materialInput(1,5));   % Transverse modulus (Pa)
LCmax = table2array(materialInput(1,6)); % Longitudinal compressive strength (Pa)
EC1 = table2array(materialInput(1,7));   % Longitudinal compressive modulus (Pa) % not used in the Analysis ***
TCmax = table2array(materialInput(1,8)); % Transverse compressive strength (Pa)
EC2 = table2array(materialInput(1,9));   % Transverse compressive modulus (Pa)   % not used in this Analysis ***
nuF = table2array(materialInput(1,10));  % Fiber Poisson's ratio
nuM = table2array(materialInput(1,11));  % Matrix Poisson's ratio
nu12 = table2array(materialInput(1,12)); % Major Poisson's ratio
xi = table2array(materialInput(1,13));   % Halpin-Tsai parameter % not used in this Analysis ***
Vf = table2array(materialInput(1,14));   % Fiber volume fraction
Ef = table2array(materialInput(1,15));   % Fiber modulus (Pa)
Em = table2array(materialInput(1,16));   % Matrix modulus (Pa)
thermalVector = [(materialInput.Thermal(1)),(materialInput.Thermal(2))];  % Thermal expansion coefficients (1/°C)


% Creating the parameter mean values:
Pmean = [Smax,LTmax,ET1,TTmax,ET2,LCmax,TCmax,nuF,nuM,nu12,Vf,Ef,Em,thermalVector(1),thermalVector(2)];

% Coefficient of variation values (Example values, refer to manufacturer data):
Pcv = [
    0.15;   % Smax - Shear strength
    0.15;   % LTmax - Longitudinal tensile strength
    0.08;   % ET1 - Longitudinal modulus
    0.20;   % TTmax - Transverse tensile strength
    0.12;   % ET2 - Transverse modulus
    0.18;   % LCmax - Longitudinal compressive strength
    0.20;   % TCmax - Transverse compressive strength
    0.03;   % nuF - Fiber Poisson's ratio
    0.05;   % nuM - Matrix Poisson's ratio
    0.06;   % nu12 - Major Poisson's ratio
    0.06;   % Vf - Fiber volume fraction
    0.04;   % Ef - Fiber modulus
    0.08;   % Em - Matrix modulus
    0.10;   % thermalVector(1) - Longitudinal CTE (alpha1), With the Correction factor for numerical issues
    0.15;   % thermalVector(2) - Transverse CTE (alpha2), With the Correction factor for numerical issues
];



% Generating the Distributions:
uncertainParams = struct();
uncertainParams.Smax = struct('mean',Pmean(1),"CV",Pcv(1));
uncertainParams.LTmax = struct('mean',Pmean(2),"CV",Pcv(2));
uncertainParams.ET1 = struct('mean',Pmean(3),"CV",Pcv(3));
uncertainParams.TTmax = struct('mean',Pmean(4),"CV",Pcv(4));
uncertainParams.ET2 = struct('mean',Pmean(5),"CV",Pcv(5));
uncertainParams.LCmax = struct('mean',Pmean(6),"CV",Pcv(6));
uncertainParams.TCmax = struct('mean',Pmean(7),"CV",Pcv(7));
uncertainParams.nuF = struct('mean',Pmean(8),"CV",Pcv(8));
uncertainParams.nuM = struct('mean',Pmean(9),"CV",Pcv(9));
uncertainParams.nu12 = struct('mean',Pmean(10),"CV",Pcv(10));
uncertainParams.Vf = struct('mean',Pmean(11),"CV",Pcv(11));
uncertainParams.Ef = struct('mean',Pmean(12),"CV",Pcv(12));
uncertainParams.Em = struct('mean',Pmean(13),"CV",Pcv(13));
uncertainParams.TV1 = struct('mean',Pmean(14),"CV",Pcv(14));
uncertainParams.TV2 = struct('mean',Pmean(15),"CV",Pcv(15));

% Generating the LHS samples:
pNames = fieldnames(uncertainParams);
nParams = length(pNames);
LHSsamples = lhsdesign(N,nParams);

% Generating the Distribution:
dist = zeros(N, nParams);

% Transform LHS samples to actual distributions
for i = 1:nParams
    paramName = pNames{i};
    
    % Extract mean and CV
    mu = uncertainParams.(paramName).mean;
    cv = uncertainParams.(paramName).CV;
    
    % Calculate standard deviation from CV
    sigma = mu * cv;
    
    % Transform from uniform [0,1] to normal distribution
    dist(:, i) = norminv(LHSsamples(:, i), mu, sigma);
end

%% Display Monte Carlo Sensitivity Analysis Setup

fprintf('\n');
fprintf('========================================\n');
fprintf('  MONTE CARLO SENSITIVITY ANALYSIS\n');
fprintf('========================================\n\n');

% Basic Setup Information
fprintf('--- Simulation Parameters ---\n');
fprintf('Number of Iterations: %d\n', N);
fprintf('Number of Parameters: %d\n', nParams);
fprintf('Ply Layup: [%s]\n', num2str(PlyAngles));
fprintf('Temperature Change: %.1f °C\n\n', deltaT);

% Parameter Statistics Table
fprintf('--- Parameter Distributions ---\n');
fprintf('%-8s | %-10s | %-12s | %-12s | %-12s\n', ...
    'Param', 'Mean', 'Std Dev', 'CV', 'Unit');
fprintf('---------|------------|--------------|--------------|------------------\n');

paramUnits = {
    'Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Pa', 'Pa', ...  % Strengths and moduli
    '-', '-', '-', '-', 'Pa', 'Pa', '1/°C', '1/°C'  % Ratios and thermal
};

for i = 1:nParams
    paramName = pNames{i};
    mu = uncertainParams.(paramName).mean;
    cv = uncertainParams.(paramName).CV;
    sigma = mu * cv;
    
    % Format output based on magnitude
    if abs(mu) >= 1e6 || abs(mu) <= 1e-3
        fprintf('%-8s | %10.3e | %12.3e | %12.2f%% | %-16s\n', ...
            paramName, mu, sigma, cv*100, paramUnits{i});
    else
        fprintf('%-8s | %10.4f | %12.4f | %12.2f%% | %-16s\n', ...
            paramName, mu, sigma, cv*100, paramUnits{i});
    end
end

fprintf('\n--- Distribution Summary ---\n');
fprintf('Sampling Method: Latin Hypercube Sampling (LHS)\n');
fprintf('Distribution Type: Normal (Gaussian)\n');
fprintf('Distribution Matrix Size: %d x %d\n\n', size(dist, 1), size(dist, 2));
fprintf('\n========================================\n\n');

%% Running the Monte Carlo Analysis

% Initialize storage structure
Data = struct('Airfoil', cell(N, 1));
TWA = zeros(N, size(Coordinates,1)-1);  % Average TW per iteration
TWM = zeros(N, 1);                       % Maximum TW per iteration

% Run Monte Carlo loop
for i = 1:N
    for j = 1:size(Coordinates,1)-1
        Data(i).Airfoil(j) = CompLT(dist(i,:), PlyAngles, CltforceMatrix(j, :)', Tau/NumPlies, deltaT);
        
        % Store maximum TW criterion value for this location
        TWA(i,j) = max(Data(i).Airfoil(j).TW_criterion(:));
    end
    % Store maximum TW across all locations for this iteration
    TWM(i) = max(TWA(i,:));
end



%% Post-Processing: Reliability Analysis and Sensitivity Indices
% Author: Based on Leander Tenbarge's Monte Carlo Analysis

%% Calculate Reliability Metrics
failureThreshold = 1.0;
failureIndicator = TWM >= failureThreshold;
nFailures = sum(failureIndicator);
Pf = nFailures / N;
R = 1 - Pf;

if Pf > 0 && Pf < 1
    beta = -norminv(Pf);
else
    beta = NaN;
end

%% Control Point Analysis
nControlPts = size(TWA, 2);
cpMean = mean(TWA, 1);
cpStd = std(TWA, 1);
cpFailRate = sum(TWA >= failureThreshold, 1) / N * 100;

%% Compute Sobol Sensitivity Indices
S_first = zeros(nParams, 1);
V_total = var(TWM);

for i = 1:nParams
    paramValues = dist(:, i);
    nBins = 25;
    edges = linspace(min(paramValues), max(paramValues), nBins+1);
    [~, binIdx] = histc(paramValues, edges);
    
    E_Y_given_Xi = zeros(nBins, 1);
    bin_counts = zeros(nBins, 1);
    
    for b = 1:nBins
        idx = binIdx == b;
        if sum(idx) > 0
            E_Y_given_Xi(b) = mean(TWM(idx));
            bin_counts(b) = sum(idx);
        end
    end
    
    valid_bins = bin_counts > 0;
    if sum(valid_bins) > 1
        weights = bin_counts(valid_bins) / N;
        V_Xi = sum(weights .* (E_Y_given_Xi(valid_bins) - mean(TWM)).^2);
        S_first(i) = max(0, min(1, V_Xi / V_total));
    end
end

S_first_normalized = S_first / sum(S_first);
[S_sorted, sortIdx] = sort(S_first_normalized, 'descend');
pNames_sorted = pNames(sortIdx);

%% Summary Statistics
fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
fprintf('║         MONTE CARLO RELIABILITY ANALYSIS RESULTS           ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

fprintf('  Simulations: %d  |  Failures: %d  |  Pf: %.3f%%  |  β: %.2f\n', ...
    N, nFailures, Pf*100, beta);
fprintf('  Mean TW: %.3f  |  Std: %.3f  |  95th%%ile: %.3f\n\n', ...
    mean(TWM), std(TWM), prctile(TWM, 95));

fprintf('  Top 5 Parameters (Sobol Index):\n');
for i = 1:5
    fprintf('    %d. %-9s: %5.2f%%\n', i, pNames_sorted{i}, S_sorted(i)*100);
end
fprintf('\n════════════════════════════════════════════════════════════\n\n');

%% Create Tabbed Figure
fig = figure('Position', [50, 50, 1600, 900], 'Color', 'w', 'Name', 'Monte Carlo Analysis');
tg = uitabgroup(fig);

%% ========================================================================
%% TAB 1: RELIABILITY ANALYSIS
%% ========================================================================
tab1 = uitab(tg, 'Title', 'Reliability Analysis');

% Plot 1: Control Point Variability
ax1 = subplot(2, 2, [1, 3], 'Parent', tab1);
boxplot(ax1, TWA, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol', '');
hold(ax1, 'on');
plot(ax1, 1:nControlPts, cpMean, 'r-', 'LineWidth', 2);
yline(ax1, failureThreshold, 'r--', 'LineWidth', 2.5, 'Alpha', 0.7);
xlabel(ax1, 'Control Point Index', 'FontSize', 11);
ylabel(ax1, 'Tsai-Wu Criterion', 'FontSize', 11);
title(ax1, 'Control Point Variability Analysis', 'FontSize', 13, 'FontWeight', 'bold');
legend(ax1, 'Mean TW', 'Failure Threshold', 'Location', 'best');
grid(ax1, 'on');
grid(ax1, 'minor');
ylim(ax1, [0, max(max(TWA))*1.1]);

% Plot 2: Failure Distribution with Normal Overlay
ax2 = subplot(2, 2, 2, 'Parent', tab1);
histogram(ax2, TWM, 40, 'Normalization', 'pdf', 'FaceColor', [0.2, 0.5, 0.8], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold(ax2, 'on');

% Overlay normal distribution
mu_TW = mean(TWM);
sigma_TW = std(TWM);
x_norm = linspace(min(TWM), max(TWM), 200);
y_norm = normpdf(x_norm, mu_TW, sigma_TW);
plot(ax2, x_norm, y_norm, 'k--', 'LineWidth', 2, 'DisplayName', 'Normal Fit');

xline(ax2, failureThreshold, 'r-', 'LineWidth', 3, 'Label', 'Failure', ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
xline(ax2, mean(TWM), 'g-', 'LineWidth', 2, 'Label', 'Mean', ...
    'LabelHorizontalAlignment', 'right', 'FontSize', 10);
xlabel(ax2, 'Max Tsai-Wu Value', 'FontSize', 11);
ylabel(ax2, 'Probability Density', 'FontSize', 11);
title(ax2, 'Failure Distribution', 'FontSize', 12, 'FontWeight', 'bold');
legend(ax2, 'Actual', 'Normal Fit', 'Location', 'best');
grid(ax2, 'on');

% Plot 3: CDF
ax3 = subplot(2, 2, 4, 'Parent', tab1);
[f, x] = ecdf(TWM);
plot(ax3, x, f, 'b-', 'LineWidth', 2.5);
hold(ax3, 'on');
xline(ax3, failureThreshold, 'r--', 'LineWidth', 2.5);
fill(ax3, [min(x), x(x<=failureThreshold)', failureThreshold], ...
    [0, f(x<=failureThreshold)', f(find(x<=failureThreshold, 1, 'last'))], ...
    'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill(ax3, [failureThreshold, x(x>=failureThreshold)', max(x)], ...
    [f(find(x>=failureThreshold, 1, 'first')), f(x>=failureThreshold)', 1], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel(ax3, 'Tsai-Wu Value', 'FontSize', 11);
ylabel(ax3, 'Cumulative Probability', 'FontSize', 11);
title(ax3, sprintf('Cumilative distribution function: R = %.1f%%, Pf = %.1f%%', R*100, Pf*100), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid(ax3, 'on');
legend(ax3, 'CDF', 'Threshold', sprintf('Safe (%.1f%%)', R*100), ...
    sprintf('Failed (%.1f%%)', Pf*100), 'Location', 'best');

%% ========================================================================
%% TAB 2: AIRFOIL SPATIAL ANALYSIS
%% ========================================================================
tab2 = uitab(tg, 'Title', 'Airfoil Spatial Analysis');

% Airfoil Uncertainty Map
ax4 = axes('Parent', tab2);
scatter(ax4, controlPoints(:,1), controlPoints(:,2), 200, cpStd, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
hold(ax4, 'on');
plot(ax4, Coordinates(:,1), Coordinates(:,2), 'k-', 'LineWidth', 2);
colormap(ax4, hot);
cb = colorbar(ax4);
ylabel(cb, 'Std Dev of TW', 'FontSize', 11);
axis(ax4, 'equal');
xlabel(ax4, 'X (m)', 'FontSize', 12);
ylabel(ax4, 'Y (m)', 'FontSize', 12);
title(ax4, 'Airfoil Uncertainty Distribution', 'FontSize', 14, 'FontWeight', 'bold');
grid(ax4, 'on');

% Add text annotation with statistics
text(ax4, 0.02, 0.98, sprintf('Mean TW: %.3f\nMax Std Dev: %.3f\nMax Fail Rate: %.1f%%', ...
    mean(cpMean), max(cpStd), max(cpFailRate)), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

%% ========================================================================
%% TAB 3: SOBOL SENSITIVITY ANALYSIS
%% ========================================================================
tab3 = uitab(tg, 'Title', 'Sobol Sensitivity Analysis');

% Plot 1: Parameter Sensitivity Rankings
ax5 = subplot(2, 2, 1, 'Parent', tab3);
topN = min(12, nParams);
y_pos = 1:topN;
barh(ax5, y_pos, S_sorted(1:topN)*100, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'k');
set(ax5, 'YTick', y_pos, 'YTickLabel', pNames_sorted(1:topN), 'YDir', 'reverse');
xlabel(ax5, 'Sobol Index (%)', 'FontSize', 11);
title(ax5, 'Parameter Sensitivity Ranking', 'FontSize', 13, 'FontWeight', 'bold');
grid(ax5, 'on');
xlim(ax5, [0, max(S_sorted(1:topN)*100)*1.15]);

% Plot 2: Most Influential Parameter Distribution
ax6 = subplot(2, 2, 2, 'Parent', tab3);
[~, mostImpIdx] = max(S_first_normalized);
scatter(ax6, dist(:, mostImpIdx), TWM, 25, TWM, 'filled', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
hold(ax6, 'on');
yline(ax6, failureThreshold, 'r--', 'LineWidth', 2.5);
xlabel(ax6, pNames{mostImpIdx}, 'FontSize', 11, 'Interpreter', 'none');
ylabel(ax6, 'Tsai-Wu Criterion', 'FontSize', 11);
title(ax6, ['Most Influential: ', pNames{mostImpIdx}], 'FontSize', 12, 'FontWeight', 'bold');
colormap(ax6, jet);
cb = colorbar(ax6);
ylabel(cb, 'TW Value', 'FontSize', 10);
grid(ax6, 'on');

% Plot 3: Second Most Influential Parameter Distribution
ax7 = subplot(2, 2, 3, 'Parent', tab3);
if nParams >= 2
    secondImpIdx = sortIdx(2);
    scatter(ax7, dist(:, secondImpIdx), TWM, 25, TWM, 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
    hold(ax7, 'on');
    yline(ax7, failureThreshold, 'r--', 'LineWidth', 2.5);
    xlabel(ax7, pNames{secondImpIdx}, 'FontSize', 11, 'Interpreter', 'none');
    ylabel(ax7, 'Tsai-Wu Criterion', 'FontSize', 11);
    title(ax7, ['2nd Most Influential: ', pNames{secondImpIdx}], 'FontSize', 12, 'FontWeight', 'bold');
    colormap(ax7, jet);
    cb = colorbar(ax7);
    ylabel(cb, 'TW Value', 'FontSize', 10);
    grid(ax7, 'on');
end

% Plot 4: Third Most Influential Parameter Distribution
ax8 = subplot(2, 2, 4, 'Parent', tab3);
if nParams >= 3
    thirdImpIdx = sortIdx(3);
    scatter(ax8, dist(:, thirdImpIdx), TWM, 25, TWM, 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
    hold(ax8, 'on');
    yline(ax8, failureThreshold, 'r--', 'LineWidth', 2.5);
    xlabel(ax8, pNames{thirdImpIdx}, 'FontSize', 11, 'Interpreter', 'none');
    ylabel(ax8, 'Tsai-Wu Criterion', 'FontSize', 11);
    title(ax8, ['3rd Most Influential: ', pNames{thirdImpIdx}], 'FontSize', 12, 'FontWeight', 'bold');
    colormap(ax8, jet);
    cb = colorbar(ax8);
    ylabel(cb, 'TW Value', 'FontSize', 10);
    grid(ax8, 'on');
end

fprintf('Analysis complete! Navigate between tabs to view different analyses.\n\n');


