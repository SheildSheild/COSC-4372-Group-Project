% Clearing workspace and close figures to avoid overloading session
close all;
clc;

% -------------------------------
% Generate or Load Data
% -------------------------------
% We'll generate the phantom3D with the new layers
disp('Generating 3D phantom leg with skin, fat, muscle, and bone layers...');

% Define dimensions
dimX = 128; dimY = 128; dimZ = 128;
outerRadiusY = 0.6; outerRadiusZ = 0.6; % Radii for the outer leg (skin)
boneRadiusY = 0.2; boneRadiusZ = 0.2;   % Radii for the bone
length = 1.0; % Length of the leg

% Generate phantom3D with layers
phantom3D = generate3Dleg(dimX, dimY, dimZ, outerRadiusY, outerRadiusZ, boneRadiusY, boneRadiusZ, length);

% Save the phantom
save('phantom_and_projection.mat', 'phantom3D');

% -------------------------------
% Simulate Fractures
% -------------------------------
disp('Simulating fractures...');

% Generate fractured phantoms with adjusted gapSize
gapSize = 1; % Gap size for fractures
phantom3DNormal = phantom3D; % No fracture
phantom3DOrthogonal = applyFracture(phantom3D, 0, gapSize); % Orthogonal fracture
phantom3DAngled = applyFracture(phantom3D, 45, gapSize); % Angled fracture

% Generate 2D projections
I0 = 1; % Initial X-ray intensity
muValues = [0.04, 0.03, 0.015, 0.005]; % [Skin, Fat, Muscle, Bone]

projection2DNormal = generate2DProjectionWithIntensity(phantom3DNormal, muValues, I0);
projection2DOrthogonal = generate2DProjectionWithIntensity(phantom3DOrthogonal, muValues, I0);
projection2DAngled = generate2DProjectionWithIntensity(phantom3DAngled, muValues, I0);

% Apply gamma correction for better visualization
gamma = 0.01; % Adjust as needed
projection2DNormal_gamma = projection2DNormal.^gamma;
projection2DOrthogonal_gamma = projection2DOrthogonal.^gamma;
projection2DAngled_gamma = projection2DAngled.^gamma;

% Save fractured data for future use
save('phantom_and_projection.mat', 'phantom3D', 'projection2DNormal', 'projection2DOrthogonal', 'projection2DAngled', '-append');
disp('Fractures and projections saved.');

% -------------------------------
% Visualizations
% -------------------------------
disp('Visualizing projections...');

% No fracture projection
figure;
imagesc(projection2DNormal_gamma);
colormap(gray);
axis equal tight;
title('Projection with No Fracture');

% Orthogonal fracture projection
figure;
imagesc(projection2DOrthogonal_gamma);
colormap(gray);
axis equal tight;
title('Projection with Orthogonal Fracture');

% Angled fracture projection
figure;
imagesc(projection2DAngled_gamma);
colormap(gray);
axis equal tight;
title('Projection with Angled Fracture');

% -------------------------------
% Analyze Signal Intensity and Contrast
% -------------------------------
disp('Analyzing intensity and contrast...');

% Analyze using the function (adjusted for new layers)
analyze_intensity_and_contrast(projection2DNormal, 'No Fracture');
analyze_intensity_and_contrast(projection2DOrthogonal, 'Orthogonal Fracture');
analyze_intensity_and_contrast(projection2DAngled, 'Angled Fracture');

% -------------------------------
% Signal Intensity Profiles
% -------------------------------
disp('Plotting signal intensity profiles...');

% Ensure only valid profiles are plotted
index = round(size(projection2DNormal, 1) / 2);
plot_intensity_profile(projection2DNormal, 'row', index, 'No Fracture');
plot_intensity_profile(projection2DOrthogonal, 'row', index, 'Orthogonal Fracture');
plot_intensity_profile(projection2DAngled, 'row', index, 'Angled Fracture');

% -------------------------------
% Function Definitions
% -------------------------------

% Generate 3D leg phantom with layers
function phantom3D = generate3Dleg(dimX, dimY, dimZ, outerRadiusY, outerRadiusZ, boneRadiusY, boneRadiusZ, length)
    % 3D grid of coordinates
    [x, y, z] = ndgrid(linspace(-1, 1, dimX), linspace(-1, 1, dimY), linspace(-1, 1, dimZ));

    % Thigh (outer cylinder - skin layer)
    leg = ((y / outerRadiusY).^2 + (z / outerRadiusZ).^2 <= 1) & (abs(x) <= length / 2);

    % Fat layer
    fatRadiusY = outerRadiusY * 0.8;
    fatRadiusZ = outerRadiusZ * 0.8;
    fat = ((y / fatRadiusY).^2 + (z / fatRadiusZ).^2 <= 1) & (abs(x) <= length / 2);

    % Muscle layer
    muscleRadiusY = fatRadiusY * 0.8;
    muscleRadiusZ = fatRadiusZ * 0.8;
    muscle = ((y / muscleRadiusY).^2 + (z / muscleRadiusZ).^2 <= 1) & (abs(x) <= length / 2);

    % Bone (inner cylinder)
    bone = ((y / boneRadiusY).^2 + (z / boneRadiusZ).^2 <= 1) & (abs(x) <= length / 2);

    % Initialize phantom
    phantom3D = zeros(dimX, dimY, dimZ);

    % Assign labels (from outer to inner layers)
    phantom3D(leg) = 1; % Skin
    phantom3D(fat) = 2; % Fat
    phantom3D(muscle) = 3; % Muscle
    phantom3D(bone) = 4; % Bone
end

% Apply fracture to a phantom
function fracturedPhantom = applyFracture(phantom3D, angle, gapSize)
    [dimX, dimY, dimZ] = size(phantom3D);

    % Create a grid for fracture simulation in the x-y plane
    [xGrid, yGrid] = ndgrid(1:dimX, 1:dimY);
    x = xGrid - (dimX + 1) / 2;
    y = yGrid - (dimY + 1) / 2;

    % Define fracture plane in the x-y plane
    angle_rad = deg2rad(angle);
    fracturePlane = abs(x * cos(angle_rad) + y * sin(angle_rad)) <= gapSize / 2;

    % Apply fracture only to the bone region (assuming bone is layer 4)
    boneRegion = (phantom3D == 4);
    fracturedPhantom = phantom3D;

    for z = 1:dimZ
        slice = phantom3D(:, :, z); % Get the z-th slice
        boneSlice = boneRegion(:, :, z); % Identify bone in the slice
        fractureMask = fracturePlane & boneSlice; % Apply fracture only to bone
        slice(fractureMask) = 0; % Remove the fractured region
        fracturedPhantom(:, :, z) = slice; % Update the fractured phantom
    end
end

% Generate 2D projection with intensity control
function projection2D = generate2DProjectionWithIntensity(phantom3D, muValues, I0)
    mu3D = zeros(size(phantom3D), 'double'); % Use double precision for accuracy
    for layer = 1:length(muValues)
        mu3D(phantom3D == layer) = muValues(layer);
    end
    % Sum mu over z-direction (depth)
    total_mu = sum(mu3D, 3);
    % Calculate projection using the Beer-Lambert law
    projection2D = I0 * exp(-total_mu);
end

% Analyze signal intensity and contrast
function analyze_intensity_and_contrast(projection2D, projectionType)
    disp(['Analyzing intensity and contrast for ', projectionType, '...']);
    
    minIntensity = min(projection2D(:));
    maxIntensity = max(projection2D(:));
    disp(['Min projection2D: ', num2str(minIntensity)]);
    disp(['Max projection2D: ', num2str(maxIntensity)]);
    
    % Define thresholds for each layer based on intensity distribution
    thresholds = linspace(minIntensity, maxIntensity, 5); % 4 intervals for 4 layers
    
    % Create masks
    skinMask = projection2D > thresholds(4);
    fatMask = (projection2D > thresholds(3)) & (projection2D <= thresholds(4));
    muscleMask = (projection2D > thresholds(2)) & (projection2D <= thresholds(3));
    boneMask = projection2D <= thresholds(2);
    
    % Extract intensities
    skinIntensity = projection2D(skinMask);
    fatIntensity = projection2D(fatMask);
    muscleIntensity = projection2D(muscleMask);
    boneIntensity = projection2D(boneMask);
    
    % Check for empty masks
    if isempty(skinIntensity)
        warning('Skin region is empty.');
    end
    if isempty(fatIntensity)
        warning('Fat region is empty.');
    end
    if isempty(muscleIntensity)
        warning('Muscle region is empty.');
    end
    if isempty(boneIntensity)
        warning('Bone region is empty.');
    end
    
    % Calculate mean intensities
    meanSkin = mean(skinIntensity);
    meanFat = mean(fatIntensity);
    meanMuscle = mean(muscleIntensity);
    meanBone = mean(boneIntensity);
    
    % Calculate contrasts
    contrastSkinFat = abs(meanSkin - meanFat) / (meanSkin + meanFat);
    contrastFatMuscle = abs(meanFat - meanMuscle) / (meanFat + meanMuscle);
    contrastMuscleBone = abs(meanMuscle - meanBone) / (meanMuscle + meanBone);
    
    % Display results with projection type
    fprintf('--- %s ---\n', projectionType);
    fprintf('Mean Signal Intensity - Skin: %.4f\n', meanSkin);
    fprintf('Mean Signal Intensity - Fat: %.4f\n', meanFat);
    fprintf('Mean Signal Intensity - Muscle: %.4f\n', meanMuscle);
    fprintf('Mean Signal Intensity - Bone: %.4f\n', meanBone);
    fprintf('Contrast (Skin vs Fat): %.4f\n', contrastSkinFat);
    fprintf('Contrast (Fat vs Muscle): %.4f\n', contrastFatMuscle);
    fprintf('Contrast (Muscle vs Bone): %.4f\n', contrastMuscleBone);
end


% Plot intensity profile
function plot_intensity_profile(projection2D, direction, index, projectionType)
    % Validate input
    if isempty(projection2D) || all(projection2D(:) == 0)
        warning('Invalid projection data. Skipping intensity profile plot.');
        return;
    end

    % Extract profile
    if strcmp(direction, 'row')
        profile = projection2D(index, :);
        disp(['Plotting row intensity profile at index: ', num2str(index), ' for ', projectionType]);
    elseif strcmp(direction, 'column')
        profile = projection2D(:, index);
        disp(['Plotting column intensity profile at index: ', num2str(index), ' for ', projectionType]);
    else
        error('Invalid direction. Use "row" or "column".');
    end

    % Debug: Check profile data
    if all(profile == 0)
        warning(['Profile data is all zeros for ', projectionType, '. Skipping plot.']);
        return;
    end

    % Plot the profile
    figure;
    plot(profile);
    xlabel('Position');
    ylabel('Signal Intensity');
    title(['Signal Intensity Profile - ', projectionType], 'Interpreter', 'none');
end
