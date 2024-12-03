% -------------------------------
% Generate or Load Data
% -------------------------------
if exist('phantom_and_projection.mat', 'file')
    load('phantom_and_projection.mat', 'phantom3D', 'projection2D');
    if ~exist('phantom3D', 'var') || ~exist('projection2D', 'var')
        error('Required data (phantom3D or projection2D) is missing in phantom_and_projection.mat.');
    end
    disp('Data loaded successfully from phantom_and_projection.mat');
else
    error('phantom_and_projection.mat file is missing. Run parts 1 and 2 first.');
end

% -------------------------------
% Simulate Fractures
% -------------------------------
disp('Simulating fractures...');

% Generate fractured phantoms with adjusted gapSize
gapSize = 1; % Gap size for horizontal split
phantom3DOrthogonal = applyFracture(phantom3D, 0, gapSize); % Horizontal split
phantom3DAngled = applyFracture(phantom3D, 45, gapSize); % Angled split

% Generate 2D projections
I0 = 1; % Initial X-ray intensity set to 1
muValues = [1, 0.6, 0.05, 0.05]; % Default mu values

projection2DOrthogonal = generate2DProjectionWithIntensity(phantom3DOrthogonal, muValues, I0);
projection2DAngled = generate2DProjectionWithIntensity(phantom3DAngled, muValues, I0);

gamma = 0.5; % Adjust this value for better brightness
projection2DOrthogonal = projection2DOrthogonal.^gamma;
projection2DAngled = projection2DAngled.^gamma;
% Save fractured data for future use
save('phantom_and_projection.mat', 'phantom3D', 'projection2D', ...
    'phantom3DOrthogonal', 'projection2DOrthogonal', ...
    'phantom3DAngled', 'projection2DAngled', '-append');
disp('Fractures and projections saved.');

% -------------------------------
% Visualizations
% -------------------------------
disp('Visualizing projections...');

% Orthogonal fracture projection
figure;
imagesc(projection2DOrthogonal);
colormap(gray);
axis equal tight;
title('Projection with Orthogonal Fracture');

% Angled fracture projection
figure;
imagesc(projection2DAngled);
colormap(gray);
axis equal tight;
title('Projection with Angled Fracture');

% -------------------------------
% Analyze Signal Intensity and Contrast
% -------------------------------
disp('Analyzing intensity and contrast...');
analyze_intensity_and_contrast(projection2DOrthogonal);
analyze_intensity_and_contrast(projection2DAngled);

% -------------------------------
% Analyze Contrast Across Fractures
% -------------------------------
disp('Analyzing contrast across fractures...');

% After gamma correction, check intensity range
minIntensityOrthogonal = min(projection2DOrthogonal(:));
maxIntensityOrthogonal = max(projection2DOrthogonal(:));
disp(['Orthogonal Projection - Min Intensity: ', num2str(minIntensityOrthogonal)]);
disp(['Orthogonal Projection - Max Intensity: ', num2str(maxIntensityOrthogonal)]);

minIntensityAngled = min(projection2DAngled(:));
maxIntensityAngled = max(projection2DAngled(:));
disp(['Angled Projection - Min Intensity: ', num2str(minIntensityAngled)]);
disp(['Angled Projection - Max Intensity: ', num2str(maxIntensityAngled)]);

% Calculate splitThreshold dynamically for each projection
splitThresholdOrthogonal = minIntensityOrthogonal + 0.5 * (maxIntensityOrthogonal - minIntensityOrthogonal);
splitThresholdAngled = minIntensityAngled + 0.5 * (maxIntensityAngled - minIntensityAngled);

% Display the calculated thresholds
disp(['Calculated splitThresholdOrthogonal: ', num2str(splitThresholdOrthogonal)]);
disp(['Calculated splitThresholdAngled: ', num2str(splitThresholdAngled)]);

% Create split masks based on the calculated thresholds
splitMaskOrthogonal = projection2DOrthogonal < splitThresholdOrthogonal;
splitMaskAngled = projection2DAngled < splitThresholdAngled;

% Analyze contrast
analyze_split_contrast(projection2DOrthogonal, splitMaskOrthogonal);
analyze_split_contrast(projection2DAngled, splitMaskAngled);


% -------------------------------
% Signal Intensity Profiles
% -------------------------------
disp('Plotting signal intensity profiles...');

% Ensure only valid profiles are plotted
plot_intensity_profile(projection2DOrthogonal, 'row', round(size(projection2DOrthogonal, 1) / 2), 'Orthogonal Fracture');
plot_intensity_profile(projection2DAngled, 'row', round(size(projection2DAngled, 1) / 2), 'Angled Fracture');

% -------------------------------
% Function Definitions
% -------------------------------

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

    % Apply fracture only to the bone region (assuming bone is layer 2)
    boneRegion = (phantom3D == 3); % Assuming '2' is the bone label in the phantom
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
    projection2D = zeros(size(phantom3D, 1), size(phantom3D, 2));
    for layer = 1:length(muValues)
        attenuation = exp(-muValues(layer) * (phantom3D == layer));
        projection2D = projection2D + sum(I0 .* attenuation, 3);
    end
end

% Analyze signal intensity and contrast
function analyze_intensity_and_contrast(projection2D)
    minIntensity = min(projection2D(:));
    maxIntensity = max(projection2D(:));
    disp(['Min projection2D: ', num2str(minIntensity)]);
    disp(['Max projection2D: ', num2str(maxIntensity)]);

    % Define masks dynamically based on intensity thresholds
    skinMask = (projection2D > minIntensity) & (projection2D <= (minIntensity + 0.3 * (maxIntensity - minIntensity)));
    boneMask = (projection2D > (minIntensity + 0.3 * (maxIntensity - minIntensity)));

    skin = projection2D(skinMask);
    bone = projection2D(boneMask);

    if isempty(skin) || isempty(bone)
        warning('Regions for skin or bone are empty. Adjust thresholds.');
        return;
    end

    meanSkin = mean(skin(:));
    meanBone = mean(bone(:));
    contrast = abs(meanSkin - meanBone) / (meanSkin + meanBone);

    fprintf('Mean Signal Intensity - Skin: %.2f\n', meanSkin);
    fprintf('Mean Signal Intensity - Bone: %.2f\n', meanBone);
    fprintf('Contrast (Skin vs Bone): %.2f\n', contrast);
end

% Analyze contrast across fractures
function analyze_split_contrast(projection2D, splitMask)
    splitRegion = projection2D(splitMask);
    nonSplitRegion = projection2D(~splitMask);

    if isempty(splitRegion) || isempty(nonSplitRegion)
        warning('No valid regions found for split contrast analysis.');
        return;
    end

    meanSplit = mean(splitRegion(:));
    meanNonSplit = mean(nonSplitRegion(:));
    contrast = abs(meanSplit - meanNonSplit) / (meanSplit + meanNonSplit);

    fprintf('Contrast across split: %.2f\n', contrast);
end

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