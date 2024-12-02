% -------------------------------
% Generate or Load Data
% -------------------------------
if exist('phantom_and_projection.mat', 'file')
    % Load pre-generated data
    load('phantom_and_projection.mat', 'phantom3D', 'projection2D');
    disp('Data loaded successfully from phantom_and_projection.mat');
else
    % Generate new data if file doesn't exist
    disp('Generating new data...');

    % Holly: Generate 3D Leg Phantom
    dimX = 128; dimY = 128; dimZ = 128;
    outerRadiusY = 0.7; outerRadiusZ = 0.7;
    boneRadiusY = 0.3; boneRadiusZ = 0.3;
    length = 1.0;

    phantom3D = generate3DLegWithLayers(dimX, dimY, dimZ, outerRadiusY, outerRadiusZ, ...
                                        boneRadiusY, boneRadiusZ, length);

    % Jess: Generate 2D Projection
    energyLevel = 60; % Energy level in keV
    I0 = 100; % Initial X-ray intensity
    muValues = [0.2, 0.15, 0.1, 0.05]; % Default mu values
    projection2D = generate2DProjectionWithIntensity(phantom3D, muValues, I0);

    % Save generated data for reuse
    save('phantom_and_projection.mat', 'phantom3D', 'projection2D');
    disp('Data saved to phantom_and_projection.mat');
end

% Confirm data validity
if ~exist('phantom3D', 'var') || ~exist('projection2D', 'var')
    error('Required data (phantom3D or projection2D) is missing.');
end

% Display data dimensions
disp(['phantom3D size: ', mat2str(size(phantom3D))]);
disp(['projection2D size: ', mat2str(size(projection2D))]);

% Quick visualization
figure;
imagesc(projection2D);
colormap(gray);
axis equal tight;
title('Loaded 2D Projection Data');

% -------------------------------
% Signal Intensity and Contrast Analysis
% -------------------------------
analyze_intensity_and_contrast(projection2D);

% -------------------------------
% Image Difference Analysis
% -------------------------------
difference = analyze_image_difference(phantom3D(:, :, round(size(phantom3D, 3) / 2)), projection2D);

% -------------------------------
% Signal Intensity Profiles
% -------------------------------
plot_intensity_profile(projection2D, 'row', round(size(projection2D, 1) / 2));

% -------------------------------
% Function Definitions
% -------------------------------

% Function to generate a multi-layered 3D phantom
function phantom3D = generate3DLegWithLayers(dimX, dimY, dimZ, outerRadiusY, outerRadiusZ, ...
                                             boneRadiusY, boneRadiusZ, length)
    % Create 3D phantom
    [x, y, z] = ndgrid(linspace(-1, 1, dimX), linspace(-1, 1, dimY), linspace(-1, 1, dimZ));

    % Layers
    skin = ((y / outerRadiusY).^2 + (z / outerRadiusZ).^2 <= 1) & (abs(x) <= length / 2);
    bone = ((y / boneRadiusY).^2 + (z / boneRadiusZ).^2 <= 1) & (abs(x) <= length / 2);
    phantom3D = skin + 2 * bone; % Multi-layer encoding: 1 for skin, 2 for bone
end

% Function to generate 2D projection with intensity control
function projection2D = generate2DProjectionWithIntensity(phantom3D, muValues, I0)
    projection2D = zeros(size(phantom3D, 1), size(phantom3D, 2));
    for layer = 1:length(muValues)
        attenuation = exp(-muValues(layer) * (phantom3D == layer));
        projection2D = projection2D + sum(I0 .* attenuation, 3);
    end
end

% Analyze signal intensity and contrast
function analyze_intensity_and_contrast(projection2D)
    % Display intensity range for debugging
    minIntensity = min(projection2D(:));
    maxIntensity = max(projection2D(:));
    disp(['Min projection2D: ', num2str(minIntensity)]);
    disp(['Max projection2D: ', num2str(maxIntensity)]);

    % Dynamically adjust thresholds based on range
    skinMask = (projection2D > minIntensity) & (projection2D <= (minIntensity + 0.3 * (maxIntensity - minIntensity)));
    boneMask = (projection2D > (minIntensity + 0.3 * (maxIntensity - minIntensity)));

    % Extract intensity values
    skin = projection2D(skinMask);
    bone = projection2D(boneMask);

    % Avoid NaN if regions are empty
    if isempty(skin) || isempty(bone)
        warning('Regions for skin or bone are empty. Adjust thresholds.');
        return;
    end

    % Calculate mean signal intensity
    meanSkin = mean(skin(:));
    meanBone = mean(bone(:));

    % Calculate contrast
    contrast = abs(meanSkin - meanBone) / (meanSkin + meanBone);

    % Display results
    fprintf('Mean Signal Intensity - Skin: %.2f\n', meanSkin);
    fprintf('Mean Signal Intensity - Bone: %.2f\n', meanBone);
    fprintf('Contrast (Skin vs Bone): %.2f\n', contrast);
end



% Analyze image difference
function difference = analyze_image_difference(originalPhantom, projection2D)
    % Resize phantom to match projection dimensions if necessary
    resizedPhantom = imresize(originalPhantom, size(projection2D));

    % Calculate absolute difference
    difference = abs(resizedPhantom - projection2D);

    % Display difference image
    figure;
    imagesc(difference);
    colormap(jet);
    colorbar;
    title('Image Difference (Original vs Projection)');
end

% Plot signal intensity profile
function plot_intensity_profile(projection2D, direction, index)
    % Extract intensity profile
    if strcmp(direction, 'row')
        profile = projection2D(index, :);
    elseif strcmp(direction, 'column')
        profile = projection2D(:, index);
    else
        error('Invalid direction. Use "row" or "column".');
    end

    % Plot the profile
    figure;
    plot(profile);
    xlabel('Position');
    ylabel('Signal Intensity');
    title('Signal Intensity Profile', 'Interpreter', 'none');
end
