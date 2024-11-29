% ---------------------------------------------
% Generate Multi-Layered Leg Phantom
% ---------------------------------------------

% Define dimensions and parameters
rows = 256; cols = 256; % 2D dimensions
dimX = 128; dimY = 128; dimZ = 128; % 3D dimensions
length = 1.0; % Length of leg in X direction

% Radii for different layers
outerRadiusX = 0.9; outerRadiusY = 0.7; % Skin
fatRadiusX = 0.7; fatRadiusY = 0.5; % Fat
muscleRadiusX = 0.5; muscleRadiusY = 0.35; % Muscle
boneRadiusX = 0.2; boneRadiusY = 0.2; % Bone

% Generate and display the 3D phantom
outerRadiusZ = 0.7; fatRadiusZ = 0.5; muscleRadiusZ = 0.35; boneRadiusZ = 0.2;
phantom3D = generate3DLegWithLayers(dimX, dimY, dimZ, outerRadiusY, outerRadiusZ, ...
         fatRadiusY, fatRadiusZ, muscleRadiusY, ...
         muscleRadiusZ, boneRadiusY, boneRadiusZ, length);

figure;
isosurface(phantom3D, 0.5);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Multi-Layer Leg Phantom');
axis equal;
light;
lighting phong;

% ---------------------------------------------
% X-Ray Signal Intensity Control
% ---------------------------------------------
energyLevel = 60; % Energy level in keV
I0 = 100; % Initial X-ray intensity

% Adjust muValues based on energy level
if energyLevel < 60
    muValues = [0.25, 0.18, 0.12, 0.08]; % Higher mu for lower energy
elseif energyLevel > 80
    muValues = [0.15, 0.12, 0.08, 0.05]; % Lower mu for higher energy
else
    muValues = [0.2, 0.15, 0.1, 0.05]; % Default values
end

% Simulate and display 2D X-ray projection with SI control
projection2D = generate2DProjectionWithIntensity(phantom3D, muValues, I0);

figure;
imagesc(projection2D);
colormap(gray);
axis equal tight;
title(['2D Projection with Intensity I_0 = ', num2str(I0), ' (Energy Level: ', num2str(energyLevel), ' keV)']);

% Simulate orthogonal and angled fractures
gapSize = 0.05;
simulateAndVisualizeFractures(phantom3D, muValues, gapSize, I0);

% ---------------------------------------------
% Functions
% ---------------------------------------------

function phantom3D = generate3DLegWithLayers(dimX, dimY, dimZ, outerRadiusY, outerRadiusZ, ...
                                             fatRadiusY, fatRadiusZ, muscleRadiusY, ...
                                             muscleRadiusZ, boneRadiusY, boneRadiusZ, length)
    % Create 3D phantom
    [x, y, z] = ndgrid(linspace(-1, 1, dimX), linspace(-1, 1, dimY), linspace(-1, 1, dimZ));

    % Layers
    skin = ((y / outerRadiusY).^2 + (z / outerRadiusZ).^2 <= 1) & (abs(x) <= length / 2);
    fat = ((y / fatRadiusY).^2 + (z / fatRadiusZ).^2 <= 1) & (abs(x) <= length / 2);
    muscle = ((y / muscleRadiusY).^2 + (z / muscleRadiusZ).^2 <= 1) & (abs(x) <= length / 2);
    bone = ((y / boneRadiusY).^2 + (z / boneRadiusZ).^2 <= 1) & (abs(x) <= length / 2);
    phantom3D = skin + 2 * fat + 3 * muscle + 4 * bone;
end

% Generate 2D projection from the 3D phantom with SI control
function projection2D = generate2DProjectionWithIntensity(phantom3D, muValues, I0)
    projection2D = zeros(size(phantom3D, 1), size(phantom3D, 2));
    for layer = 1:length(muValues)
        attenuation = exp(-muValues(layer) * (phantom3D == layer));
        projection2D = projection2D + sum(I0 .* attenuation, 3);
    end
end

function simulateAndVisualizeFractures(phantom3D, muValues, gapSize, I0)
    % Orthogonal fracture
    orthogonalPhantom = applyFracture(phantom3D, 90, gapSize);
    orthogonalProjection = generate2DProjectionWithIntensity(orthogonalPhantom, muValues, I0);
    figure;
    imagesc(orthogonalProjection);
    colormap(gray);
    axis equal tight;
    title('X-Ray Projection of Leg Phantom with Orthogonal Fracture');

    % Angled fracture
    angledPhantom = applyFracture(phantom3D, 45, gapSize);
    angledProjection = generate2DProjectionWithIntensity(angledPhantom, muValues, I0);
    figure;
    imagesc(angledProjection);
    colormap(gray);
    axis equal tight;
    title('X-Ray Projection of Leg Phantom with Angled Fracture');
end

% Applies fracture to the phantom
function fracturedPhantom = applyFracture(phantom3D, angle, gapSize)
    % Get phantom dimensions
    [dimX, dimY, dimZ] = size(phantom3D);

    % Create X and Z grids for fracture calculation
    [x, z] = ndgrid(1:dimX, 1:dimZ);

    % Center X and Z coordinates
    x = x - dimX / 2;
    z = z - dimZ / 2;

    % Define fracture plane
    if angle == 90
        % Orthogonal fracture
        fracturePlane = abs(x) <= gapSize / 2;
    else
        % Angled fracture
        fracturePlane = abs(x - tan(deg2rad(angle)) * z) <= gapSize / 2;
    end

    % Initialize fractured phantom
    fracturedPhantom = phantom3D;

    % Apply fracture across all Y-slices
    for y = 1:dimY
        currentSlice = squeeze(fracturedPhantom(:, y, :));
        currentSlice(fracturePlane) = 0;
        fracturedPhantom(:, y, :) = currentSlice;
    end
end
