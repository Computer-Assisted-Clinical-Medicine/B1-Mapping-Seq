function mask = gaussian2D_mask(nx, ny, sigma_x, sigma_y, acc)
%GAUSSIAN2D_MASK Create a random 2D Gaussian sampling mask with independent sigmas
%   mask = gaussian2D_mask(nx, ny, sigma_x, sigma_y, acc)
%
%   nx, ny    - matrix size
%   sigma_x   - std dev of Gaussian in x direction (rows)
%   sigma_y   - std dev of Gaussian in y direction (columns)
%   acc       - desired acceleration factor (e.g., 4 means ~25% samples kept)
%
% Example:
%   mask = gaussian2D_mask(128, 128, 30, 10, 4);
%   imagesc(mask); axis image; colormap gray;

% Coordinate grids centered at zero
[x, y] = meshgrid( (-floor(ny/2)):(ceil(ny/2)-1), ...
                   (-floor(nx/2)):(ceil(nx/2)-1) );

% 2D Gaussian with independent sigmas
density = exp(-(x.^2 / (2*sigma_y^2) + y.^2 / (2*sigma_x^2)));
density = density / max(density(:)); % normalize to [0,1]

% Determine threshold for target samples based on acceleration
target_samples = round(nx * ny / acc);
flat_density = density(:);
[~, idx] = sort(flat_density, 'descend');
thresh = flat_density(idx(target_samples));

% Generate sampling mask
mask = rand(nx, ny) < density * (1/thresh);

% Always include the DC point (center)
mask(round(nx/2), round(ny/2)) = 1;
end
