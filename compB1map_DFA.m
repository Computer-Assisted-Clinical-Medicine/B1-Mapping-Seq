function B1map_smoothed = compB1DFAmap(filepath, abs_map)
% -------------------------------------------------------------------------
% compB1DFAmap — B1 Map Reconstruction via Double Flip Angle Method
% Author: Valentin Jost (2025)
%
% Reconstructs a 3D B1 efficiency map from two spoiled GRE acquisitions
% at different flip angles (alpha1, alpha2 = 2*alpha1). The local flip
% angle is estimated from the signal ratio using:
%
%   r = acos( S(2*alpha) / (2 * S(alpha)) )
%   B1_rel = r / alpha_nominal
%
% Optionally returns an absolute B1 map in Tesla.
%
% Inputs:
%   filepath  : Path to the .seq/.dat files (with or without extension)
%   abs_map   : If true, returns absolute B1 in Tesla; otherwise relative
%               B1 efficiency (default: false)
%
% Output:
%   B1map_smoothed : 3D B1 map, Gaussian smoothed (sigma = 1 voxel)
%                    Relative [a.u.] or absolute [T] depending on abs_map
%
% Dependencies:
%   mapVBVD, ifft3c_new
%
% This code is provided for research and educational use only.
% Provided "as is" without warranty of any kind.
% -------------------------------------------------------------------------

%% Input Handling

if nargin < 2
    abs_map = 0;
end

%% Read Sequence Parameters

% Strip file extension and resolve paths for .seq and .dat files
filepath = filepath(1:end-4);
seqpath  = strcat(filepath, '.seq');
datapath = strcat(filepath, '.dat');

seq = mr.Sequence();
seq.read(seqpath, 'detectRFuse');

matsize = seq.getDefinition('mat');       % Acquisition matrix [Nx, Ny, Nz]
pattern = seq.getDefinition('pattern');   % Sampled lines [len_y, len_z]
alphas  = seq.getDefinition('alphas');    % Flip angles [alpha1, alpha2]

% Fall back to full matrix if no undersampling pattern was stored
if numel(pattern) == 0
    pattern = matsize;
end

% Identify the smaller flip angle as the nominal reference angle
if alphas(1) < alphas(2)
    alpha = deg2rad(alphas(1));
else
    alpha = deg2rad(alphas(2));
end

%% Load Raw K-Space Data

% mapVBVD may return a cell array for multi-RAID files; always use the last entry
twix_obj = mapVBVD(datapath);
if iscell(twix_obj)
    data_unsorted = twix_obj{end}.image.unsorted();
else
    data_unsorted = twix_obj.image.unsorted();
end

clear twix_obj;

%% Sort ADC Stream into Structured Array
% Data layout: [kx, ky_sampled, kz_sampled, flip_angle, average]
% The two flip-angle volumes are acquired sequentially in the .seq file.

data_size = size(reshape(data_unsorted, matsize(1), pattern(1), pattern(2), 2, []));
if length(data_size) < 5
    data_size = [data_size, 1];   % Ensure 5-D even without averages
end

data_unsorted = squeeze(data_unsorted);
data_sorted   = zeros(data_size);

adc_i = 1;
for avg_i = 1:data_size(5)
    for alpha_i = 1:2
        for iZ = 1:pattern(2)
            for iY = 1:pattern(1)
                data_sorted(:, iY, iZ, alpha_i, avg_i) = data_unsorted(:, adc_i);
                adc_i = adc_i + 1;
            end
        end
    end
end

% Average over repeated acquisitions if present
if data_size(5) > 1
    data_sorted = squeeze(mean(data_sorted, 5));
else
    data_sorted = squeeze(data_sorted);
end

%% Reconstruct

% 3D iFFT reconstruction of each flip-angle volume
S1 = ifft3c_new(squeeze(data_sorted(:, :, :, 1)));   % Image at smaller flip angle
S2 = ifft3c_new(squeeze(data_sorted(:, :, :, 2)));   % Image at larger flip angle

% Quick preview of central slice for both flip angles
% N = round(matsize(3) / 2);
% figure;
% subplot(1, 2, 1); imagesc(abs(S1(:, :, N))); title(sprintf('S1 (\\alpha = %d°)', alphas(1))); axis image; colorbar;
% subplot(1, 2, 2); imagesc(abs(S2(:, :, N))); title(sprintf('S2 (\\alpha = %d°)', alphas(2))); axis image; colorbar;

%% B1 Map Computation
% DFA formula: r = acos( S2 / (2*S1) )
% eps prevents division by zero in noise-only voxels.

r      = acos(S2 ./ (2 * S1 + eps));
B1_rel = r ./ alpha;   % Relative B1 efficiency [a.u.]; 1.0 = nominal flip angle

if abs_map
    % Convert relative B1 to absolute B1 amplitude in Tesla:
    % B1_nominal = alpha / (gamma * tau), where tau is the pulse duration
    tau        = 0.5e-3;        % Block pulse duration [s]
    gamma_na   = 11.262e6;      % Gyromagnetic ratio for 23Na [Hz/T]
    B1_nominal = alpha / (gamma_na * tau);
    B1_abs     = B1_rel .* B1_nominal;
    B1map_smoothed = imgaussfilt(abs(B1_abs), 1);
else
    B1map_smoothed = imgaussfilt(abs(B1_rel), 1);
end

%% Display Three Orthogonal Central Slices

% figure;
% subplot(1, 3, 1); imagesc(abs(squeeze(B1map_smoothed(:, :, N)))); title('Axial');    axis image; colorbar;
% subplot(1, 3, 2); imagesc(abs(squeeze(B1map_smoothed(:, N, :)))); title('Coronal');  axis image; colorbar;
% subplot(1, 3, 3); imagesc(abs(squeeze(B1map_smoothed(N, :, :)))); title('Sagittal'); axis image; colorbar;

end