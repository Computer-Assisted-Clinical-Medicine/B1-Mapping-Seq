function B1Map = compB1BSLesch(filepath, output_size)
% -------------------------------------------------------------------------
% compB1BSLesch — B1 Map Reconstruction via Bloch-Siegert Shift Method
% Author: Valentin Jost (2025)
%
% Reconstructs a 3D relative B1 efficiency map from two 
% acquisitions with the Bloch-Siegert (BS) pulse applied at +/- bs_offset.
% Reconstruction uses the TGV-regularised BS framework by Lesch et al.
%
% The BS phase shift accrued per voxel is:
%
%   phi_BS = Kbs * B1_local^2
%
% B1_local is recovered by BlochSiegReco() and normalised to the nominal
% B1 amplitude to yield a relative B1 efficiency map.
%
% Inputs:
%   filepath    : Path to the .seq/.dat files (with or without extension)
%   output_size : (Optional) Scalar. If provided, the B1 map is isotropically
%                 resampled to [output_size x output_size x output_size]
%                 using Lanczos-3 interpolation.
%
% Output:
%   B1Map       : 3D relative B1 efficiency map [a.u.]; 1.0 = nominal flip angle
%
% Dependencies:
%   mapVBVD               — raw data loading (Siemens TWIX format)
%   ifft3c_new            — 3D centred iFFT
%   BlochSiegReco         — TGV-regularised BS reconstruction
%                           (Lesch et al., MRM 2019, doi:10.1002/mrm.27434)
%                           https://github.com/IMTtugraz/BSReconFramework
%   BSS_TGV_vals.mat      — lookup table for TGV regularisation parameters
%                           (lambda, mu) indexed by undersampling pattern size;
%                           provided in this repository
%
% Reference:
%   Lesch, A. et al. (2019). Accelerated 3D Bloch-Siegert B1+ mapping
%   using variational spatial regularization. Magnetic Resonance in Medicine,
%   81(2), 839–851. https://doi.org/10.1002/mrm.27434
%
% This code is provided for research and educational use only.
% Provided "as is" without warranty of any kind.
% -------------------------------------------------------------------------

%% Read Sequence Parameters

% Strip file extension; resolve .seq and .dat paths separately
% (the .dat file sits in the same folder as the .seq file)
filepath   = filepath(1:end-4);
path_split = split(filepath, '\');
seqname    = path_split(end);
seqpath    = strcat(filepath, '.seq');
folderpath = join(path_split(1:end-1), '\');
datapath   = strcat(folderpath{1,1}, '\', seqname{1,1}, '.dat');

seq = mr.Sequence();
seq.read(seqpath);

% Recover all acquisition and BS pulse parameters stored at write time
fov       = seq.getDefinition('FOV');
matsize   = seq.getDefinition('mat');        % [Nx, Ny, Nz]
Kbs       = seq.getDefinition('Kbs');        % BS phase constant [rad/G²]
T         = seq.getDefinition('len_fermi');  % Fermi pulse duration [s]
t_res     = seq.getDefinition('fermi_res'); % RF raster time [s]
rf_amp    = seq.getDefinition('Fermi_shape'); % Normalised Fermi envelope
bs_offset = seq.getDefinition('off_freq');   % Off-resonance frequency [Hz]
pat_param = seq.getDefinition('pat_param');  % Undersampling block size [len_y, len_z]
mask_us   = seq.getDefinition('mask');       % Serialised undersampling mask

% Restore 2D undersampling mask from the serialised vector
mask_us = reshape(mask_us, matsize(2), matsize(3));

% Gyromagnetic ratio for 23Na
gamma         = 11.262e6;               % [Hz/T]
gamma_rad     = gamma * 2 * pi;         % [rad/s/T]
gamma_rad_gau = gamma_rad * 1e-4;       % [rad/s/G]

%% Load Raw K-Space Data

% mapVBVD may return a cell array for multi-RAID files; always use the last entry
twix_obj = mapVBVD(datapath);
if iscell(twix_obj)
    data_unsorted = twix_obj{end}.image.unsorted();
else
    data_unsorted = twix_obj.image.unsorted();
end

data_unsorted = squeeze(data_unsorted);
clear twix_obj;

%% Sort ADC Stream into Structured Array
% Data layout: [kx, ky_sampled, kz_sampled, off-resonance sign, average]
% Only lines present in mask_us contain ADC data; unsampled lines remain zero.

% Infer number of averages from total ADC count and sampled k-space lines
avg = size(data_unsorted, 2) / (sum(mask_us, 'all') * 2);
disp(['avgs: ' num2str(avg)]);

data_size   = [matsize(1), matsize(2), matsize(3), 2, avg];
disp(data_size);
data_sorted = zeros(data_size);

adc_i = 1;
for avg_i = 1:data_size(5)
    for offres_i = 1:2
        for iZ = 1:matsize(3)
            for iY = 1:matsize(2)
                % Only fill k-space lines that were acquired
                if mask_us(iY, iZ)
                    data_sorted(:, iY, iZ, offres_i, avg_i) = data_unsorted(:, adc_i);
                    adc_i = adc_i + 1;
                end
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

%% Reconstruct Preview Images

% Separate +/- off-resonance k-space volumes
kspace_plus  = data_sorted(:, :, :, 1);
kspace_minus = data_sorted(:, :, :, 2);

% Binary mask: acquired lines are non-zero
mask = data_sorted(:, :, :, :) ~= 0;

% 3D iFFT reconstruction of both BS acquisitions
img_plus  = ifft3c_new(kspace_plus);
img_minus = ifft3c_new(kspace_minus);

% Magnitude preview — central slice
figure;
subplot(1, 2, 1); imagesc(abs(img_plus(:, :, 20)));  title('+\Deltaf magnitude'); axis image; colorbar;
subplot(1, 2, 2); imagesc(abs(img_minus(:, :, 20))); title('-\Deltaf magnitude'); axis image; colorbar;

% Phase preview — the phase difference between the two encodes the B1 map
figure;
subplot(1, 2, 1); imagesc(angle(img_plus(:, :, 20)));  title('+\Deltaf phase'); axis image; colorbar;
subplot(1, 2, 2); imagesc(angle(img_minus(:, :, 20))); title('-\Deltaf phase'); axis image; colorbar;

%% Set Up Lesch TGV Toolbox Options
% All fields follow the BSReconFramework convention.
% See https://github.com/IMTtugraz/BSReconFramework for full documentation.

tgv_opt.Kbs          = Kbs;           % BS phase constant [rad/G²]
tgv_opt.pulse        = rf_amp;        % Normalised Fermi pulse shape
tgv_opt.Tpulse       = T;             % BS pulse duration [s]
tgv_opt.deltaf       = bs_offset;     % Off-resonance frequency [Hz]
tgv_opt.deltaOmega   = 2*pi * tgv_opt.deltaf;  % Angular off-resonance [rad/s]
tgv_opt.gamma        = gamma_rad;     % Gyromagnetic ratio [rad/s/T]
tgv_opt.deltat       = t_res;         % RF raster time [s]
tgv_opt.impl         = 'CPU';         % 'CPU' or 'GPU'
tgv_opt.pattern      = mask;          % Undersampling pattern (logical)
tgv_opt.NC           = 1;             % Number of receive channels
tgv_opt.coilSens     = 'singleCoil';  % Coil sensitivity model
tgv_opt.dimY         = matsize(2);    % Phase-encoding lines
tgv_opt.dimX         = matsize(1);    % Frequency-encoding samples
tgv_opt.dimSlice     = matsize(3);    % Partition-encoding steps
tgv_opt.y            = data_sorted;   % Undersampled k-space data
tgv_opt.PhaseEncDir  = 'Lin';         % Phase-encoding direction
tgv_opt.dx           = 1;             % Resolution in x [voxel units]
tgv_opt.dy           = 1;             % Resolution in y [voxel units]
tgv_opt.dz           = 1;             % Resolution in z [voxel units]
tgv_opt.traj         = 'cart';        % Trajectory type (Cartesian only)
tgv_opt.maxitTGV     = 200;           % Max TGV iterations
tgv_opt.maxitH1      = 200;           % Max H1 iterations

%% Load TGV Regularisation Parameters
% lambda and mu are pattern-size-dependent and stored in BSS_TGV_vals.mat.
% The lookup table is indexed by the block undersampling dimensions
% (pat_param(1) x pat_param(2)) and must be provided with this repository.

load('\\BSS_TGV_vals.mat');  % Loads: BSS_TGV_vals struct array

for xi = 2:length(BSS_TGV_vals)
    for yi = 2:length(BSS_TGV_vals)
        if (BSS_TGV_vals(xi).xpattern == pat_param(1) && ...
            BSS_TGV_vals(yi).ypattern == pat_param(2) && xi == yi)
            tgv_opt.lambda = BSS_TGV_vals(xi).lambda;
            tgv_opt.mu     = BSS_TGV_vals(xi).mu;
        end
    end
end

%% Reconstruction
% BlochSiegReco returns the absolute B1+ amplitude [T or G].
% Normalise to the nominal B1 amplitude to obtain relative B1 efficiency.

B1plus    = BlochSiegReco(tgv_opt);
B1nominal = (pi/2) / (gamma_rad * 500e-6);   % Nominal B1 for a 90° block pulse of 500 µs [T]
B1Map     = B1plus ./ B1nominal;              % Relative B1 efficiency [a.u.]

%% Display Three Orthogonal Central Slices

% figure;
% subplot(1, 3, 1); imagesc(squeeze(abs(B1Map(:, :, 20)))); title('Axial');    axis image; colorbar;
% subplot(1, 3, 2); imagesc(squeeze(abs(B1Map(:, 20, :)))); title('Coronal');  axis image; colorbar;
% subplot(1, 3, 3); imagesc(squeeze(abs(B1Map(20, :, :)))); title('Sagittal'); axis image; colorbar;

%% Optional Resampling
% If an output size is requested, resample the B1 map isotropically using
% Lanczos-3 interpolation (antialiasing disabled to preserve B1 scaling).

if nargin > 1
    B1Map = imresize3(B1Map, [output_size, output_size, output_size], ...
                      'lanczos3', 'Antialiasing', false);
end

disp('B1comp DONE.');

end