% -------------------------------------------------------------------------
% 3D Cartesian B1 Mapping — Bloch-Siegert Shift Method
% Author: Valentin Jost (2025)
%
% Generates a Pulseq .seq file for 3D Cartesian B1 mapping using the
% Bloch-Siegert (BS) shift method. A 90° block pulse excites the
% magnetisation; a subsequent off-resonance Fermi pulse imparts a phase
% shift proportional to the local B1 amplitude squared:
%
%   phi_BS = Kbs * B1_local^2
%
% Two complete 3D acquisitions are written into a single .seq file with
% the Fermi pulse applied at +bs_offset and -bs_offset, respectively.
% The B1 map is recovered during reconstruction from the phase difference
% between the two acquisitions.
%
% Reference: Sacolick et al., MRM 2010. doi:10.1002/mrm.22357
%
% This code is provided for research and educational use only.
% Users must verify all safety limits (SAR, gradients, duty cycle, etc.)
% before execution on any scanner. Provided "as is" without warranty.
% -------------------------------------------------------------------------

clear; clc; close all;

%% System Configuration
% Adapt MaxGrad, MaxSlew, and RF/ADC dead times to your scanner!

gamma = 11.262e6;   % Gyromagnetic ratio for 23Na [Hz/T]

sys = mr.opts('B0', 3, 'MaxGrad', 45, 'GradUnit', 'mT/m', ...
              'MaxSlew', 200, 'SlewUnit', 'T/m/s', ...
              'rfRingdownTime', 10e-6, 'rfDeadTime', 100e-6, ...
              'adcDeadTime', 10e-6, 'gamma', gamma);

seq = mr.Sequence(sys);

%% Sequence Parameters

rf_dur = 500e-6;                        % Excitation pulse duration [s]
alpha  = 90;                            % Excitation flip angle [degrees]

fov = [240e-3, 240e-3, 240e-3];         % Field of view [m]
Nx = 40;  Ny = 40;  Nz = 40;           % Matrix size [freq, phase, partition]

Tread = 3e-3;                           % Readout flat-top duration [s]
Tpre  = 0.8e-3;                         % Pre/re-phaser duration [s]
TR    = 150e-3;                         % Repetition time [s]

%% Gradient and ADC Events

deltak = 1 ./ fov;   % k-space step size [1/m]

gx    = mr.makeTrapezoid('x', sys, 'FlatArea', Nx*deltak(1), 'FlatTime', Tread);
adc   = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
gxPre = mr.makeTrapezoid('x', sys, 'Area', -gx.area/2, 'Duration', Tpre);  % Dephaser to -kx_max
gxRe  = mr.makeTrapezoid('x', sys, 'Area', -gx.area/2, 'Duration', Tpre);  % Rewinder after readout

areaY = ((0:Ny-1) - Ny/2) * deltak(2);   % Phase-encode areas [1/m]
areaZ = ((0:Nz-1) - Nz/2) * deltak(3);   % Partition-encode areas [1/m]

%% Excitation Pulse

rf = mr.makeBlockPulse(deg2rad(alpha), 'Duration', rf_dur, 'System', sys,'use','excitation');

%% Bloch-Siegert Fermi Pulse
% A Fermi-shaped off-resonance pulse is used as the BS pulse because its
% smooth edges reduce spectral leakage. The pulse is applied at +/- bs_offset
% in the two acquisitions.

T     = 2e-3;              % Fermi pulse duration [s]
t_res = sys.rfRasterTime;  % RF raster time [s]
phi   = 0;                 % Initial phase of BS pulse [rad]

% Fermi window shape; endpoints forced to zero to avoid hard edges
rf_amp    = fermiwin(T/t_res, T);
rf_amp(1) = 0;
rf_amp(end) = 0;

rf_phase   = phi * ones(size(rf_amp));       % Constant phase across pulse
rf_complex = rf_amp .* exp(1i * rf_phase);

% Scale peak B1 amplitude: 7.41 is the empirical multiplier for a 90° pulse
% on 23Na with a 500 µs block pulse and a 4 ms Fermi envelope.
frac      = 1;
frac_peak = frac * (T / 4e-3);
B1peak    = frac_peak * 7.41;   % Peak B1 of the Fermi pulse [µT]

rf_blosi = mr.makeArbitraryRf(rf_complex, B1peak, 'System', sys,'use','other');

%% Bloch-Siegert Phase Constant (Kbs)
% Kbs encodes the relationship between local B1 amplitude and the accrued
% BS phase: phi_BS = Kbs * B1_local^2.
% Computed in Gauss units to match the conventional BS formulation.

bs_offset   = 1e3;               % Off-resonance frequency [Hz]
gamma_rad   = gamma * 2 * pi;    % Gyromagnetic ratio [rad/s/T]
gamma_gauss = gamma_rad * 1e-4;  % Gyromagnetic ratio [rad/s/G]  (1 T = 1e4 G)

% Kbs = integral(B1_norm^2 dt) * gamma^2 / (4*pi*off_res)
% Normalised to peak amplitude so Kbs is independent of absolute B1 scaling.
Kbs = sum(rf_amp.^2) * (gamma_gauss)^2 * t_res / (4 * pi * bs_offset) / (max(rf_amp))^2;

%% Delay Preparation

TRrest  = TR - mr.calcDuration(rf) - mr.calcDuration(rf_blosi) ...
             - mr.calcDuration(gx) - mr.calcDuration(gxPre)...
             -mr.calcDuration(gxRe) - 100e-6;
delayTR = mr.makeDelay(TRrest);

%% K-Space Sampling Mask
% Select undersampling pattern. For fully sampled acquisitions set pattern
% to anything other than 'block' or 'gaussian'. Custom mask functions
% (block2D_mask, gaussian2D_mask) must be available on the MATLAB path.

pattern = 'full';   % Options: 'block', 'gaussian', 'full'

if strcmp(pattern, 'block')
    len_x = 14;
    len_y = 14;
    sampling_mask = block2D_mask(Nx, Ny, len_x, len_y);
elseif strcmp(pattern, 'gaussian')
    sampling_mask = gaussian2D_mask(Nx, Ny, 5, 2, 50);
else
    sampling_mask = ones(Nx, Ny, Nz);   % Fully sampled
    len_x = Ny;
    len_y = Nz;
end

%% Sequence Loops
% The outer loop runs the full 3D acquisition twice:
%   offres_i = 1 : Fermi pulse at +bs_offset Hz
%   offres_i = 2 : Fermi pulse at -bs_offset Hz
% The B1 map is recovered from the phase difference in reconstruction.

for offres_i = 1:2

    % Alternate the frequency offset of the Bloch-Siegert pulse
    if offres_i == 1
        rf_blosi.freqOffset =  bs_offset;
    else
        rf_blosi.freqOffset = -bs_offset;
    end

    for iZ = 1:Nz

        % Partition-encode pre- and re-phasers for current Z step
        gzPre = mr.makeTrapezoid('z', 'Area',  areaZ(iZ), 'Duration', Tpre);
        gzRe  = mr.makeTrapezoid('z', 'Area', -areaZ(iZ), 'Duration', Tpre);

        for iY = 1:Ny

            % Skip lines excluded by the undersampling mask
            if sampling_mask(iY, iZ)

                % Phase-encode pre- and re-phasers for current Y step
                gyPre = mr.makeTrapezoid('y', 'Area',  areaY(iY), 'Duration', Tpre);
                gyRe  = mr.makeTrapezoid('y', 'Area', -areaY(iY), 'Duration', Tpre);

                % Excitation
                seq.addBlock(rf);

                % Short delay between excitation and BS pulse to avoid RF overlap
                seq.addBlock(mr.makeDelay(100e-6));

                % Off-resonance Bloch-Siegert pulse
                seq.addBlock(rf_blosi);

                % Simultaneous x/y/z pre-phasers
                seq.addBlock(gxPre, gyPre, gzPre);

                % Readout
                seq.addBlock(gx, adc);

                % Simultaneous x/y/z re-phasers
                seq.addBlock(gxRe, gyRe, gzRe);

                % TR padding delay
                seq.addBlock(delayTR);

            end
        end
    end
end

%% Timing Check

fprintf('Timing check...');
[ok, error_report] = seq.checkTiming;
if ok
    fprintf(' passed successfully\n');
else
    fprintf(' failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% seq.plot();
return

%% Report

rep = seq.testReport;
fprintf([rep{:}]);
return

%% Write Sequence File

seq.setDefinition('FOV',         fov);
seq.setDefinition('Name',        'B1map_BS');
seq.setDefinition('mat',         [Nx, Ny, Nz]);
seq.setDefinition('Kbs',         Kbs);
seq.setDefinition('mask',        sampling_mask);
seq.setDefinition('pattern',     pattern);
seq.setDefinition('len_fermi',   T);
seq.setDefinition('fermi_res',   t_res);
seq.setDefinition('Fermi_shape', rf_amp);
seq.setDefinition('off_freq',    bs_offset);

% Store undersampling pattern parameters
if strcmp(pattern, 'block') || strcmp(pattern, 'gaussian')
    seq.setDefinition('pat_param', [len_x, len_y]);
end

savestr = '\\B1_BSSacc.seq';
seq.write(savestr);

return

% =========================================================================
%% Helper Function
% =========================================================================

function [W, varargout] = fermiwin(L, t_length, alpha)
% fermiwin — Generate a Fermi window function
%
% Produces a smooth, bell-shaped amplitude envelope suitable for RF pulse
% shaping. The Fermi function is defined as:
%
%   F(x) = 1 / ( exp( (|x| - x0) / alpha ) + 1 )
%
% where x0 controls the flat-top width and alpha controls the edge
% steepness. Smaller alpha → sharper edges; larger alpha → smoother taper.
%
% Inputs:
%   L        : Number of samples
%   t_length : Total duration of the window [s] (default: 1)
%   alpha    : Edge steepness parameter (default: 1/33.81)
%
% Outputs:
%   W           : Fermi window amplitude [L x 1]
%   varargout{1}: Time axis [0, t_length] with L points (if t_length given)

    if nargin < 3
        alpha = 1/33.81;   % Default edge steepness
    end
    if nargin < 2
        t_length = 1;      % Default duration [s]
    end

    ep = linspace(0, 1, L) - 0.5;   % Normalised symmetric time axis [-0.5, 0.5]
    t0 = alpha * 10;                  % Flat-top half-width (scales with alpha)

    W = 1 ./ (exp((abs(ep) - t0) / alpha) + 1);   % Fermi envelope

    if nargin > 1
        varargout{1} = linspace(0, t_length, L);   % Optional time axis output
    end

end