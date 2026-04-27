% -------------------------------------------------------------------------
% 3D Cartesian B1 Mapping — Double Flip Angle Method
% Author: Valentin Jost (2025)
%
% Generates a Pulseq .seq file for 3D Cartesian B1 mapping using the
% Double Flip Angle (DFA) method. Two complete 3D acquisitions
% are written into a single .seq file: the first at flip angle alpha1, the
% second at alpha2 = 2*alpha1. The local B1 efficiency is estimated during
% reconstruction from the signal ratio S(alpha2) / S(alpha1).
%
% Supports three k-space sampling patterns:
%   'full'     — fully sampled (default)
%   'block'    — rectangular undersampling mask
%   'gaussian' — variable-density Gaussian undersampling mask
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

rf_dur = 500e-6;                        % RF pulse duration [s]
alpha1 = 60;                            % First flip angle [degrees]
alpha2 = 120;                           % Second flip angle [degrees] — typically 2*alpha1

fov = [240e-3, 240e-3, 240e-3];         % Field of view [m]
Nx = 40;  Ny = 40;  Nz = 40;           % Matrix size [freq, phase, partition]

Tread = 3e-3;                           % Readout duration [s]
Tpre  = 0.8e-3;                         % Pre/re-phaser duration [s]
TR    = 250e-3;                         % Repetition time [s]

%% Gradient and ADC Events

deltak = 1 ./ fov;   % k-space step size [1/m]

gx    = mr.makeTrapezoid('x', sys, 'FlatArea', Nx*deltak(1), 'FlatTime', Tread);
adc   = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
gxPre = mr.makeTrapezoid('x', sys, 'Area', -gx.area/2, 'Duration', Tpre);  % Dephaser to -kx_max
gxRe  = mr.makeTrapezoid('x', sys, 'Area', -gx.area/2, 'Duration', Tpre);  % Rewinder after readout

areaY = ((0:Ny-1) - Ny/2) * deltak(2);   % Phase-encode areas [1/m]
areaZ = ((0:Nz-1) - Nz/2) * deltak(3);   % Partition-encode areas [1/m]

%% Delay Preparation

TRrest = TR - rf_dur - sys.rfDeadTime - sys.rfRingdownTime ...
            - mr.calcDuration(gxPre) ...
            - mr.calcDuration(gx) ...
            - mr.calcDuration(gxRe);

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
%   meas_i = 1 : flip angle alpha1
%   meas_i = 2 : flip angle alpha2
% Both measurements are written sequentially into a single .seq file.

for meas_i = 1:2

    % Select flip angle for this measurement
    if meas_i == 1
        rf = mr.makeBlockPulse(deg2rad(alpha1), 'Duration', rf_dur, 'System', sys, 'use', 'excitation');
    else
        rf = mr.makeBlockPulse(deg2rad(alpha2), 'Duration', rf_dur, 'System', sys, 'use', 'excitation');
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

                % Simultaneous x/y/z pre-phasers
                seq.addBlock(gxPre, gyPre, gzPre);

                % Readout
                seq.addBlock(gx, adc);

                % Simultaneous x/y/z re-phasers (rewind all encoding gradients)
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

seq.setDefinition('FOV',     fov);
seq.setDefinition('Name',    'B1map_DFA');
seq.setDefinition('mat',     [Nx, Ny, Nz]);
seq.setDefinition('pattern', [len_x, len_y]);
seq.setDefinition('alphas',  [alpha1, alpha2]);

savestr = '\\B1_DFA.seq';
seq.write(savestr);

return