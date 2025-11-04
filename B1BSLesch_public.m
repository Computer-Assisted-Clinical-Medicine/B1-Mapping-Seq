% 3D Cartesian B1map with Bloch siegert method and optional TV
% regularization
% -------------------------------------------------------------------------
% 3D Cartesian B1map with Bloch siegert method and undersampling
% Author: Valentin Jost (2025)
%
% This code is provided for research and educational use only.
% The author assumes no responsibility for its use or any consequences
% arising from running it on MRI hardware.
%
% Users must verify all safety limits (SAR, gradients, duty cycle, etc.)
% before execution on any scanner.
%
% Provided "as is" without warranty of any kind. Use at your own risk.
%
% Based on: Sacolick et al., MRM 2010, DOI: 10.1002/mrm.28284
%           Lesch tel al., MRM 2018, DOI: 10.1002/mrm.27434
% -------------------------------------------------------------------------

clear; clc; close all;

gamma = 11.262e6; %MHz/T %Sodium

sys = mr.opts('B0', 3,'MaxGrad', 45, 'GradUnit', 'mT/m', ...
    'MaxSlew', 200, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'gamma',gamma);

seq=mr.Sequence(sys);           % Create a new sequence object

%Parameters
rf_dur = 500e-6;             
alpha=90; % flip angle
rf = mr.makeBlockPulse(deg2rad(alpha),'Duration',rf_dur,'System',sys);

fov=[240e-3 240e-3 240e-3];  % Define FOV
Nx=40; Ny=40; Nz=40;            %kspace dim
Tread=3e-3;
Tpre=0.8e-3;

TR = 200e-3; %atm do ~+10ms to what you want to achieve

% % Define other gradients and ADC events
deltak=1./fov;
gx = mr.makeTrapezoid('x',sys,'FlatArea',Nx*deltak(1),'FlatTime',Tread);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre);
gxRe = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre);

areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);

% Define 90 degree RF pulse
T = 2e-3;                     
t_res = sys.rfRasterTime;
phi = 0;

%Fermi pulse modulation
rf_amp = fermiwin(T/t_res,T);
rf_amp(1) = 0;
rf_amp(end) = 0;
rf_phase = phi * ones(size(rf_amp)); %constant phase
rf_complex = rf_amp .* exp(1i * rf_phase);

%Fermi amplitude calibration
frac = 1;
frac_peak = frac*(T/4e-3);
B1peak = frac_peak*7.41; %7.41 is the multiplier for a 90° pulse for sodium in 500e-6s and 4ms fermi

rf_blosi = mr.makeArbitraryRf(rf_complex,B1peak,'System', sys,'use','other'); 

% Define frequency offset for Bloch-Siegert
bs_offset = 1e3; % 4 kHz off-resonance

gamma_rad = gamma*2*pi; % rad/s/T
gamma_gauss = gamma_rad*1e-4; %rad/s/gauss %1T=10^4 gauss
Kbs = sum(rf_amp.^2) * (gamma_gauss)^2 * t_res / (4*pi*bs_offset) / (max(rf_amp))^2;

%prepare delays
TRrest = (TR - rf_dur - T - 2*Tpre - Tread);
delayTR=mr.makeDelay(TRrest);

%define subsampling mask
len_x = 12;
len_y = 12;
sampling_mask = block2D_mask(Nx,Ny,len_x,len_y);

for offres_i = 1:2
    if offres_i == 1
        rf_blosi.freqOffset = bs_offset;
    else
        rf_blosi.freqOffset = -bs_offset;
    end

    for iZ=1:Nz
        gzPre = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',Tpre);
        gzRe = mr.makeTrapezoid('z','Area',-areaZ(iZ),'Duration',Tpre);

        for iY=1:Ny
            if sampling_mask(iY,iZ)
                gyPre = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);
                gyRe = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',Tpre);
    
                seq.addBlock(rf)
                seq.addBlock(mr.makeDelay(100e-6));
                seq.addBlock(rf_blosi)
                
                % Encoding
                seq.addBlock(gxPre,gyPre,gzPre);
                
                seq.addBlock(gx,adc);
    
                seq.addBlock(gxRe,gyRe,gzRe);
                    
                seq.addBlock(delayTR)
            end    
        end
    end
end

%% check whether the timing of the sequence is correct
fprintf('Timing check...')
[ok, error_report]=seq.checkTiming;

if (ok)
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
%% write sequence
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'B1map_BS');
seq.setDefinition('mat', [Nx Ny Nz]);
seq.setDefinition('Kbs',Kbs);
seq.setDefinition('mask',sampling_mask);
seq.setDefinition('pattern',pattern);
seq.setDefinition('pat_param',[len_x,len_y]);

savestr = sprintf('\\BSSacc_%dms_%dx%d_%dTR_%dFoV.seq',Tread*1e3,len_x,len_y,TR*1e3,fov(1)*1e3);
seq.write(savestr);

return


% methods
function mask = block2D_mask(Nx,Ny,len_x,len_y)
    x_center = floor(Nx/2);
    y_center = floor(Ny/2);
    halflen_x = floor(len_x/2);
    halflen_y = floor(len_y/2);

    mask = zeros(Nx,Ny);
    mask(x_center-halflen_x+1:x_center+halflen_x, ...
        y_center-halflen_y+1:y_center+halflen_y) = 1;
end

function [W,varargout] = fermiwin(L,t_length,alpha)
    if nargin < 3
        alpha = 1/33.81;
    end

    if nargin < 2
        t_length = 1;
    end

    ep = linspace(0,1,L)-0.5;
    t0 = alpha*10;
    F1 = 1 ./ (exp((abs(ep)-t0)/alpha) + 1);
    W = [F1];
    if nargin > 1
        varargout{1} = linspace(0,t_length,L);
    end
end