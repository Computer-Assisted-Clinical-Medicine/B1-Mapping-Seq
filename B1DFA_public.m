% -------------------------------------------------------------------------
% 3D Cartesian B1map with Double Flip Angle method
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
% Based on: Stollberger et al., MRM 1996, DOI: 10.1002/mrm.1910350217
% -------------------------------------------------------------------------

clear; clc; close all;

gamma = 11.262e6; %MHz/T %Sodium

sys = mr.opts('B0', 3,'MaxGrad', 45, 'GradUnit', 'mT/m', ...
    'MaxSlew', 200, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'gamma',gamma);

seq=mr.Sequence(sys); % Create a new sequence object

%Parameters
rf_dur = 500e-6;             
alpha1 = 60; %degree
alpha2 = 120; %degree

fov=[240e-3 240e-3 240e-3];  % Define FOV
Nx=40; Ny=40; Nz=40; %kspace dims
Tread=3e-3;
Tpre=0.8e-3;

%Full relaxation is necessary, even after 120° pulse!
TR = 250e-3; %atm do +10ms to what you want to achieve

% % Define other gradients and ADC events
deltak=1./fov;
gx = mr.makeTrapezoid('x',sys,'FlatArea',Nx*deltak(1),'FlatTime',Tread);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre);
gxRe = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre);

areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);


%prepare delays
TRrest = (TR - rf_dur - 2*Tpre - Tread);
delayTR=mr.makeDelay(TRrest);

len_x = 14;
len_y = 14;
sampling_mask = block2D_mask(Nx,Ny,len_x,len_y);

for meas_i = 1:2
    if meas_i == 1
        rf = mr.makeBlockPulse(deg2rad(alpha1),'Duration',rf_dur,'System',sys,'use','excitation');
    else
        rf = mr.makeBlockPulse(deg2rad(alpha2),'Duration',rf_dur,'System',sys,'use','excitation');
    end

    for iZ=1:Nz
        gzPre = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',Tpre);
        gzRe = mr.makeTrapezoid('z','Area',-areaZ(iZ),'Duration',Tpre);

        for iY=1:Ny
            if sampling_mask(iY,iZ)
                gyPre = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);
                gyRe = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',Tpre);
    
                seq.addBlock(rf)
                
                % Encoding
                seq.addBlock(gxPre,gyPre,gzPre);
                
                seq.addBlock(gx,adc);
    
                seq.addBlock(gxRe,gyRe,gzRe);
                    
                seq.addBlock(delayTR)
            end
        end
    end
end

return
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
seq.setDefinition('Name', 'B1map_DFA');
seq.setDefinition('mat', [Nx Ny Nz]);
seq.setDefinition('pattern',[len_x, len_y]);
seq.setDefinition('alphas',[alpha1,alpha2]);

savestr = sprintf('B1DFA_%dms_%dx%d_%dTR_%dFoV.seq',Tread*1e3,len_x,len_y,TR*1e3,fov(1)*1e3);
seq.write(savestr);

return

%% methods
function mask = block2D_mask(Nx,Ny,len_x,len_y)
    x_center = floor(Nx/2);
    y_center = floor(Ny/2);
    halflen_x = floor(len_x/2);
    halflen_y = floor(len_y/2);

    mask = zeros(Nx,Ny);
    mask(x_center-halflen_x+1:x_center+halflen_x, ...
        y_center-halflen_y+1:y_center+halflen_y) = 1;
end