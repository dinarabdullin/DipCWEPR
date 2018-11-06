clear
clc

%% Constants
Const.deg2rad = pi / 180;
Const.rad2deg = 180 / pi;
Const.Hz2MHz = 1e-6;
Const.MHz2Hz = 1e6;
Const.Hz2GHz = 1e-9;
Const.mT2T = 1e-3;
Const.T2mT = 1e3;
Const.nm2m = 1e-9;
Const.MHz2mT = 0.1 / 2.8; % in mT/MHz
Const.plankConstant = 6.626070040e-34; % in J*s
Const.bohrMagneton = 9.274009994e-24; % in J/T
Const.ge = 2.0023; % free electron g-factor
Const.vacuumPermeabilityOver4Pi = 1e-7; % in T*m/A
Const.Fez = Const.plankConstant * Const.MHz2Hz * Const.T2mT / (Const.ge * Const.bohrMagneton); % in mT/MHz
Const.Fdd = Const.Hz2MHz * Const.vacuumPermeabilityOver4Pi * Const.ge^2 * ...
            Const.bohrMagneton^2 / (Const.plankConstant * Const.nm2m^3); % in MHz

%% Load the spectrum       
[filename, path] = uigetfile('*.DTA','Load the spectrum');
if ~isequal(filename,0)
    fnSample = strcat(path, filename);
    [xs, ys, ps] = eprload(fnSample);
else
	error('No file was selected!')
end

%% Load the reference spectrum       
[filename, path] = uigetfile('*.DTA','Load the reference spectrum');
if ~isequal(filename,0)
    fnRef = strcat(path, filename);
    [xr, yr, pr] = eprload(fnRef);
else
	error('No file was selected!')
end

%% Optimization parameters
title = 'Enter optimization parameters';
lines = 1;
prompt = {
    'Simulation (0) or fitting (1)?', ...
    'Strong coupling (0) or weak coupling (1)?', ...
    'Mean distance / nm', ...
    'Standard deviation / nm', ...
    'Weigth', ...
    'Mumber of MC samples', ...
    'Magnetic field increment / mT', ...
    'Microwave frequency / GHz'};
def = {'0', '0', '0.63', '0.04', '0.25', '1000000', '0.01', '9.440000'};
answer = inputdlg(prompt, title, lines, def);
mode = str2num([answer{1}]);
regime = str2num([answer{2}]);
rMean = str2num([answer{3}]);
rStd = str2num([answer{4}]);
weight = str2num([answer{5}]);
nSamples = str2num([answer{6}]);
dx = str2num([answer{7}]);
mwFreq = str2num([answer{8}]);

%% Some pre-calculations
% Extract the real parts of the experimental spectra
ys = real(ys);
yr = real(yr);
% Normalize the experimental spectra
ys = ys / max(ys);
yr = yr / max(yr);
% Convert the units of the field from Gauss to mT
xs = 0.1 * xs;
xr = 0.1 * xr;
% Shift the spectra of the sample and reference such that they correspond to the same m.w. frequency
xs = xs .* mwFreq / (ps.MWFQ * Const.Hz2GHz);
xr = xr .* mwFreq / (pr.MWFQ * Const.Hz2GHz);
% Set the min and max values of the field
xMin = xs(1);
xMax = xs(end);
if xr(1) > xMin
    xMin = xr(1);
end
if xr(end) < xMax
    xMax = xr(end);
end
% Set the linear field grid
xc = xMin:dx:xMax;
Nx = size(xc,2);
% Interpolate the spectra to fit the new field grid
ysc = interp1(xs,ys,xc,'linear');
yrc = interp1(xr,yr,xc,'linear');
% Shift the centrum of spectra to 0
intg = cumtrapz(xc,ysc);
[intgMax, intgIdx] = max(intg);
xshift = xc(intgIdx); 
xc0 = xc - xshift;
if (intgIdx <= 0.5*Nx)
    xc0 = xc0(1:(2*intgIdx - 1));
    ysc0 = ysc(1:(2*intgIdx - 1));
    yrc0 = yrc(1:(2*intgIdx - 1));
else
    xc0 = xc0((2*intgIdx - 1 - Nx):end);
    ysc0 = ysc((2*intgIdx - 1 - Nx):end);
    yrc0 = yrc((2*intgIdx - 1 - Nx):end);
end

%% Simulate the spectrum
if mode == 0
    v = [rMean, rStd, weight];
    pf.x = xc0;
    pf.ys = ysc0;
    pf.yr = yrc0;
    pf.regime = regime;
    pf.nSamples = nSamples;
    rmsd = fit_spectrum(v, pf, Const);
    fprintf('RMSD: %.4f\n', rmsd);
end

%% Fit the spectrum
if mode == 1
    nopt = 2;
    vopt = zeros([nopt, 3]);
    rmsdopt = zeros([nopt, 1]);
    for i=1:nopt
        v = [rMean, rStd, weight];
        pf.x = xc0;
        pf.ys = ysc0;
        pf.yr = yrc0;
        pf.regime = regime;
        pf.nSamples = nSamples;
        options = optimset('Display','iter', 'FunValCheck', 'on', ... 
                           'MaxFunEvals', 200*size(v,2), 'MaxIter', 200*size(v,2), ... 
                           'OutputFcn', [], 'PlotFcns', @optimplotfval, ...
                           'TolFun', 1e-4, 'TolX', 1e-4);
        [a,rmsd] = fminsearch(@(v) fit_spectrum(v, pf, Const), v, options);
        fprintf('Mean distance: %.2f nm\n', a(1));
        fprintf('Standard deviation: %.2f nm\n', a(2));
        fprintf('Weight: %.2f\n', a(3));
        vopt(i,:) = a;
        rmsdopt(i) = rmsd;
    end
    % Find best fit
    [M,I] = min(rmsdopt);
    fprintf('\nBest mean inter-spin distance: %.2f nm\n', vopt(I(1),1));
    fprintf('Best standard deviation: %.2f nm\n', vopt(I(1),2));
    fprintf('Best relative weight: %.2f\n', vopt(I(1),3));
    v = vopt(I(1),:);
    % Save the best fit
    save_fit(v, pf, Const, xshift, fnSample);
end   