function f = save_fit(v, pf, Const, xshift, fn)
    
    % Read input parameters
    x = pf.x;
    ys = pf.ys;
    yr = pf.yr;
    regime = pf.regime;
    nSamples = pf.nSamples;
    rMean = v(1);
    rStd = v(2);
    weight = v(3);
    
    % Simulate the normal distribution of inter-spin distances
    r = rMean + rStd .* randn(nSamples,1);    
    [pr, rv] = hist(r, 0:0.01:3);
    pr(1) = 0; pr(end) = 0;
    pr = pr / max(pr);

    % Simulate the cos-weighted distribution of theta
    theta = acos(rand(nSamples,1));
    
    % Calculate the dipolar spectrum
    if (regime == 1)
        fdd = (1/2 * Const.Fdd) .* (1 - 3 * cos(theta).^2) ./ r.^3;
    elseif (regime == 0)
        fdd = (3/4 * Const.Fdd) .* (1 - 3 * cos(theta).^2) ./ r.^3;
    end
    bdd = Const.Fez .* fdd;
    [ydd, xdd] = hist(bdd, x);
    ydd(1) = 0; ydd(end) = 0;
    ydd = ydd + fliplr(ydd);
    ydd = ydd / max(ydd);
    
    % Convolute the reference spectrum with the dipolar spectrum
    yc = conv(yr,ydd,'same');
    yc = yc / max(yc);
    
    % Calculate the sum of unbroadened and broadened spectra
    ir = cumtrapz(x,cumtrapz(x,yr));
    ic = cumtrapz(x,cumtrapz(x,yc));
    yf = yr * (1-weight) / ir(end) + yc * weight / ic(end);
    yf = yf / max(yf);
    
    % Save the data
    x = x + xshift;
    M = [x.' yr.' ys.' yf.'];
    dlmwrite(strcat(fn(1:end-4),'_fit.dat'),M,'delimiter','\t','precision',6);
end

