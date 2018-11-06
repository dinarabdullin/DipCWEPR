<<<<<<< HEAD
function f = fit_spectrum(v, pf, Const)
    
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
    figure(1);
    subplot(1,3,1);
    plot(rv,pr,'black');
    xlabel('Inter-spin distance / nm'); 
    ylabel('Probability density'); 
    xlim([0 3]); 
    ylim([-0.1 1.2]);
    title('Inter-spin distance distribution');
    %hist(r, min(r):0.01:max(r))

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
    subplot(1,3,2);
    plot(xdd,ydd,'black')
    xlabel('Magnetic field offset / mT'); 
    ylabel('Relative Intensity');
    xlim([min(x) max(x)]); 
    ylim([-0.1 1.2]);
    title('Dipolar spectrum');
    
    % Convolute the reference spectrum with the dipolar spectrum
    yc = conv(yr,ydd,'same');
    yc = yc / max(yc);
    
    % Calculate the sum of unbroadened and broadened spectra
    ir = cumtrapz(x,cumtrapz(x,yr));
    ic = cumtrapz(x,cumtrapz(x,yc));
    yf = yr * (1-weight) / ir(end) + yc * weight / ic(end);
    yf = yf / max(yf);
    subplot(1,3,3);
    plot(x,yr,'blue',x,ys,'black',x,yf,'red')
    xlabel('Magnetic field offset / mT'); 
    ylabel('Relative Intensity'); 
    xlim([min(x) max(x)]);
    ylim([-0.1 0.1]);
    legend('reference','sample','fit');
    title('Fitting of dipolar-broadened spectrum');
    
    % Calculate RMSD between the sample spectrum and the convolution
    f = sqrt(sum((yf - ys).^2) / size(ys,1));
end

=======
function f = fit_spectrum(v, pf, Const)
    
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
    figure(1);
    subplot(1,3,1);
    plot(rv,pr,'black');
    xlabel('Inter-spin distance / nm'); 
    ylabel('Probability density'); 
    xlim([0 3]); 
    ylim([-0.1 1.2]);
    title('Inter-spin distance distribution');
    %hist(r, min(r):0.01:max(r))

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
    subplot(1,3,2);
    plot(xdd,ydd,'black')
    xlabel('Magnetic field offset / mT'); 
    ylabel('Relative Intensity');
    xlim([min(x) max(x)]); 
    ylim([-0.1 1.2]);
    title('Dipolar spectrum');
    
    % Convolute the reference spectrum with the dipolar spectrum
    yc = conv(yr,ydd,'same');
    yc = yc / max(yc);
    
    % Calculate the sum of unbroadened and broadened spectra
    ir = cumtrapz(x,cumtrapz(x,yr));
    ic = cumtrapz(x,cumtrapz(x,yc));
    yf = yr * (1-weight) / ir(end) + yc * weight / ic(end);
    yf = yf / max(yf);
    subplot(1,3,3);
    plot(x,yr,'blue',x,ys,'black',x,yf,'red')
    xlabel('Magnetic field offset / mT'); 
    ylabel('Relative Intensity'); 
    xlim([min(x) max(x)]);
    ylim([-0.1 0.1]);
    legend('reference','sample','fit');
    title('Fitting of dipolar-broadened spectrum');
    
    % Calculate RMSD between the sample spectrum and the convolution
    f = sqrt(sum((yf - ys).^2) / size(ys,1));
end

>>>>>>> 85f612edeb1d7cf8e7fa5aa956a3cdc503a58eec
