function [signal_corrected]= TDDR(signal,sample_rate)
%% retirado de https://github.com/frankfishburn/TDDR/blob/master/TDDR.m
%% Artigo: Fishburn, F. A., Ludlum, R. S., Vaidya, C. J., & Medvedev, A. V. (2019).
%Temporal Derivative Distribution Repair (TDDR): A motion correction method for fNIRS. NeuroImage, 184, 171–179.
%https://doi.org/10.1016/j.neuroimage.2018.09.025
%% Iterate over each channel
nch = size(signal,2);
if nch>1
    signal_corrected = zeros(size(signal));
    for ch = 1:nch
        signal_corrected(:,ch) = TDDR( signal(:,ch) , sample_rate );
    end
    return
end

%% Preprocess: Separate high and low frequencies
filter_cutoff = .5;
filter_order = 3;
Fc = filter_cutoff * 2/sample_rate;
signal_mean = mean(signal);
signal = signal - signal_mean;
if Fc<1
    [fb,fa] = butter(filter_order,Fc);
    signal_low = filtfilt(fb,fa,signal);
else
    signal_low = signal;
end
signal_high = signal - signal_low;

%% Initialize
tune = 4.685;
D = sqrt(eps(class(signal)));
mu = inf;
iter = 0;

%% Step 1. Compute temporal derivative of the signal
deriv = diff(signal_low);

%% Step 2. Initialize observation weights
w = ones(size(deriv));

%% Step 3. Iterative estimation of robust weights
while iter < 50
    
    iter = iter + 1;
    mu0 = mu;
    
    % Step 3a. Estimate weighted mean
    mu = sum( w .* deriv ) / sum( w );
    
    % Step 3b. Calculate absolute residuals of estimate
    dev = abs(deriv - mu);

    % Step 3c. Robust estimate of standard deviation of the residuals
    sigma = 1.4826 * median(dev);

    % Step 3d. Scale deviations by standard deviation and tuning parameter
    r = dev / (sigma * tune);
    
    % Step 3e. Calculate new weights accoring to Tukey's biweight function
    w = ((1 - r.^2) .* (r < 1)) .^ 2;

    % Step 3f. Terminate if new estimate is within machine-precision of old estimate
    if abs(mu-mu0) < D*max(abs(mu),abs(mu0))
        break;
    end

end

%% Step 4. Apply robust weights to centered derivative
new_deriv = w .* (deriv-mu);

%% Step 5. Integrate corrected derivative
signal_low_corrected = cumsum([0; new_deriv]);

%% Postprocess: Center the corrected signal
signal_low_corrected = signal_low_corrected - mean(signal_low_corrected);

%% Postprocess: Merge back with uncorrected high frequency component
signal_corrected = signal_low_corrected + signal_high + signal_mean;

%% Plot da derivada, derivada corrigida e ponderações
figure;
subplot(3,1,1)
plot(deriv(:,1));
ylim([-2 2])
xlabel('Tempo');
ylabel('Derivada');
title('Derivada ');
subplot(3,1,2)
plot(new_deriv(:,1));
ylim([-2 2])
xlabel('Tempo');
ylabel('Derivada');
title('Derivada corrigida ');
subplot(3,1,3)
plot(w(:,1))
ylim([0 1.5])
xlabel('Tempo');
ylabel('Ponderação');
title('Ponderações');
  
end

