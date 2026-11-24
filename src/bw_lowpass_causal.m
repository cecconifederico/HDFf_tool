function y = bw_lowpass_causal(x, fc, Fs, order)
% BW_LOWPASS_CAUSAL Causal Butterworth low-pass filter (single pass).
%
%   y = bw_lowpass_causal(x, fc, Fs, order)
%
% Implemented as cascaded second-order sections (SOS) for numerical
% stability at cutoff frequencies much lower than Fs. The filter is causal
% (no filtfilt): y(n) depends only on x(n), x(n-1), ...
% The initial state is set to steady state for a constant input equal to
% x(1), which removes the start-up transient without using future samples.
%
% INPUT
%   x     signal (vector)
%   fc    cutoff frequency (-3 dB), in Hz
%   Fs    sampling frequency, in Hz
%   order filter order (default 7)

    if nargin < 4 || isempty(order), order = 7; end
    x = double(x(:));
    Wn = fc / (Fs/2);
    if ~(Wn > 0 && Wn < 1)
        error('bw_lowpass_causal: invalid fc=%g Hz for Fs=%g Hz.', fc, Fs);
    end

    [z, p, kg] = butter(order, Wn, 'low');
    [sos, g] = zp2sos(z, p, kg);

    y = g * x;
    for s = 1:size(sos, 1)
        b = sos(s, 1:3);
        a = sos(s, 4:6);
        y0  = y(1);
        yss = y0 * sum(b) / sum(a);               % steady-state output for constant input y0
        zi  = filtic(b, a, [yss yss], [y0 y0]);   % matching initial state
        y   = filter(b, a, y, zi);
    end
end
