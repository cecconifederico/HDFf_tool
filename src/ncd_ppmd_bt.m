function [Lk, info] = ncd_ppmd_bt(X, ks, n_bits, Fs, fc_Hz, reportBase)
% NCD_PPMD_BT Normalized Compression Distance NCD(P,Qk) using the PPMd
% compressor (variant H), with causal Butterworth anti-aliasing filtering
% before decimation.
%
%   [Lk, info] = ncd_ppmd_bt(X, ks, n_bits, Fs, fc_Hz, reportBase)
%
% For each decimation factor k in ks:
%   1) Y = 7th-order causal Butterworth low-pass (cutoff fc_Hz(i)) applied
%      to the full signal X
%   2) q = Y(1:k:end)                         (decimation)
%   3) if reportBase: Q = each sample of q repeated k times and truncated to
%      the length of X, so that P and Q share the same time base and length;
%      otherwise Q = q
%   4) P and Q are quantized to n_bits (z-score, clipped at +/-4 SD) and
%      stored as bytes
%   5) NCD(P,Q) = [C(P+Q) - min(C(P),C(Q))] / max(C(P),C(Q)),
%      where C(.) is the size in bytes compressed by PPMd (variant H, order 6)
% P is the quantized, unfiltered signal X.
%
% No fit and no fractal dimension: the NCD value is returned for each k.
%
% INPUT (all optional except X)
%   X            signal (vector)
%   ks           decimation factors            (default [2 5 10 20 33 67])
%   n_bits       quantization levels 2^n_bits  (default 3)
%   Fs           sampling frequency in Hz      (default 200)
%   fc_Hz        cutoff for each k, in Hz      (default [35 14 7 3.5 2 1])
%   reportBase   true/false                    (default false)
%
% OUTPUT
%   Lk    1 x numel(ks) vector: Lk(i) = NCD for ks(i)
%   info  struct with the parameters used
%
% Requires: bw_lowpass_causal.m, and a pyenv with Python 3.9-3.12 and
%           pyppmd installed (pip install pyppmd).

    if nargin < 2 || isempty(ks),         ks = [2 5 10 20 33 67]; end
    if nargin < 3 || isempty(n_bits),     n_bits = 3; end
    if nargin < 4 || isempty(Fs),         Fs = 200; end
    if nargin < 5 || isempty(fc_Hz),      fc_Hz = [35 14 7 3.5 2 1]; end
    if nargin < 6 || isempty(reportBase), reportBase = false; end

    if numel(fc_Hz) ~= numel(ks)
        error('ncd_ppmd_bt: fc_Hz must have the same length as ks.');
    end

    variant   = 'H';
    max_order = int32(6);

    X = X(:);
    X = X(isfinite(X));
    N = numel(X);

    info = struct('N', N, 'ks', ks, 'fc_Hz', fc_Hz, 'Fs', Fs, 'n_bits', n_bits, ...
                  'reportBase', reportBase, 'variant', variant, 'max_order', double(max_order));

    Lk = NaN(1, numel(ks));
    if N < 20
        return;
    end

    % Quantize and compress the full signal once
    Pq = local_quantize(X, n_bits);
    cP = local_ppmd_csize(Pq, variant, max_order);
    info.cP = cP;

    for i = 1:numel(ks)
        k = ks(i);
        if N < 2*k
            continue;
        end

        Y = bw_lowpass_causal(X, fc_Hz(i), Fs, 7);
        q = Y(1:k:end);

        if reportBase
            Q = repelem(q, k);
            Q = Q(1:N);
        else
            Q = q;
        end

        Qq  = local_quantize(Q, n_bits);
        cQ  = local_ppmd_csize(Qq, variant, max_order);
        cPQ = local_ppmd_csize([Pq; Qq], variant, max_order);

        Lk(i) = (cPQ - min(cP, cQ)) / max(cP, cQ);
    end
end


function b = local_quantize(x, n_bits)
% z-score, clip at +/-4 SD, scale to 2^n_bits levels (uint8, 0..2^n_bits-1)
    n_std = 4;
    x  = double(x(:));
    sd = std(x);
    if sd == 0, sd = 1; end
    z  = (x - mean(x)) / sd;
    z  = max(min(z, n_std), -n_std);
    b  = uint8(round((z + n_std) / (2*n_std) * (2^n_bits - 1)));
end


function c = local_ppmd_csize(uint8vec, variant, max_order)
% Compressed size (bytes) of a uint8 vector using PPMd via pyppmd
    pyBytesObj = py.bytes(uint8vec);
    comp = py.pyppmd.compress(pyBytesObj, pyargs('variant', variant, 'max_order', max_order));
    c = double(py.len(comp));
end
