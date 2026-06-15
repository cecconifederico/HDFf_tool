function [hfd, Lk, fitInfo] = fractaldim_frechet(X, kmax, alpha)
%FRACTALDIM_FRECHET Estimate HFDf using the discrete Frechet distance.

%   This function computes a modified Higuchi-like fractal dimension by
%   replacing the Euclidean distance with the discrete Frechet distance.

    if nargin < 3 || isempty(alpha)
        alpha = 2;
    end

    kmin = 2;
    minPtsFit = 8;

    X = X(:);
    X = X(isfinite(X));
    N = numel(X);

    if N < 20
        hfd = NaN;
        Lk = [];
        fitInfo = struct();
        return;
    end

    kmax = min(kmax, floor(N/2));
    Lk = NaN(1, kmax);

    P = X;

    for k = kmin:kmax

        Qk = X(1:k:end);

        if numel(Qk) < 2
            continue;
        end

        dF = DiscreteFrechetDist(P, Qk);
        Lk(k) = dF / (k^alpha);

        if Lk(k) <= 0
            Lk(k) = NaN;
        end
    end

    kk = (kmin:kmax)';
    yy = Lk(kk)';

    valid = isfinite(yy) & (yy > 0);

    if nnz(valid) < minPtsFit
        hfd = NaN;
        fitInfo = struct();
        return;
    end

    xfit = log(1 ./ kk(valid));
    yfit = log(yy(valid));

    coeff = polyfit(xfit, yfit, 1);
    hfd = coeff(1);

    yhat = polyval(coeff, xfit);

    SSres = sum((yfit - yhat).^2);
    SStot = sum((yfit - mean(yfit)).^2);

    R2 = 1 - SSres / (SStot + eps);

    fitInfo = struct( ...
        'coeff', coeff, ...
        'R2', R2, ...
        'alpha', alpha, ...
        'k_used', kk(valid), ...
        'xfit', xfit, ...
        'yfit', yfit);

end


function d = DiscreteFrechetDist(P, Q)

    P = P(:);
    Q = Q(:);

    nP = numel(P);
    nQ = numel(Q);

    if nP == 0 || nQ == 0
        d = 0;
        return;
    end

    CA = inf(nP, nQ);

    CA(1,1) = abs(P(1) - Q(1));

    for i = 2:nP
        CA(i,1) = max(CA(i-1,1), abs(P(i) - Q(1)));
    end

    for j = 2:nQ
        CA(1,j) = max(CA(1,j-1), abs(P(1) - Q(j)));
    end

    for i = 2:nP
        Pi = P(i);

        for j = 2:nQ
            dij = abs(Pi - Q(j));

            CA(i,j) = max( ...
                min([CA(i-1,j), CA(i-1,j-1), CA(i,j-1)]), ...
                dij);
        end
    end

    d = CA(nP, nQ);

end
