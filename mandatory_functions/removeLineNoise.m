function cleaned = removeLineNoise(signal, Fs, lineFreqs, replaceHalfBW, refHalfBW)
% REMOVELINENOISE  Remove mains line noise by spectral interpolation.
%
%   cleaned = removeLineNoise(signal, Fs)
%   cleaned = removeLineNoise(signal, Fs, lineFreqs, replaceHalfBW, refHalfBW)
%
%  Removes narrow line-noise peaks (default 50/100/150 Hz, European mains +
%  harmonics). For each line frequency it REPLACES the FFT magnitude in a
%  narrow band around the line with a value linearly interpolated between the
%  local baseline on each side, keeping the original phase, then transforms
%  back. Only the narrow contaminated band is altered; the surrounding neural
%  power (low/high gamma) is preserved — no hole carved, unlike a hard notch.
%
%  The baseline on each side is estimated as the MEAN magnitude over a small
%  reference band just outside the replaced band (not a single FFT bin, which
%  is far too noisy in a long recording and produces dips below baseline).
%
%  INPUTS
%    signal        : LFP. A vector (one channel) OR a matrix
%                    [nSamples x nChan] (columns = channels, e.g.
%                    getLFP('all').data). NaNs not allowed.
%    Fs            : sampling rate in Hz (e.g. 1250).
%    lineFreqs     : line frequencies to remove. Default [50 100 150].
%    replaceHalfBW : half-width of the REPLACED band, in Hz. Default 2
%                    (i.e. f0-2 .. f0+2). Widen if the spike is broad.
%    refHalfBW     : width of each side's baseline reference band, in Hz,
%                    sitting just outside the replaced band. Default 2.
%
%  OUTPUT
%    cleaned       : same size/orientation as 'signal', real-valued.
%
%  NOTE: operates on the whole trace at once (global FFT). Assumes the line
%  frequency is stable across the recording, which holds for mains noise.

if nargin < 3 || isempty(lineFreqs),     lineFreqs     = [50 100 150]; end
if nargin < 4 || isempty(replaceHalfBW), replaceHalfBW = 2;            end
if nargin < 5 || isempty(refHalfBW),     refHalfBW     = 2;            end

% ── Matrix input: clean each channel (column) independently ───────────────
if ~isvector(signal)
    cleaned = zeros(size(signal));
    for c = 1:size(signal, 2)
        cleaned(:, c) = removeLineNoise(double(signal(:, c)), Fs, ...
                                        lineFreqs, replaceHalfBW, refHalfBW);
    end
    return
end

wasRow = isrow(signal);
x = double(signal(:));
N = numel(x);

X = fft(x);
f = (0:N-1).' * (Fs / N);     % frequency of each FFT bin (0 .. Fs)

for k = 1:numel(lineFreqs)
    f0 = lineFreqs(k);
    if f0 >= Fs/2, continue; end                      % at/above Nyquist: skip

    fLo = f0 - replaceHalfBW;
    fHi = f0 + replaceHalfBW;

    idx = find(f >= fLo & f <= fHi);                  % bins to replace
    if isempty(idx), continue; end

    % Baseline reference bands just outside the replaced band
    refLoIdx = find(f >= fLo - refHalfBW & f < fLo);
    refHiIdx = find(f >  fHi & f <= fHi + refHalfBW);
    if isempty(refLoIdx), refLoIdx = max(idx(1)   - 1, 1); end
    if isempty(refHiIdx), refHiIdx = min(idx(end) + 1, N); end

    magLeft  = mean(abs(X(refLoIdx)));                % stable baseline, left
    magRight = mean(abs(X(refHiIdx)));                % stable baseline, right

    % Linear interpolation across the replaced band; keep original phase
    tt     = (f(idx) - fLo) / (fHi - fLo);            % 0 at left edge, 1 at right
    newMag = magLeft + tt .* (magRight - magLeft);
    X(idx) = newMag .* exp(1i * angle(X(idx)));
end

cleaned = ifft(X, 'symmetric');   % 'symmetric' guarantees a real result
if wasRow, cleaned = cleaned.'; end
end
