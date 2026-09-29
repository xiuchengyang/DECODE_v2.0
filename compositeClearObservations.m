function Tout = compositeClearObservations(T, maxN)
% compositeClearObservations
% Select at most maxN clear observations per year using greedy
% day-of-year spacing.
%
% Duplicate observations on the same year + day-of-year are reduced to
% one randomly selected record.
%
% Inputs:
%   T     - input table with variable 'date'
%   maxN  - maximum number of observations per year, e.g. 18
%
% Output:
%   Tout  - selected table

    if nargin < 2 || isempty(maxN)
        maxN = 18;
    end

    if ~isdatetime(T.date)
        T.date = datetime(T.date);
    end

    % Add temporary year and day-of-year columns
    T.tmpYear = year(T.date);
    T.tmpDOY = day(T.date, 'dayofyear');

    % Randomize row order so duplicate records on the same day
    % are randomly reduced to one record
    T = T(randperm(height(T)), :);

    % Keep one random observation per year + day-of-year
    [~, ia] = unique(T(:, {'tmpYear', 'tmpDOY'}), 'stable');
    T = T(ia, :);

    uniqueYears = unique(T.tmpYear);

    keepIdxAll = [];

    for i = 1:numel(uniqueYears)

        thisYear = uniqueYears(i);
        idx = find(T.tmpYear == thisYear);
        Ty = T(idx, :);

        n = height(Ty);

        if n <= maxN
            keepIdxAll = [keepIdxAll; idx];
            continue
        end

        doy = double(Ty.tmpDOY(:));

        selected = false(n, 1);

        % Select earliest and latest observation first
        [~, iMin] = min(doy);
        [~, iMax] = max(doy);

        selected(iMin) = true;
        selected(iMax) = true;

        % Greedy selection:
        % repeatedly select the observation farthest from
        % the nearest already-selected observation
        while sum(selected) < maxN

            remainingIdx = find(~selected);
            selectedIdx = find(selected);

            score = zeros(numel(remainingIdx), 1);

            for j = 1:numel(remainingIdx)

                r = remainingIdx(j);

                distToSelected = abs(doy(r) - doy(selectedIdx));
                score(j) = min(distToSelected);

            end

            [~, bestJ] = max(score);
            selected(remainingIdx(bestJ)) = true;

        end

        keepIdxAll = [keepIdxAll; idx(selected)];

    end

    Tout = T(keepIdxAll, :);

    % Remove temporary columns
    Tout.tmpYear = [];
    Tout.tmpDOY = [];

    % Sort output by date
    Tout = sortrows(Tout, 'date');

end