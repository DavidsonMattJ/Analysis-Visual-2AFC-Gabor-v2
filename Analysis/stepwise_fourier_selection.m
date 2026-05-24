function results = stepwise_fourier_selection(X, xvec, min_series, hz_range, predropped)
% STEPWISE_FOURIER_SELECTION
% Finds the best subset of time-series for a fourier1 group-level fit
%
% INPUTS:
%   X            - [40 x T] matrix of time series (rows = series, cols = time)
%   xvec         - time vector, e.g. 1:10
%   min_series   - minimum series to retain (default: 30)
%   hz_range     - [lo hi] acceptable Hz window (default: [1.5 2.5])
%   predropped   - indices to exclude before selection begins (default: [])
%                  e.g. [1 2] removes series 1 and 2 from all consideration
%
% OUTPUTS:
%   results      - struct with best subset, R² trajectory, Hz estimates, flags
%
% EXAMPLE:
%   results = stepwise_fourier_selection(X, 1:10, 30, [1.5 2.5], [1 2]);

    % --- Defaults ---
    if nargin < 3, min_series  = 30;        end
    if nargin < 4, hz_range    = [1.5 2.5]; end
    if nargin < 5, predropped  = [];        end

    [n_series, ~] = size(X);

    % Validate predropped indices
    if ~isempty(predropped)
        predropped = unique(predropped(:)');  % ensure row vector, no duplicates
        assert(all(predropped >= 1 & predropped <= n_series), ...
            'predropped contains indices outside valid range [1, %d]', n_series);
        fprintf('Pre-dropped series (excluded from selection): %s\n', ...
            mat2str(predropped));
    end

    % Pool available to the algorithm — everything except pre-dropped
    available  = setdiff(1:n_series, predropped);
    n_avail    = numel(available);
    n_steps    = n_avail - min_series;  % how many the algorithm can drop

    assert(n_avail >= min_series, ...
        'Only %d series available after pre-dropping; cannot retain %d.', ...
        n_avail, min_series);

    fprintf('Available series: %d | Pre-dropped: %d | Steps to run: %d\n\n', ...
        n_avail, numel(predropped), n_steps);

    % --- Initialise ---
    retained      = available;   % start with all available
    eliminated    = [];          % dropped by the algorithm (swap-eligible)
    dropped_order = zeros(1, n_steps);

    r2_trajectory  = zeros(1, n_steps + 1);
    hz_trajectory  = zeros(1, n_steps + 1);
    hz_flag        = false(1,  n_steps + 1);
    fit_trajectory = cell(1,   n_steps + 1);

    % Baseline fit — all available series
    [r2_trajectory(1), hz_trajectory(1), fit_trajectory{1}] = ...
        compute_fit(X, retained, xvec);
    hz_flag(1) = is_out_of_range(hz_trajectory(1), hz_range);

    fprintf('Step 0 (%d available series): R² = %.4f | Hz = %.3f%s\n', ...
        n_avail, r2_trajectory(1), hz_trajectory(1), flag_str(hz_flag(1)));

    % ================================================================
    % PHASE 1: Backward elimination
    % ================================================================
    for step = 1:n_steps

        n_retained    = numel(retained);
        r2_candidates = nan(1, n_series);
        hz_candidates = nan(1, n_series);

        % Try removing each retained series one at a time
        for k = 1:n_retained
            candidate = retained([1:k-1, k+1:n_retained]);
            [r2_candidates(retained(k)), hz_candidates(retained(k))] = ...
                compute_fit(X, candidate, xvec);
        end

        % Best removal = highest R²
        [best_r2_step, best_drop_idx] = max(r2_candidates);
        best_hz_step                  = hz_candidates(best_drop_idx);

        % Update lists
        retained(retained == best_drop_idx) = [];
        eliminated(end+1)   = best_drop_idx;            %#ok<AGROW>
        dropped_order(step) = best_drop_idx;

        r2_trajectory(step+1) = best_r2_step;
        hz_trajectory(step+1) = best_hz_step;
        hz_flag(step+1)       = is_out_of_range(best_hz_step, hz_range);
        [~, ~, fit_trajectory{step+1}] = compute_fit(X, retained, xvec);

        fprintf('Step %d: Dropped series %2d | R² = %.4f | Hz = %.3f%s | Retained: %d\n', ...
            step, best_drop_idx, best_r2_step, best_hz_step, ...
            flag_str(hz_flag(step+1)), numel(retained));
    end

    % ================================================================
    % PHASE 2: Swap refinement
    % ================================================================
    % Note: only algorithm-eliminated series are swap candidates.
    % Pre-dropped series are never reintroduced.
    fprintf('\n--- Swap refinement (algorithm-eliminated series only) ---\n');
    improved   = true;
    swap_count = 0;

    while improved
        improved = false;

        for di = 1:numel(eliminated)
            for ri = 1:numel(retained)

                candidate = retained;
                candidate(ri) = eliminated(di);

                [r2_swap, hz_swap] = compute_fit(X, candidate, xvec);
                [r2_curr, ~]       = compute_fit(X, retained,  xvec);

                if r2_swap > r2_curr
                    out_flag = is_out_of_range(hz_swap, hz_range);
                    fprintf('  Swapped IN %2d, OUT %2d | R² %.4f -> %.4f | Hz = %.3f%s\n', ...
                        eliminated(di), retained(ri), r2_curr, r2_swap, ...
                        hz_swap, flag_str(out_flag));

                    old_retained     = retained(ri);
                    retained(ri)     = eliminated(di);
                    eliminated(di)   = old_retained;
                    swap_count       = swap_count + 1;
                    improved         = true;
                end
            end
        end
    end

    if swap_count == 0
        fprintf('  No swaps improved fit.\n');
    else
        fprintf('  %d swap(s) accepted.\n', swap_count);
    end

    % ================================================================
    % Final fit on best subset
    % ================================================================
    [best_r2, best_hz, best_fit, best_gof] = compute_fit(X, retained, xvec);
    best_flag = is_out_of_range(best_hz, hz_range);

    fprintf('\n--- Final Result ---\n');
    fprintf('Pre-dropped series : %s\n', mat2str(predropped));
    fprintf('Retained series    : %s\n', mat2str(sort(retained)));
    fprintf('Eliminated series  : %s\n', mat2str(sort(eliminated)));
    fprintf('Drop order         : %s\n', mat2str(dropped_order));
    fprintf('R²                 : %.4f\n', best_r2);
    fprintf('Estimated Hz       : %.3f%s\n', best_hz, flag_str(best_flag));

    % ================================================================
    % Package results
    % ================================================================
    results.best_subset    = sort(retained);
    results.predropped     = predropped;       % excluded before algorithm
    results.eliminated     = sort(eliminated); % dropped by algorithm
    results.dropped_order  = dropped_order;
    results.r2_trajectory  = r2_trajectory;
    results.hz_trajectory  = hz_trajectory;
    results.hz_flag        = hz_flag;
    results.fit_objects    = fit_trajectory;
    results.best_fit       = best_fit;
    results.best_gof       = best_gof;
    results.best_r2        = best_r2;
    results.best_hz        = best_hz;
    results.hz_range       = hz_range;

    % ================================================================
    % Plot
    % ================================================================
    plot_results(X, xvec, retained, predropped, eliminated, dropped_order, ...
                 r2_trajectory, hz_trajectory, hz_flag, ...
                 best_fit, best_r2, best_hz, hz_range);
end


% ================================================================
% LOCAL FUNCTIONS
% ================================================================

function [r2, hzapp, f, gof] = compute_fit(X, subset_idx, xvec)
    gM = mean(X(subset_idx, :), 1);
    [f, gof] = fit(xvec', gM', 'fourier1');
    hzapp = xvec(end) / (2*pi / f.w);
    r2    = gof.rsquare;
end


function out = is_out_of_range(hz, hz_range)
    out = hz < hz_range(1) || hz > hz_range(2);
end


function s = flag_str(is_flagged)
    if is_flagged
        s = '  *** OUT OF RANGE ***';
    else
        s = '';
    end
end


function plot_results(X, xvec, retained, predropped, eliminated, dropped_order, ...
                      r2_traj, hz_traj, hz_flag, ...
                      best_fit, best_r2, best_hz, hz_range)

    n_steps = numel(dropped_order);
    steps   = 0:n_steps;
    figure(8);clf;
    set(gcf,'Position', [100 100 1300 850]);

    % --- Panel 1: R² trajectory ---
    subplot(2,2,1);
    plot(steps, r2_traj, 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    hold on;
    flagged_steps = steps(hz_flag);
    if ~isempty(flagged_steps)
        plot(flagged_steps, r2_traj(hz_flag), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    [max_r2, max_idx] = max(r2_traj);
    plot(steps(max_idx), max_r2, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
    xlabel('Series dropped'); ylabel('R²');
    title('R² trajectory');
    legend('R²', 'Hz out of range', 'Best subset', 'Location', 'best');
    xticks(steps);
    grid on;

    % --- Panel 2: Hz trajectory ---
    subplot(2,2,2);
    plot(steps, hz_traj, 's-', 'LineWidth', 2, 'MarkerFaceColor', [0.4 0.7 0.4]);
    hold on;
    yline(hz_range(1), 'r--', 'LineWidth', 1.5);
    yline(hz_range(2), 'r--', 'LineWidth', 1.5);
    patch([steps(1) steps(end) steps(end) steps(1)], ...
          [hz_range(1) hz_range(1) hz_range(2) hz_range(2)], ...
          [0.8 1 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    if ~isempty(flagged_steps)
        plot(flagged_steps, hz_traj(hz_flag), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    xlabel('Series dropped'); ylabel('Hz');
    title('Estimated Hz trajectory');
    legend('Hz estimate', sprintf('Range [%.1f %.1f]', hz_range(1), hz_range(2)), ...
           'Location', 'best');
    xticks(steps);
    grid on;

    % --- Panel 3: Individual series + best fit ---
    % Three categories: predropped (grey), eliminated (light red), retained (light blue)
    subplot(2,2,3);
    for i = predropped
        plot(xvec, X(i,:), 'Color', [0.82 0.82 0.82], 'LineWidth', 0.5); hold on;
    end
    for i = eliminated
        plot(xvec, X(i,:), 'Color', [1 0.75 0.75], 'LineWidth', 0.5); hold on;
    end
    for i = retained
        plot(xvec, X(i,:), 'Color', [0.75 0.85 1], 'LineWidth', 0.5); hold on;
    end
    gM_best = mean(X(retained,:), 1);
    plot(xvec, gM_best, 'b', 'LineWidth', 2.5);
    x_fine = linspace(xvec(1), xvec(end), 500);
    plot(x_fine, best_fit(x_fine), 'r-', 'LineWidth', 2);
    xlabel('xvec'); ylabel('Amplitude');
    title('Series: grey=pre-dropped, red=eliminated, blue=retained');
    legend('Pre-dropped', 'Eliminated', 'Retained', 'Group mean', 'fourier1 fit', ...
           'Location', 'best');
    grid on;

    % --- Panel 4: Drop order ---
    subplot(2,2,4);
    bar(1:n_steps, dropped_order, 'FaceColor', [0.85 0.4 0.4]);
    xlabel('Elimination step'); ylabel('Series index dropped');
    title('Algorithm drop order');
    xticks(1:n_steps);
    grid on;

    sgtitle(sprintf(['Stepwise Fourier Selection | Best R² = %.4f | Hz = %.3f | ' ...
        '%d retained | %d pre-dropped | %d eliminated'], ...
        best_r2, best_hz, numel(retained), numel(predropped), numel(eliminated)));
end