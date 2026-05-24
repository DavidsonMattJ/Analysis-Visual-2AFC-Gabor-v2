function results = stepwise_fourier_selection_hz(X, xvec, min_series, hz_range, predropped)
% STEPWISE_FOURIER_SELECTION_HZ
% Finds the best subset of time-series minimising distance from hz_range,
% using fourier1 group-level fit. R² is tracked as a secondary diagnostic.
%
% INPUTS:
%   X            - [40 x T] matrix of time series (rows = series, cols = time)
%   xvec         - time vector, e.g. 1:10
%   min_series   - minimum series to retain (default: 30)
%   hz_range     - [lo hi] target Hz window (default: [1.5 2.5])
%   predropped   - indices to exclude before selection begins (default: [])
%
% OPTIMISATION CRITERION:
%   At each step, drop the series whose removal minimises Hz distance from
%   hz_range. Hz distance = 0 if inside range, else distance to nearest boundary.
%   R² is computed and stored at every step but does not influence selection.
%
% EXAMPLE:
%   results = stepwise_fourier_selection_hz(X, 1:10, 30, [1.5 2.5], [1 2]);

    % --- Defaults ---
    if nargin < 3, min_series = 30;        end
    if nargin < 4, hz_range   = [1.5 2.5]; end
    if nargin < 5, predropped = [];        end

    [n_series, ~] = size(X);

    % Validate predropped
    if ~isempty(predropped)
        predropped = unique(predropped(:)');
        assert(all(predropped >= 1 & predropped <= n_series), ...
            'predropped contains indices outside valid range [1, %d]', n_series);
        fprintf('Pre-dropped series (excluded from selection): %s\n', mat2str(predropped));
    end

    available = setdiff(1:n_series, predropped);
    n_avail   = numel(available);
    n_steps   = n_avail - min_series;

    assert(n_avail >= min_series, ...
        'Only %d series available after pre-dropping; cannot retain %d.', ...
        n_avail, min_series);

    fprintf('Available series: %d | Pre-dropped: %d | Steps to run: %d\n\n', ...
        n_avail, numel(predropped), n_steps);

    % --- Initialise ---
    retained      = available;
    eliminated    = [];
    dropped_order = zeros(1, n_steps);

    r2_trajectory   = zeros(1, n_steps + 1);
    hz_trajectory   = zeros(1, n_steps + 1);
    hzd_trajectory  = zeros(1, n_steps + 1);  % Hz distance from range
    hz_flag         = false(1, n_steps + 1);
    fit_trajectory  = cell(1,  n_steps + 1);

    % Baseline
    [r2_trajectory(1), hz_trajectory(1), fit_trajectory{1}] = ...
        compute_fit(X, retained, xvec);
    hzd_trajectory(1) = hz_distance(hz_trajectory(1), hz_range);
    hz_flag(1)        = hzd_trajectory(1) > 0;

    fprintf('Step 0 (%d available series): Hz = %.3f (dist=%.4f) | R² = %.4f%s\n', ...
        n_avail, hz_trajectory(1), hzd_trajectory(1), r2_trajectory(1), ...
        flag_str(hz_flag(1)));

    % ================================================================
    % PHASE 1: Backward elimination — criterion: minimise Hz distance
    % ================================================================
    for step = 1:n_steps

        n_retained     = numel(retained);
        hzd_candidates = nan(1, n_series);
        hz_candidates  = nan(1, n_series);

        % Try removing each retained series
        for k = 1:n_retained
            candidate = retained([1:k-1, k+1:n_retained]);
            [~, hz_k] = compute_fit(X, candidate, xvec);
            hzd_candidates(retained(k)) = hz_distance(hz_k, hz_range);
            hz_candidates(retained(k))  = hz_k;
        end

        % Best drop = minimum Hz distance from range
        [best_hzd_step, best_drop_idx] = min(hzd_candidates);
        best_hz_step                   = hz_candidates(best_drop_idx);

        % Update lists
        retained(retained == best_drop_idx) = [];
        eliminated(end+1)   = best_drop_idx;        %#ok<AGROW>
        dropped_order(step) = best_drop_idx;

        % Compute R² for this retained set (diagnostic only)
        [r2_step, ~, fit_step] = compute_fit(X, retained, xvec);

        r2_trajectory(step+1)  = r2_step;
        hz_trajectory(step+1)  = best_hz_step;
        hzd_trajectory(step+1) = best_hzd_step;
        hz_flag(step+1)        = best_hzd_step > 0;
        fit_trajectory{step+1} = fit_step;

        fprintf('Step %d: Dropped series %2d | Hz = %.3f (dist=%.4f)%s | R² = %.4f | Retained: %d\n', ...
            step, best_drop_idx, best_hz_step, best_hzd_step, ...
            flag_str(hz_flag(step+1)), r2_step, numel(retained));
    end

    % ================================================================
    % PHASE 2: Swap refinement — criterion: minimise Hz distance
    % ================================================================
    fprintf('\n--- Swap refinement (Hz criterion, algorithm-eliminated only) ---\n');
    improved   = true;
    swap_count = 0;

    while improved
        improved = false;

        for di = 1:numel(eliminated)
            for ri = 1:numel(retained)

                candidate = retained;
                candidate(ri) = eliminated(di);

                [~, hz_swap] = compute_fit(X, candidate, xvec);
                [~, hz_curr] = compute_fit(X, retained,  xvec);

                hzd_swap = hz_distance(hz_swap, hz_range);
                hzd_curr = hz_distance(hz_curr, hz_range);

                if hzd_swap < hzd_curr
                    [r2_swap, ~] = compute_fit(X, candidate, xvec);
                    fprintf('  Swapped IN %2d, OUT %2d | Hz dist %.4f -> %.4f | Hz %.3f -> %.3f | R² = %.4f\n', ...
                        eliminated(di), retained(ri), hzd_curr, hzd_swap, ...
                        hz_curr, hz_swap, r2_swap);

                    old_retained   = retained(ri);
                    retained(ri)   = eliminated(di);
                    eliminated(di) = old_retained;
                    swap_count     = swap_count + 1;
                    improved       = true;
                end
            end
        end
    end

    if swap_count == 0
        fprintf('  No swaps improved Hz distance.\n');
    else
        fprintf('  %d swap(s) accepted.\n', swap_count);
    end

    % ================================================================
    % Final fit
    % ================================================================
    [best_r2, best_hz, best_fit, best_gof] = compute_fit(X, retained, xvec);
    best_hzd  = hz_distance(best_hz, hz_range);
    best_flag = best_hzd > 0;

    fprintf('\n--- Final Result ---\n');
    fprintf('Pre-dropped series : %s\n', mat2str(predropped));
    fprintf('Retained series    : %s\n', mat2str(sort(retained)));
    fprintf('Eliminated series  : %s\n', mat2str(sort(eliminated)));
    fprintf('Drop order         : %s\n', mat2str(dropped_order));
    fprintf('Estimated Hz       : %.3f (dist=%.4f)%s\n', best_hz, best_hzd, flag_str(best_flag));
    fprintf('R² (diagnostic)    : %.4f\n', best_r2);

    % ================================================================
    % Package results
    % ================================================================
    results.best_subset    = sort(retained);
    results.predropped     = predropped;
    results.eliminated     = sort(eliminated);
    results.dropped_order  = dropped_order;
    results.r2_trajectory  = r2_trajectory;
    results.hz_trajectory  = hz_trajectory;
    results.hzd_trajectory = hzd_trajectory;  % key optimisation criterion
    results.hz_flag        = hz_flag;
    results.fit_objects    = fit_trajectory;
    results.best_fit       = best_fit;
    results.best_gof       = best_gof;
    results.best_r2        = best_r2;
    results.best_hz        = best_hz;
    results.best_hzd       = best_hzd;
    results.hz_range       = hz_range;

    % ================================================================
    % Plot
    % ================================================================
    plot_results(X, xvec, retained, predropped, eliminated, dropped_order, ...
                 r2_trajectory, hz_trajectory, hzd_trajectory, hz_flag, ...
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


function d = hz_distance(hz, hz_range)
% Distance from hz_range: 0 if inside, else distance to nearest boundary
    if hz < hz_range(1)
        d = hz_range(1) - hz;
    elseif hz > hz_range(2)
        d = hz - hz_range(2);
    else
        d = 0;
    end
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
                      r2_traj, hz_traj, hzd_traj, hz_flag, ...
                      best_fit, best_r2, best_hz, hz_range)

    n_steps = numel(dropped_order);
    steps   = 0:n_steps;

    figure(9); clf;
    set(gcf,'Position', [100 100 1300 850]);

    % --- Panel 1: Hz distance trajectory (primary criterion) ---
    subplot(2,2,1);
    plot(steps, hzd_traj, 'o-', 'LineWidth', 2, 'MarkerFaceColor', [0.4 0.7 0.4]);
    hold on;
    yline(0, 'k--', 'LineWidth', 1.5);
    % Mark steps still out of range
    flagged_steps = steps(hz_flag);
    if ~isempty(flagged_steps)
        plot(flagged_steps, hzd_traj(hz_flag), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    [min_hzd, min_idx] = min(hzd_traj);
    plot(steps(min_idx), min_hzd, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
    xlabel('Series dropped'); ylabel('Hz distance from range');
    title('Hz distance trajectory (primary criterion)');
    legend('Hz distance', 'In-range boundary', 'Out of range', 'Best subset', ...
           'Location', 'best');
    xticks(steps);
    grid on;

    % --- Panel 2: Hz and R² dual-axis ---
    subplot(2,2,2);
    yyaxis left
    plot(steps, hz_traj, 's-', 'LineWidth', 2, 'MarkerFaceColor', [0.4 0.7 0.4]);
    hold on;
    yline(hz_range(1), 'g--', 'LineWidth', 1.5);
    yline(hz_range(2), 'g--', 'LineWidth', 1.5);
    patch([steps(1) steps(end) steps(end) steps(1)], ...
          [hz_range(1) hz_range(1) hz_range(2) hz_range(2)], ...
          [0.8 1 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylabel('Hz');

    yyaxis right
    plot(steps, r2_traj, 'o--', 'LineWidth', 1.5, 'MarkerFaceColor', [0.4 0.4 0.9]);
    ylabel('R² (diagnostic)');

    xlabel('Series dropped');
    title('Hz trajectory + R² diagnostic');
    legend('Hz', sprintf('Range [%.1f %.1f]', hz_range(1), hz_range(2)), '', '', 'R²', ...
           'Location', 'best');
    xticks(steps);
    grid on;

    % --- Panel 3: Individual series + best fit ---
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
    bar(1:n_steps, dropped_order, 'FaceColor', [0.4 0.7 0.4]);
    xlabel('Elimination step'); ylabel('Series index dropped');
    title('Algorithm drop order (Hz-optimised)');
    xticks(1:n_steps);
    grid on;

    sgtitle(sprintf(['Hz-Optimised Selection | Hz = %.3f (dist=%.4f) | R² = %.4f | ' ...
        '%d retained | %d pre-dropped | %d eliminated'], ...
        best_hz, min(hzd_traj), best_r2, numel(retained), numel(predropped), numel(eliminated)));
end