function results = stepwise_fourier_selection_composite(X, xvec, min_series, hz_range, predropped)
% STEPWISE_FOURIER_SELECTION_COMPOSITE
% Finds the best subset of time-series optimising R² within hz_range.
%
% SELECTION CRITERION (applied at every step):
%   Primary   : Hz must stay within hz_range (hard gate)
%   Secondary : Among in-range candidates, pick highest R²
%   Deadlock  : If ALL candidate drops move Hz out of range,
%               fall back to least-bad Hz distance (flagged in output)
%
% INPUTS:
%   X            - [40 x T] matrix of time series (rows = series, cols = time)
%   xvec         - time vector, e.g. 1:10
%   min_series   - minimum series to retain (default: 30)
%   hz_range     - [lo hi] acceptable Hz window (default: [1.5 2.5])
%   predropped   - indices to exclude before selection begins (default: [])
%
% EXAMPLE:
%   results = stepwise_fourier_selection_composite(X, 1:10, 30, [1.5 2.5], [1 2]);

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

    % Warn if baseline Hz already out of range
    [r2_base, hz_base] = compute_fit(X, available, xvec);
    if is_out_of_range(hz_base, hz_range)
        warning('Baseline Hz (%.3f) is already outside hz_range [%.1f %.1f]. Deadlock fallback may be used from step 1.', ...
            hz_base, hz_range(1), hz_range(2));
    end

    fprintf('Available series: %d | Pre-dropped: %d | Steps to run: %d\n\n', ...
        n_avail, numel(predropped), n_steps);

    % --- Initialise ---
    retained      = available;
    eliminated    = [];
    dropped_order = zeros(1, n_steps);
    deadlock_flag = false(1, n_steps + 1);  % true = fallback used at this step

    r2_trajectory  = zeros(1, n_steps + 1);
    hz_trajectory  = zeros(1, n_steps + 1);
    hzd_trajectory = zeros(1, n_steps + 1);
    hz_flag        = false(1, n_steps + 1);
    fit_trajectory = cell(1,  n_steps + 1);

    % Baseline
    r2_trajectory(1)   = r2_base;
    hz_trajectory(1)   = hz_base;
    hzd_trajectory(1)  = hz_distance(hz_base, hz_range);
    hz_flag(1)         = hzd_trajectory(1) > 0;
    [~, ~, fit_trajectory{1}] = compute_fit(X, retained, xvec);

    fprintf('Step 0 (%d available series): R² = %.4f | Hz = %.3f (dist=%.4f)%s\n', ...
        n_avail, r2_trajectory(1), hz_trajectory(1), hzd_trajectory(1), ...
        flag_str(hz_flag(1)));

    % ================================================================
    % PHASE 1: Backward elimination
    %   Primary:  keep Hz in range (hard gate)
    %   Secondary: maximise R² among in-range candidates
    %   Deadlock:  fall back to min Hz distance if no in-range candidate
    % ================================================================
    for step = 1:n_steps

        n_retained     = numel(retained);
        r2_candidates  = nan(1, n_series);
        hz_candidates  = nan(1, n_series);
        hzd_candidates = nan(1, n_series);

        % Evaluate every possible drop
        for k = 1:n_retained
            candidate = retained([1:k-1, k+1:n_retained]);
            [r2_candidates(retained(k)), hz_candidates(retained(k))] = ...
                compute_fit(X, candidate, xvec);
            hzd_candidates(retained(k)) = hz_distance(hz_candidates(retained(k)), hz_range);
        end

        % Identify in-range candidates (Hz distance == 0)
        in_range_mask = (hzd_candidates == 0);

        if any(in_range_mask)
            % Normal path: pick highest R² among in-range drops
            r2_masked = r2_candidates;
            r2_masked(~in_range_mask) = -inf;
            [best_r2_step, best_drop_idx] = max(r2_masked);
            deadlock_flag(step+1) = false;
            mode_str = 'normal';
        else
            % Deadlock: no drop keeps Hz in range — pick least-bad Hz distance
            [~, best_drop_idx] = min(hzd_candidates);
            best_r2_step = r2_candidates(best_drop_idx);
            deadlock_flag(step+1) = true;
            mode_str = 'DEADLOCK FALLBACK';
        end

        best_hz_step  = hz_candidates(best_drop_idx);
        best_hzd_step = hzd_candidates(best_drop_idx);

        % Update lists
        retained(retained == best_drop_idx) = [];
        eliminated(end+1)   = best_drop_idx;        %#ok<AGROW>
        dropped_order(step) = best_drop_idx;

        r2_trajectory(step+1)  = best_r2_step;
        hz_trajectory(step+1)  = best_hz_step;
        hzd_trajectory(step+1) = best_hzd_step;
        hz_flag(step+1)        = best_hzd_step > 0;
        [~, ~, fit_trajectory{step+1}] = compute_fit(X, retained, xvec);

        fprintf('Step %d: Dropped %2d | R² = %.4f | Hz = %.3f (dist=%.4f)%s | Retained: %d [%s]\n', ...
            step, best_drop_idx, best_r2_step, best_hz_step, best_hzd_step, ...
            flag_str(hz_flag(step+1)), numel(retained), mode_str);
    end

    % ================================================================
    % PHASE 2: Swap refinement
    %   Same composite logic: Hz in range is gate, then best R²
    %   Deadlock swaps: only accept if Hz distance improves
    % ================================================================
    fprintf('\n--- Swap refinement (composite: Hz gate + R² criterion) ---\n');
    improved   = true;
    swap_count = 0;

    while improved
        improved = false;

        for di = 1:numel(eliminated)
            for ri = 1:numel(retained)

                candidate = retained;
                candidate(ri) = eliminated(di);

                [r2_swap, hz_swap] = compute_fit(X, candidate, xvec);
                [r2_curr, hz_curr] = compute_fit(X, retained,  xvec);

                hzd_swap = hz_distance(hz_swap, hz_range);
                hzd_curr = hz_distance(hz_curr, hz_range);

                % Determine if swap is an improvement under composite logic
                swap_accepted = false;
                swap_mode     = '';

                if hzd_swap == 0 && hzd_curr == 0
                    % Both in range: accept if R² improves
                    if r2_swap > r2_curr
                        swap_accepted = true;
                        swap_mode     = 'both in range, R² improved';
                    end
                elseif hzd_swap == 0 && hzd_curr > 0
                    % Swap brings Hz into range: always accept
                    swap_accepted = true;
                    swap_mode     = 'Hz restored to range';
                elseif hzd_swap > 0 && hzd_curr > 0
                    % Both out of range (deadlock territory): accept if Hz distance improves
                    if hzd_swap < hzd_curr
                        swap_accepted = true;
                        swap_mode     = 'deadlock: Hz distance improved';
                    end
                end
                % hzd_swap > 0 && hzd_curr == 0: never accept (would leave range)

                if swap_accepted
                    fprintf('  Swapped IN %2d, OUT %2d | R² %.4f->%.4f | Hz %.3f->%.3f (dist %.4f->%.4f) [%s]\n', ...
                        eliminated(di), retained(ri), r2_curr, r2_swap, ...
                        hz_curr, hz_swap, hzd_curr, hzd_swap, swap_mode);

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
        fprintf('  No swaps improved fit.\n');
    else
        fprintf('  %d swap(s) accepted.\n', swap_count);
    end

    % ================================================================
    % Final fit
    % ================================================================
    [best_r2, best_hz, best_fit, best_gof] = compute_fit(X, retained, xvec);
    best_hzd  = hz_distance(best_hz, hz_range);
    best_flag = best_hzd > 0;
    n_deadlocks = sum(deadlock_flag);

    fprintf('\n--- Final Result ---\n');
    fprintf('Pre-dropped series : %s\n',  mat2str(predropped));
    fprintf('Retained series    : %s\n',  mat2str(sort(retained)));
    fprintf('Eliminated series  : %s\n',  mat2str(sort(eliminated)));
    fprintf('Drop order         : %s\n',  mat2str(dropped_order));
    fprintf('R²                 : %.4f\n', best_r2);
    fprintf('Estimated Hz       : %.3f (dist=%.4f)%s\n', best_hz, best_hzd, flag_str(best_flag));
    if n_deadlocks > 0
        fprintf('Deadlock steps     : %d (fallback to min Hz distance)\n', n_deadlocks);
    end

    % ================================================================
    % Package results
    % ================================================================
    results.best_subset    = sort(retained);
    results.predropped     = predropped;
    results.eliminated     = sort(eliminated);
    results.dropped_order  = dropped_order;
    results.deadlock_flag  = deadlock_flag;   % which steps used fallback
    results.r2_trajectory  = r2_trajectory;
    results.hz_trajectory  = hz_trajectory;
    results.hzd_trajectory = hzd_trajectory;
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
                 r2_trajectory, hz_trajectory, hzd_trajectory, ...
                 hz_flag, deadlock_flag, best_fit, best_r2, best_hz, hz_range);
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
    if is_flagged, s = '  *** OUT OF RANGE ***'; else, s = ''; end
end


function plot_results(X, xvec, retained, predropped, eliminated, dropped_order, ...
                      r2_traj, hz_traj, hzd_traj, hz_flag, deadlock_flag, ...
                      best_fit, best_r2, best_hz, hz_range)

    n_steps = numel(dropped_order);
    steps   = 0:n_steps;

    figure('Position', [100 100 1300 850]);

    % --- Panel 1: R² trajectory — primary optimisation output ---
    subplot(2,2,1);
    plot(steps, r2_traj, 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    hold on;
    % Mark deadlock steps in orange
    dl_steps = steps(deadlock_flag);
    if ~isempty(dl_steps)
        plot(dl_steps, r2_traj(deadlock_flag), 'o', 'MarkerSize', 10, ...
            'LineWidth', 2, 'Color', [0.9 0.5 0], 'MarkerFaceColor', [0.9 0.5 0]);
    end
    % Mark out-of-range steps in red
    flagged_steps = steps(hz_flag);
    if ~isempty(flagged_steps)
        plot(flagged_steps, r2_traj(hz_flag), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    [max_r2, max_idx] = max(r2_traj);
    plot(steps(max_idx), max_r2, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
    xlabel('Series dropped'); ylabel('R²');
    title('R² trajectory (primary criterion)');
    legend('R²', 'Deadlock step', 'Hz out of range', 'Best R²', 'Location', 'best');
    xticks(steps); grid on;

    % --- Panel 2: Hz trajectory with range shading ---
    subplot(2,2,2);
    plot(steps, hz_traj, 's-', 'LineWidth', 2, 'MarkerFaceColor', [0.4 0.7 0.4]);
    hold on;
    yline(hz_range(1), 'r--', 'LineWidth', 1.5);
    yline(hz_range(2), 'r--', 'LineWidth', 1.5);
    patch([steps(1) steps(end) steps(end) steps(1)], ...
          [hz_range(1) hz_range(1) hz_range(2) hz_range(2)], ...
          [0.8 1 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    if ~isempty(dl_steps)
        plot(dl_steps, hz_traj(deadlock_flag), 'o', 'MarkerSize', 10, ...
            'LineWidth', 2, 'Color', [0.9 0.5 0], 'MarkerFaceColor', [0.9 0.5 0]);
    end
    if ~isempty(flagged_steps)
        plot(flagged_steps, hz_traj(hz_flag), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    xlabel('Series dropped'); ylabel('Hz');
    title('Hz trajectory (hard gate)');
    legend('Hz', sprintf('Range [%.1f %.1f]', hz_range(1), hz_range(2)), '', '', ...
           'Deadlock', 'Out of range', 'Location', 'best');
    xticks(steps); grid on;

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

    % --- Panel 4: Drop order, coloured by deadlock ---
    subplot(2,2,4);
    bar_colours = repmat([0.3 0.6 1], n_steps, 1);  % default blue
    for s = 1:n_steps
        if deadlock_flag(s+1)
            bar_colours(s,:) = [0.9 0.5 0];  % orange = deadlock step
        end
    end
    b = bar(1:n_steps, dropped_order, 'FaceColor', 'flat');
    b.CData = bar_colours;
    xlabel('Elimination step'); ylabel('Series index dropped');
    title('Drop order (orange = deadlock fallback step)');
    xticks(1:n_steps); grid on;

    n_deadlocks = sum(deadlock_flag);
    dl_note = '';
    if n_deadlocks > 0
        dl_note = sprintf(' | %d deadlock step(s)', n_deadlocks);
    end

    sgtitle(sprintf(['Composite Selection (R² within Hz range) | R² = %.4f | Hz = %.3f | ' ...
        '%d retained | %d pre-dropped | %d eliminated%s'], ...
        best_r2, best_hz, numel(retained), numel(predropped), numel(eliminated), dl_note));
end