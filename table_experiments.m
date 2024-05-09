%% Post-process results of numerical_experiments

clc
clearvars

% Load results
load results/cubic_min_max_experiment.mat

it_times = cellfun(@(x) x.it_times, stats, 'UniformOutput', false);
total_times = cellfun(@(x) x.total_time, stats);
norm_vals = cellfun(@(x) x.F_norms, stats, 'UniformOutput', false);

for ind_n = 1:nn

    %% Table
    fprintf('\n\n:: n=%d ::\n', nvals(ind_n));
    fprintf('+------------------+-------+-----+--------+------+-----+-----+-----+\n');
    fprintf('|                  |       |     |  final |linear|  F  |  DF |total|\n');
    fprintf('|      Method      |time(s)|iters| norm F |solves|evals|evals|inner|\n');
    fprintf('+------------------+-------+-----+--------+------+-----+-----+-----+\n')
    for i=1:nalgs
        fprintf('|%18s', Algs{i, 1});
        fprintf('| %6.4g', stats{i, ind_n}.total_time);
        fprintf('| %3d ', stats{i, ind_n}.n_iters);
        fprintf('|%.2e', stats{i, ind_n}.F_norms(end))
        fprintf('|  %3d ', stats{i, ind_n}.linsolve_count);
        fprintf('| %3d ', stats{i, ind_n}.F_evals);
        fprintf('| %3d ', stats{i, ind_n}.DF_evals);
        if isfield(stats{i, ind_n}, 'total_inner')
            fprintf('|%4d |\n', stats{i, ind_n}.total_inner);
        else
            fprintf('|     |\n');
        end 
        fprintf('+------------------+-------+-----+--------+------+-----+-----+-----+\n')
    end
    
end
