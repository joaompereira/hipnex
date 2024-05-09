%% Post-process results of numerical_experiments

clc
clearvars

% Load results
load results/cubic_min_max_experiment.mat

it_times = cellfun(@(x) x.it_times, stats, 'UniformOutput', false);
total_times = cellfun(@(x) x.total_time, stats);
norm_vals = cellfun(@(x) x.F_norms, stats, 'UniformOutput', false);

styles = ["-o", "--<", "-diamond", ":>", '-.square'];

xlims = [.5, 10, 15];

legends = cell(nalgs, 1);

close(figure(1))
hf = figure(1); clf
hf.Position = [100 100 900 300];


for ind_n = 1:nn

    %% Plot
    subplot('Position', [.07 + .33 * (ind_n - 1), .27, .25, .62])


    for i=1:nalgs
        hp = plot(cumsum(it_times{i, ind_n}), norm_vals{i, ind_n}, styles(i), MarkerSize=5);
        if i==1
            hold on
        end
        set(hp, 'MarkerEdgeColor', 'k', 'markerfacecolor', get(hp, 'color'));
    end
    hold off
    set(gca, 'YScale', 'log', 'FontSize', 10, 'TickLabelInterpreter', "latex")
    
    xlim([0, mean(cellfun(@sum, it_times(1:2, ind_n))) * 3])
    ylim([hyperparameters.tol / 2, stats{1, ind_n}.F_norms(1) * 2])
    
    title(sprintf('$n=%d$', nvals(ind_n)), Interpreter="latex")
    xlabel('time (s)', Interpreter="latex")
    ylabel('$\|F(x_k)\|$', Interpreter="latex")
    
end

hl = legend(Algs(:, 1), Location="northeast",Interpreter="latex",Orientation='horizontal');
hl.Position = [0.2    0.04    0.6    0.05];
pdf_path = sprintf("results/cubic_plot");

set(gcf,'Units','Inches');
pos = get(gcf ,'Position');
set(gcf ,'PaperPositionMode','Auto',...
            'PaperUnits','Inches',...
            'PaperSize',[pos(3), pos(4)]);
print(gcf, pdf_path,'-dpdf','-r0');
