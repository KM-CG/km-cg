for a = 1 : numel(algos_in_legend)
    b = algos_in_legend(a);
    plot(-1, -1, lines{b}, 'Color', colors(b, :), 'LineWidth', line_widths(b)); 
end
l = legend(algo_names(algos_in_legend));
set(l, 'Position', [0, 0, 1, 1]);
xlim([0 1]);
ylim([0 1]);
axis off;
