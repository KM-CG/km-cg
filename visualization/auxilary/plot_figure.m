function plot_figure(figure_handle, x_axis, values, v_min, v_max, linestyle, color, linewidth, x_clip_value, y_clip_value, crashes)
% Auxiliary plot function

        set(0, 'currentfigure', figure_handle); % open the figure that we consider

        too_large = find(x_axis > x_clip_value); % find x values that are above the clip value
        if ~isempty(too_large)
            x_clip_value = x_axis(too_large(1)); % take the next element so we do not cut off the line
        end
        idx = v_max > y_clip_value; % find elements above the y clip value
        idx = idx & (x_axis <= x_clip_value); % remove clipped indices
        plot(x_axis(idx), values(idx), linestyle, 'Color', color, 'LineWidth', linewidth); % plot mean
        % plot min and max
        fill([x_axis(idx), fliplr(x_axis(idx))], [v_max(idx), fliplr(v_min(idx))], color, 'FaceAlpha', 0.1, 'LineStyle', linestyle, 'EdgeColor', color);
        
        % get crashes that are not clipped
        crashes = crashes(crashes <= length(idx));
        
        % plot crashes
        plot(x_axis(crashes), values(crashes), 'x', 'Color', color);
end

