classdef StructuredMatrixFigureFormater < DenseMatrixFigureFormater
    methods
        function format_tikz_xaxis_steps(obj, show_x_axis, show_y_axis, ax1, N)            
            fig = gcf;
            if fig.Number == 9 || fig.Number == 6      
                obj.format_cg_steps_x_axis()
                obj.format_y_axis(show_y_axis);
                set(gca, 'color', 'none'); % make axis transparent
            end
            % use original axis to draw
            fig = gcf; fig.CurrentAxes = ax1;
            obj.set_common_x_axis_props();
            xticklabels(cell(size(xticks())));
            obj.format_y_axis(show_y_axis);
        end

        function format_y_axis(obj, show_y_axis)
            set(gca, 'yscale', 'log');
            ylim([10^(-3) 10^3]);
            yticks([10^-2, 0.1, 1, 10, 100]);
            if ~show_y_axis
                yticklabels(cell(size(yticks())));
            else
                yticklabels({'$10^{-2}$', '$10^{-1}$', '$10^0$', '$10^1$', '$10^2$'});
            end
        end
        
        function format_tikz_xaxis_time(obj, show_x_axis, show_y_axis, ax1, N)
            fig = gcf;
            show_x_axis = fig.Number == 9 || fig.Number == 6;
            if ~show_x_axis, xticks([]); end
            obj.format_y_axis(show_y_axis);
        end
    end
end

