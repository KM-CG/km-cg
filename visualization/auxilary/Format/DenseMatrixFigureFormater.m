classdef DenseMatrixFigureFormater
    % auxiliary class to format figures
    methods
        function format_figure(obj, TIKZ, show_time, show_x_axis, show_y_axis, data_set, N, D)    
            axis tight;
            set(gca, 'yscale', 'log');
            ax1 = gca;
            if TIKZ
                text(1, 1, ['\datasetCaption{' data_set '}\datasetSize{' sprintf('%i}', N) '\datasetDim{' sprintf('%i}', D)], 'units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
                if ~show_time
                    obj.format_tikz_xaxis_steps(show_x_axis, show_y_axis, ax1, N);
                else
                    obj.format_tikz_xaxis_time(show_x_axis, show_y_axis, ax1, N);
                end
            end
        end
        
        function format_tikz_xaxis_steps(obj, show_x_axis, show_y_axis, ax1, N)
            if show_x_axis
                obj.format_cg_steps_x_axis();
                obj.format_y_axis(show_y_axis); % we have to format the y-axis again
                set(gca, 'color', 'none'); % make axis transparent
            end
            % use original axis to draw
            fig = gcf; fig.CurrentAxes = ax1;
            hold on;
            obj.format_inducing_input_x_axis(N);
            obj.format_y_axis(show_y_axis);
        end
        
        function format_cg_steps_x_axis(obj)
            ax = axes('XAxisLocation', 'bottom');
            hold on;
            obj.set_common_x_axis_props();
            lbls = num2cell(xticks());
            xticklabels(ax, lbls);
        end
        
        function format_inducing_input_x_axis(obj, N)
            axl = gca;
            set(axl,'xaxisLocation','top');
            obj.set_common_x_axis_props();
            lbls = num2cell(ceil(sqrt(N*xticks())));
            xticklabels(axl, {'$M=$', lbls{2:end}});
        end

        function set_common_x_axis_props(obj)
            xlim([1 99]);
            xticks([20 40 60 80]);
        end

        function format_y_axis(obj, show_y_axis)
            set(gca, 'yscale', 'log');
            ylim([10^(-3) 10^2]);
            yticks([10^-2, 0.1, 1, 10]);
            if ~show_y_axis
                yticklabels(cell(size(yticks())));
            else
                yticklabels({'$10^{-2}$', '$10^{-1}$', '$10^0$', '$10^1$'});
            end
        end
        
        function format_tikz_xaxis_time(obj, show_x_axis, show_y_axis, ax1, N)
            set(ax1,'xaxisLocation','bottom');
            % x-limits are defined in show.m after reading in the results
            obj.format_y_axis(show_y_axis);
        end
    end
end

