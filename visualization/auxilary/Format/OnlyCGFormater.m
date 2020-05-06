classdef OnlyCGFormater < DenseMatrixFigureFormater
    % formater used when using only CG algorithms (i.e. no inducing input
    % methods)
    methods                
        function format_inducing_input_x_axis(obj, N)
            axl = gca;
            set(axl,'xaxisLocation','top');
            obj.set_common_x_axis_props();
            xticklabels(axl, cell(size(xticks())));
        end
    end
end

