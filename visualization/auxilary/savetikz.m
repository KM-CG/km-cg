function savetikz( figure_handles, prefix, seed, cleanFigure )
    if nargin < 4, cleanFigure = false; end
    for fig = 1:numel(figure_handles)
        fh = figure_handles{fig};
        set(0, 'currentfigure', fh);
        h = get(gca,'Title');
        fig_title = get(h, 'String');
        if isempty(fig_title) || strcmp(fig_title, '')
            error('Figure number %i has no title!', fh.Number);
        end
        title('');
        if cleanFigure
            cleanfigure('minimumPointsDistance',1e-8);
        end
        file_name = sprintf('%s__title_%s__seed_%i', prefix, fig_title, seed);
        file_name(~ismember(file_name,['A':'Z' 'a':'z' '0':'9'])) = '_';
        matlab2tikz(['../../tikz/' file_name '.tikz'], ...
            'figurehandle', fh, ...
            'width','\figwidth','height','\figheight',...
            'parseStrings', false, 'relativeDataPath', './tikz/', 'extraAxisOptions', {'mystyle'});%, 'ylabel style={yshift=-2ex}'});
    end
end

