function savelegend( fig_handle, legend_title)
        set(0, 'currentfigure', fig_handle);
        title(''); % remove title
        file_name = sprintf('legend__%s', legend_title);
        file_name(~ismember(file_name,['A':'Z' 'a':'z' '0':'9'])) = '_';
        matlab2tikz(['../../tikz/' file_name '.tikz'], ...
       'width','\figwidth','height','\figheight',...
       'parseStrings', false, 'relativeDataPath', './tikz/', 'extraAxisOptions', {'mystyle', sprintf('legend to name=%s',legend_title)});%, 'ylabel style={yshift=-2ex}'});
end

