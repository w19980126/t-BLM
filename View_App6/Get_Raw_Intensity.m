switch names{end}
    case strcmp(names{end},'jdoe') == 0
        names(end) = []
    otherwise
        return
end
strcmp(names,'jdoe')