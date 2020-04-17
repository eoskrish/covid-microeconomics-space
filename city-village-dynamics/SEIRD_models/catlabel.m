function text_output = catlabel(ind_a,ind_e)

    part1 = [];
    part2 = [];

    if ind_a == 1
        part1 = 'children';
    elseif ind_a == 2
        part1 = 'young';
    elseif ind_a == 3
        part1 = 'old';
    end

    if ind_e == 1
        part2 = 'immobile poor';
    elseif ind_e == 2
        part2 = 'mobile poor';
    elseif ind_e == 3
        part2 = 'rich';
    end

    if ind_a == 0 && ind_e == 0 
        text_output = [];
    end

    if ind_a == 0 && ind_e ~= 0 
        text_output = part2;
    end

    if ind_a ~= 0 && ind_e == 0 
        text_output = part1;
    end

    if ind_a ~= 0 && ind_e ~= 0 
        text_output = [part1 ', ' part2];
    end

end

