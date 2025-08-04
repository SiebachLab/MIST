% Main script for formula output
% Organizes elements in decreasing order of abundance
% Created 2-9-2022 by Ellie Moreland

function [formula_out, total_out] = Formula_Output(values, elements)
        
    input = [elements; values];

    col_zeros = ~(cell2mat(input(2,:)) > 0.01);

    if(all(col_zeros))
        formula_out = '';
        total_out = '';
        return;
    end

    input(:,col_zeros) = [];
    
    [~, idx] = sort([input{2,:}], 'descend');
    input_sorted = input(:, idx);
    formula = convertCharsToStrings(sprintf('%s_%.2f', input_sorted{:}));
    formula_out = strcat('(', formula, ')');

    total = num2str(round(sum(cell2mat(input(2,:))),2));
    total_out = strcat('∑={', total, '}');
    
end