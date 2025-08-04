% Function for formatting input file by adding any left-out columns and
% filling any empty cells in the input with 0

function [T] = Input_Formatting(T)

    % List of column names requred by the code
    req_column_names = {'Point', 'Input_SiO2', 'Input_TiO2', 'Input_Al2O3', 'Input_Cr2O3', 'Input_FeO', 'Input_NiO', 'Input_MnO', 'Input_MgO', 'Input_CaO', 'Input_Na2O', 'Input_K2O',...
        'Input_P2O5', 'Input_SO3', 'Input_F', 'Input_Cl', 'Input_V2O3', 'Input_ZnO', 'Input_CoO', 'Input_BaO', 'Input_SrO', 'Input_B2O3', ...
        'Input_PbO', 'Input_CuO', 'Input_Sb2O3', 'Input_As2O5', 'Input_ThO2', 'Input_ZrO2', 'Input_HfO2', 'Input_Ag2O', 'Input_Y2O3', 'Input_La2O3', 'Input_Ce2O3', 'Input_Nd2O3'};
    
    % Check if each required column is present in the input file
    for i = 1:length(req_column_names)
        column_name = req_column_names{i};
    
        if ~ismember(column_name, T.Properties.VariableNames)
            % If the column name is not present in the input file, add that
            % column and fill it with zeros
            T.(column_name) = zeros(height(T), 1);
    
        else
            % If the column name is present in the input file, replace blank spaces with zeros in the existing column
            col_data = T.(column_name);
    
            if(iscell(col_data) && ~strcmp(column_name, 'Point'))
                % If the column happens to be a cell array, convert it to numeric and replace blanks
                col_data = str2double(col_data);
                
            end
         
            if isnumeric(col_data)
                col_data(isnan(col_data)) = 0; % Replace NaNs resulting from blank spaces
            end
            
            
            T.(column_name) = col_data;
            
        end
    end

end