%READING_DATA_FROM_EXCEL.m

% --- IMPORTANT: Make sure the Excel file is in the same folder as this script! ---

raw_data_filename = 'Raw_Project_Data.xlsx';

try
    % Use xlsread to load data from the first sheet, first column ('A:A').
    % The output 'Data' will be a column vector.
    [Data, ~, ~] = xlsread(raw_data_filename, 1, 'A:A'); 
    
    % Remove any potential NaN values (often added by xlsread if the range is too big)
    Data(isnan(Data)) = [];
    
    % Confirm successful load
    fprintf('Data successfully loaded from %s.\n', raw_data_filename);
    fprintf('The loaded vector "Data" contains N = %d observations.\n', length(Data));
    
catch ME
    fprintf('\nERROR: Could not load data from Excel using xlsread.\n');
    fprintf('Please check the file name and ensure the file is closed.\n');
    fprintf('MATLAB Error: %s\n', ME.message);
    
end
% Now the variable 'Data' is available for the rest of the project code.