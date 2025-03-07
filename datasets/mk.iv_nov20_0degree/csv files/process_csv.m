% data processing function to combine column csv files into a single csv
% file
function combined_csv = process_csv()
    % find csv files in current directory
    files = dir('*.csv');
    % create empty cell to hold csv data
    combined_csv = {};
    % iterate over all csv files to make one data cell
    for i = 1:length(files)
        myfile = sprintf('col%u.csv', i);
        temp_cell = readcell(myfile);
        newrows = size(temp_cell,1); 
        newcols = size(temp_cell,2);
        combined_csv(1:newrows, end+1:end+newcols) = temp_cell;
    end
    % create an array of headers
    headers = [];
    for i = 1:length(files)
        time_header = sprintf('col%u: Time (s)', i);
        header_temp = {time_header, "Bx-5mT (mT)", "By-5mT (mT)", "Bz-5mT (mT)", "Bx-130mT (mT)", "By-130mT (mT)", "Bz-130mT (mT)"};
        headers = [headers header_temp];
    end
    headers = num2cell(headers);
    headers = vertcat(headers{:});
    % add headers to the combined csv file
    combined_csv(1, :) = headers;
    % write file to the current directory
    writecell(combined_csv, 'mapping_data.csv')
end