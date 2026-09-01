function Write_Data_whole_bin(Data, output_file)
%WRITE_DATA_WHOLE_BIN Write a MATLAB matrix in the HR-GGM binary format.
%
%   Write_Data_whole_bin(Data)
%   Write_Data_whole_bin(Data, output_file)
%
%   Input matrix orientation:
%       rows    = local features
%       columns = samples
%
%   Data must be a real MATLAB double matrix. If output_file is omitted,
%   Data_whole.bin is written to the current folder.

    if nargin < 2 || isempty(output_file)
        output_file = fullfile(pwd, 'Data_whole.bin');
    end

    if ~isa(Data, 'double') || ~ismatrix(Data) || ~isreal(Data)
        error('Data must be a real MATLAB double matrix.');
    end

    if isempty(Data)
        error('Data must not be empty.');
    end

    if any(~isfinite(Data(:)))
        error('Data must contain only finite values.');
    end

    output_file = char(output_file);
    output_folder = fileparts(output_file);

    if ~isempty(output_folder) && ~isfolder(output_folder)
        error('Output folder does not exist: %s', output_folder);
    end

    fid = fopen(output_file, 'wb');

    if fid == -1
        error('Cannot open binary file: %s', output_file);
    end

    file_cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fwrite(fid, uint64(size(Data, 1)), 'uint64');
    fwrite(fid, uint64(size(Data, 2)), 'uint64');

    % MATLAB is column-major; transposition writes Data in row-major order.
    fwrite(fid, Data.', 'double');
end
