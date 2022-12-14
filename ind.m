function index = ind(arr, val)
% ind - Find the index of an array element of a given value
%
% index = ind(arr, val)
%
% ind searches for the value val in the array arr
% and returns the index of the closest matching array element

if isempty(arr)
    index = [];
else
    if isempty(val)
        error('val cannot be empty.')
    end
    if ~isvector(val)
        error('val must be a scalar or a vector')
    end
    index = zeros(1, numel(val));
    for k = 1:numel(val);
        dif = abs(arr-val(k));
        [~, index(k)] = min(dif(:));
    end
end