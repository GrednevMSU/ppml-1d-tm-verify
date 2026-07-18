function arr = cell2mat_struct(c)
% Concatenate a cell array of scalar structs (identical fields) into a struct array.
arr = c{1};
for i = 2:numel(c), arr(i) = c{i}; end
end
