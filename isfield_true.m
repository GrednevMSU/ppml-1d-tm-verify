function t = isfield_true(s, fld)
t = isfield(s, fld) && ~isempty(s.(fld)) && s.(fld);
end
