function s = engine_string()
if exist('OCTAVE_VERSION','builtin') ~= 0
    s = ['Octave ' OCTAVE_VERSION];
else
    v = ver('MATLAB'); s = ['MATLAB ' v.Version];
end
end
