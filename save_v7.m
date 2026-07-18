function save_v7(fn, fx) %#ok<INUSD>
% Cross-engine -v7 save (option position differs MATLAB vs Octave).
if exist('OCTAVE_VERSION','builtin') ~= 0
    save('-v7', fn, 'fx');
else
    save(fn, 'fx', '-v7');
end
end
