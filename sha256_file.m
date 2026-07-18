function h = sha256_file(fn)
% SHA-256 hex digest of a file's raw bytes.
% Cross-engine: Octave uses builtin hash(); MATLAB uses java.security.MessageDigest.
% Byte-faithful for binary .mat files (verified on high bytes 0x80..0xFF).

fid = fopen(fn, 'r');
if fid < 0, error('sha256_file:open', 'cannot open %s', fn); end
bytes = fread(fid, Inf, '*uint8');
fclose(fid);

use_octave_hash = (exist('hash', 'builtin') == 5) || (exist('OCTAVE_VERSION', 'builtin') ~= 0);

if use_octave_hash
    h = hash('sha256', char(bytes(:).'));   % Octave builtin, byte-faithful
else
    md = java.security.MessageDigest.getInstance('SHA-256');   % MATLAB path
    md.update(bytes);
    dig = typecast(md.digest(), 'uint8');    % java int8 -> uint8
    h = lower(reshape(dec2hex(dig, 2).', 1, []));
end
end
