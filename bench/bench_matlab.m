% bench_matlab.m — wall-clock timing of the ORIGINAL RTA_1d_tm in MATLAB/Octave.
% Pair with bench/benchmark.py to fill the ratio column in BENCHMARK.md.
% Run from the repo root:  matlab -batch "addpath('bench'); bench_matlab"
addpath(fullfile(pwd,'matlab_src','1d_tm'));
addpath(fullfile(pwd,'matlab_src','general'));

epsAu = -4000 + 300i; epsbg = 10.05; k0 = 150*(2*pi/1240);
theta = 0.1;
mkargs = @(hp) {3.2, 3, 1.0, 1.0, ...
    [epsbg epsbg epsbg], [epsbg epsbg epsbg], ...
    [epsbg epsAu epsbg], [epsbg epsAu epsbg], ...
    [0 0 0 0], [0 0.73 0], [0 0.5 0.05 1.3 0], ...
    hp, k0, k0*sin(theta*pi/180)};

fprintf('  N     halfnpw   calls   mean [ms]   median [ms]\n');
for hp = [5 10 20 30]
    a = mkargs(hp);
    RTA_1d_tm(a{:});                 % warm up
    if hp <= 20, reps = 200; else, reps = 80; end
    ts = zeros(reps,1);
    for i = 1:reps
        t0 = tic; RTA_1d_tm(a{:}); ts(i) = toc(t0)*1e3;
    end
    fprintf('  %3d    %3d      %4d    %8.3f    %8.3f\n', 2*hp+1, hp, reps, mean(ts), median(ts));
end
