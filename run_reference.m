%% run_reference.m  —  PPML v3.0 1-D TM  REFERENCE RUNNER  (Phase 3)
% Runs the ORIGINAL, UNMODIFIED 1-D TM code (matlab_src/) over every fixture and
% writes reference_outputs/<engine>/<id>.mat (-v7) + reference_outputs/<engine>/env.json.
% Runs on BOTH MATLAB (reference of record) and Octave (cross-check); each writes to
% its own engine subdir so Phase 5 can diff them (cross-engine tier, VERIFICATION_NOTES
% §10.3).
%
% OUTCOME MODEL (§9a): every call wrapped try/catch (+ soft wall-clock cap). Recorded
% outcome is 'ok' | 'error' | 'timeout' with err_id/err_msg/elapsed. For 'ok' runs the
% harness also logs redheffer_cond (max_l cond(I - S2*S3)) -> §10.2 resonance tier.
% Out-of-domain fixtures' actual outcome is checked against their declared `expect`.

clear; clc;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(here);                                   % harness helpers + diagnostic
addpath(fullfile(here,'matlab_src','1d_tm'));    % RTA_1d_tm, SM_1d_tm, field_1d_tm
addpath(fullfile(here,'matlab_src','general'));  % sqrt_whittaker, smpropag_*_cond

isoct = (exist('OCTAVE_VERSION','builtin') ~= 0);
engine = ternary(isoct, 'octave', 'matlab');
outdir = fullfile(here, 'reference_outputs', engine);
if exist(outdir,'dir') ~= 7, mkdir(outdir); end

CAP = 60;    % soft wall-clock cap (s). NOTE: in-process; cannot hard-interrupt a hung
             % call. Adequate here (max npw=61 -> all calls are milliseconds). A hard
             % timeout would need subprocess isolation (documented limitation).

files = dir(fullfile(here,'fixtures','*.mat'));
names = sort({files.name});
nfx = numel(names);
fprintf('run_reference [%s]: %d fixtures -> %s\n', engine, nfx, outdir);

nok = 0; nerr = 0; ntmo = 0; nmis = 0;

for k = 1:nfx
    S = load(fullfile(here,'fixtures',names{k}));
    fx = S.fx;

    ref = struct();
    ref.id = fx.id; ref.func = fx.func; ref.engine = engine;
    ref.expect = fx.expect; ref.expect_error_id = fx.expect_error_id;
    ref.err_id = ''; ref.err_msg = ''; ref.elapsed_s = NaN;
    ref.warn_id = ''; ref.warn_msg = '';   % singular-matrix etc. warnings are part of the outcome
    ref.redheffer_cond = NaN;

    % ---- main call, wrapped -------------------------------------------------
    lastwarn('', '');       % clear so we capture only THIS call's warning
    t0 = tic;
    try
        switch fx.func
            case 'RTA'
                [RR,TT,AA] = RTA_1d_tm(fx.a,fx.L,fx.epssup,fx.epssub,...
                    fx.epsxA,fx.epszA,fx.epsxB,fx.epszB,fx.sigma,fx.f,fx.d,...
                    fx.halfnpw,fx.k0,fx.kpar);
                ref.RR = RR; ref.TT = TT; ref.AA = AA;
            case 'SM'
                [rl,rr,tlr,trl] = SM_1d_tm(fx.a,fx.L,fx.epssup,fx.epssub,...
                    fx.epsxA,fx.epszA,fx.epsxB,fx.epszB,fx.sigma,fx.f,fx.d,...
                    fx.halfnpw,fx.k0,fx.kpar);
                ref.rl = rl; ref.rr = rr; ref.tlr = tlr; ref.trl = trl;
            case 'field'
                [x,z,Ex,Ez,Sz] = field_1d_tm(fx.a,fx.L,fx.epssup,fx.epssub,...
                    fx.epsxA,fx.epszA,fx.epsxB,fx.epszB,fx.sigma,fx.f,fx.d,...
                    fx.halfnpw,fx.k0,fx.kpar,fx.nx,fx.nz);
                ref.x = x; ref.z = z; ref.Ex = Ex; ref.Ez = Ez; ref.Sz = Sz;
            otherwise
                error('runref:badfunc','unknown target function %s', fx.func);
        end
        ref.elapsed_s = toc(t0);
        [wmsg, wid] = lastwarn();          % e.g. 'MATLAB:singularMatrix' / Octave equiv
        ref.warn_msg = wmsg; ref.warn_id = wid;
        if ref.elapsed_s > CAP
            ref.outcome = 'timeout';
        else
            ref.outcome = 'ok';
        end
    catch err
        ref.elapsed_s = toc(t0);
        ref.outcome = 'error';
        ref.err_id  = err.identifier;
        ref.err_msg = err.message;
    end

    % ---- redheffer_cond diagnostic (independent of main outcome) ------------
    % Silence the mirror's expected near-cutoff "matrix singular" warnings (diagnostic
    % only; the main call above already ran with warnings live).
    ws = warning('off','all');
    ref.redheffer_cond = redheffer_cond_1d_tm(fx.a,fx.L,fx.epssup,fx.epssub,...
        fx.epsxA,fx.epszA,fx.epsxB,fx.epszB,fx.sigma,fx.f,fx.d,fx.halfnpw,fx.k0,fx.kpar);
    warning(ws);

    % ---- outcome vs declared expectation -----------------------------------
    switch fx.expect
        case 'ok',          exp_kind = 'ok';       % completes, physical
        case 'nonphysical', exp_kind = 'ok';       % completes, but garbage (no guard)
        case 'error',       exp_kind = 'error';
        otherwise,          exp_kind = 'ok';
    end
    actual_kind = ref.outcome;
    if strcmp(actual_kind,'timeout'), actual_kind = 'ok'; end  % ran, just slow
    match = strcmp(actual_kind, exp_kind);
    if match && strcmp(exp_kind,'error') && ~isempty(fx.expect_error_id)
        match = strcmp(ref.err_id, fx.expect_error_id);
    end
    ref.expect_match = match;

    % ---- tally + save ------------------------------------------------------
    switch ref.outcome
        case 'ok',      nok = nok + 1;
        case 'error',   nerr = nerr + 1;
        case 'timeout', ntmo = ntmo + 1;
    end
    if ~match, nmis = nmis + 1; end

    outfile = fullfile(outdir, [fx.id '.mat']);
    if isoct, save('-v7', outfile, 'ref'); else, save(outfile, 'ref', '-v7'); end

    if ~match
        fprintf('  [expect-mismatch] %-40s expect=%s actual=%s err_id=%s\n', ...
                fx.id, fx.expect, ref.outcome, ref.err_id);
    end
end

fprintf('\n== reference run [%s] ==\n', engine);
fprintf('  ok=%d  error=%d  timeout=%d  expect-mismatch=%d  total=%d\n', ...
        nok, nerr, ntmo, nmis, nfx);

% ---- env.json ---------------------------------------------------------------
try, bl = version('-blas'); catch, bl = ''; end
try, la = version('-lapack'); catch, la = ''; end
env = struct('schema_version',1,'phase',3,'engine',engine_string(),...
             'engine_key',engine,'platform',computer(),...
             'blas',bl,'lapack',la,'date',datestr(now,31),...
             'n_fixtures',nfx,'counts',struct('ok',nok,'error',nerr,'timeout',ntmo,...
             'expect_mismatch',nmis),...
             'vacuum_impedance_constant',376.730,...
             'note','Original PPML v3.0 1-D TM, unmodified. Outcome model per VERIFICATION_NOTES 9a.');
fid = fopen(fullfile(outdir,'env.json'),'w');
fprintf(fid,'%s\n', jsonencode(env));
fclose(fid);
fprintf('  wrote %d reference_outputs + env.json to %s\n', nfx, outdir);
