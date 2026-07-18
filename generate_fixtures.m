%% generate_fixtures.m  —  PPML v3.0 1-D TM differential-test corpus
% MATLAB/Octave compatible. Writes fixtures/*.mat (-v7) + fixtures/manifest.json.
%
% Corpus = scenario replication (from examples/) + branch/edge cases (every entry
% in the BRANCH MAP) + synthetic conducting-interface (sigma~=0) + randomized sweep
% (N>=200, fixed disclosed seed, WIDE/STRESS ranges per user decision).
%
% Fixtures store INPUTS only; reference outputs are produced by run_reference.m
% (Phase 3). Each fixture .mat contains one struct `fx` with all call arguments and
% bookkeeping tags. manifest.json records id/function/ranges/rationale/seed +
% SHA-256 of every fixture file (auditability: corpus cannot be cherry-picked later).
%
% Run once; commit the resulting fixtures/ + manifest.json. Cross-engine RNG identity
% is NOT required — the frozen files + hashes are the ground truth both runners load.

clear; clc;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
fixdir = fullfile(here, 'fixtures');
if exist(fixdir, 'dir') ~= 7, mkdir(fixdir); end
addpath(here);   % sha256_file, order_metrics

SEED = 20260717;
if (exist('rng', 'builtin') == 5) || (exist('rng', 'file') == 2)
    rng(SEED, 'twister');
else
    rand('state', SEED); randn('state', SEED);   %#ok<RAND>  Octave fallback
end

% Documented sampling ranges (WIDE / stress-mode) -> echoed into manifest.
RANGES = struct( ...
   'lambda0_um',   [0.5 40], ...     % vacuum wavelength -> k0 = 2*pi/lambda0
   'a_um',         [0.2 40], ...     % period: subwavelength through several orders
   'theta_deg',    [0 80], ...       % incl. normal (0) and steep
   'L_set',        [0 1 2 3 4], ...
   'f',            [0.05 0.95], ...
   'd_internal_um',[0.02 5], ...     % thick metal branch pushes to [2 8]
   'eps_diel_re',  [1.5 13], 'eps_diel_im', [0 0.5], ...
   'eps_metal_re', [-4000 -2], 'eps_metal_im', [0.1 500], ...
   'epssup_set',   [1 2.25 4 11], ...
   'halfnpw_ladder',[0 1 2 3 5 8 10 15 20 30], ...
   'sigma_nonzero_frac', 0.20, 'sigma_mag', [1e-5 1e-2], ...
   'anomaly_threshold', 1e-2);       % |Re(kz)|/k0 < this  ==> near-anomaly tier

F = {};   % cell array of fixture structs

% ---- small inline material sampler -------------------------------------------
mk_diel  = @() (RANGES.eps_diel_re(1)  + diff(RANGES.eps_diel_re)*rand)  + 1i*(RANGES.eps_diel_im(1)  + diff(RANGES.eps_diel_im)*rand);
mk_metal = @() (RANGES.eps_metal_re(1) + diff(RANGES.eps_metal_re)*rand) + 1i*(RANGES.eps_metal_im(1) + diff(RANGES.eps_metal_im)*rand);
pickrng  = @(r) r(1) + diff(r)*rand;

% ============================================================================ %
%  A. SCENARIO REPLICATION  (anchored to examples/, real dispersive materials)
% ============================================================================ %

% ---- S1: OpticalCritical_APL2013 (RTA), a=30, L=4, halfnpw=20, oblique MIM ----
for pt = [struct('nu',3.0,'th',30); struct('nu',3.5,'th',45)].'
    nu = pt.nu; th = pt.th; h = 4.67;
    wn = nu*100/3;
    epsAu   = 1 - (7.25e4)^2/(wn^2 + 1i*216*wn);
    epsGaAs = 11*(1 + (292^2 - 268^2)/(268^2 - wn^2 - 1i*2.4*wn));
    om = 2*pi*nu*1e12;
    epsGaAs_dop = 11 - (2.6e24*1.6e-19^2/(8.8e-12*0.067*9.1e-31*(om^2 + 1i*om/100e-15)));
    fx = base_fx('RTA', sprintf('scenario_S1_OpticalCritical_nu%0.1f_th%d', nu, th), 'APL2013 fig3a: oblique MIM, lossy Au sub, doped GaAs');
    fx.a = 30; fx.L = 4; fx.halfnpw = 20;
    fx.f     = [.77 .77 1 1];
    fx.d     = [30 .2 h (9-h) 1 30];
    fx.epsxA = [1 1 1 1];             fx.epszA = [1 1 1 1];
    fx.epsxB = [epsAu epsGaAs epsGaAs epsGaAs_dop];
    fx.epszB = fx.epsxB;
    fx.sigma = zeros(1,5);
    fx.epssup = 1; fx.epssub = epsAu;
    fx.lambda0 = []; fx.k0 = 2*pi*nu*0.01/3; fx.theta_deg = th;
    fx.kpar = fx.k0*sin(th*pi/180);
    F{end+1} = fx; %#ok<*SAGROW>
end

% ---- S2/S3/S4: UniversalLineshapes_SciRep2016, a=3.2, L=3, halfnpw=10 ---------
epsAu10 = -4000 + 300*1i; epsbg = 10.05;
for cfg = [struct('om',150,'n2',3e11,'fn','RTA'); ...
           struct('om',130,'n2',3e11,'fn','RTA'); ...
           struct('om',150,'n2',0,   'fn','RTA'); ...
           struct('om',150,'n2',3e11,'fn','SM')].'
    omv = cfg.om; k0 = omv*(2*pi/1240); theta = 0.1;
    omp = sqrt(4*pi*4.8e-10^2*cfg.n2/(epsbg*0.067*9.1e-28*21.8e-7))*6.58e-13;
    om12 = 150; gamma12 = 3;
    epszn = epsbg/(1 - omp^2/(om12^2 - omv^2 - 1i*2*gamma12*omv));
    tag = sprintf('scenario_S%s_UniversalLineshapes_%s_om%d_n%g', ...
                  ternary(strcmp(cfg.fn,'SM'),'4','2'), cfg.fn, omv, cfg.n2);
    fx = base_fx(cfg.fn, tag, 'SciRep2016: near-normal thin metal grating, anisotropic eps_z (MQW)');
    fx.a = 3.2; fx.L = 3; fx.halfnpw = 10;
    fx.f     = [0 0.73 0];
    fx.d     = [0 0.5 0.05 1.3 0];
    fx.epsxA = [epsbg epsbg epsbg];   fx.epszA = [epsbg epsbg epszn];
    fx.epsxB = [epsbg epsAu10 epsbg]; fx.epszB = [epsbg epsAu10 epszn];
    fx.sigma = zeros(1,4);
    fx.epssup = 1; fx.epssub = 1;
    fx.lambda0 = []; fx.k0 = k0; fx.theta_deg = theta;
    fx.kpar = k0*sin(theta*pi/180);
    F{end+1} = fx;
end

% ============================================================================ %
%  B. BRANCH / EDGE FIXTURES  (every BRANCH MAP entry gets >=1)
% ============================================================================ %
% Simple lossless dielectric base so energy R+T=1 exactly (independent check).
DA = 12.0; DB = 2.25;   % lossless dielectrics (Si-ish / SiO2-ish, real)

% B1a: L=0 bare interface (skip internal-layer block)
fx = base_fx('RTA', 'branch_B1_L0_bare', 'L=0: bare super/sub interface, internal-layer block skipped');
fx.a = 2; fx.L = 0; fx.halfnpw = 5;
fx.f = []; fx.d = [1 1]; fx.epsxA = []; fx.epszA = []; fx.epsxB = []; fx.epszB = [];
fx.sigma = 0; fx.epssup = 1; fx.epssub = DB;
fx.k0 = 2*pi/1.5; fx.theta_deg = 20; fx.kpar = fx.k0*sin(20*pi/180); fx.lambda0 = 1.5;
F{end+1} = fx;

% B2: halfnpw=0 -> npw=1, TMM degenerate (only Fourier diagonal path)
fx = base_fx('RTA', 'branch_B2_halfnpw0_TMM', 'halfnpw=0: 1x1, reduces to ordinary 2x2 multilayer TMM');
fx.a = 1; fx.L = 1; fx.halfnpw = 0;
fx.f = 0.5; fx.d = [1 0.4 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
fx.sigma = [0 0]; fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.0; fx.theta_deg = 15; fx.kpar = fx.k0*sin(15*pi/180); fx.lambda0 = 1.0;
F{end+1} = fx;

% B2b: halfnpw=1 -> npw=3 minimal patterned
fx = base_fx('RTA', 'branch_B2b_halfnpw1', 'halfnpw=1 (npw=3): minimal patterned grating');
fx.a = 1.2; fx.L = 1; fx.halfnpw = 1;
fx.f = 0.4; fx.d = [1 0.5 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
fx.sigma = [0 0]; fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.0; fx.theta_deg = 10; fx.kpar = fx.k0*sin(10*pi/180); fx.lambda0 = 1.0;
F{end+1} = fx;

% B3a: normal incidence kpar=0 (theta=0) — special path (docs warn "never 0")
fx = base_fx('RTA', 'branch_B3_normal_incidence', 'theta=0 (kpar=0): normal-incidence path (docstring warns against it)');
fx.a = 1.5; fx.L = 1; fx.halfnpw = 6;
fx.f = 0.5; fx.d = [1 0.6 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
fx.sigma = [0 0]; fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.0; fx.theta_deg = 0; fx.kpar = 0; fx.lambda0 = 1.0;
F{end+1} = fx;

% B3b: single-order (subwavelength a<lambda): only 0 order propagates
fx = base_fx('RTA', 'branch_singleorder_subwl', 'subwavelength a<lambda: single propagating order');
fx.a = 0.4; fx.L = 1; fx.halfnpw = 8;
fx.f = 0.5; fx.d = [1 0.3 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
fx.sigma = [0 0]; fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/2.0; fx.theta_deg = 12; fx.kpar = fx.k0*sin(12*pi/180); fx.lambda0 = 2.0;
F{end+1} = fx;

% B3c: multi-order (a>>lambda): several propagating orders
fx = base_fx('RTA', 'branch_multiorder', 'a>>lambda: several propagating diffraction orders');
fx.a = 8.0; fx.L = 1; fx.halfnpw = 15;
fx.f = 0.5; fx.d = [1 0.5 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
fx.sigma = [0 0]; fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.0; fx.theta_deg = 25; fx.kpar = fx.k0*sin(25*pi/180); fx.lambda0 = 1.0;
F{end+1} = fx;

% B4: energy-conservation guard exercised by every RTA fixture; add a purely
% lossless multi-layer where R+T must equal 1 to machine precision.
fx = base_fx('RTA', 'branch_B4_lossless_energy', 'lossless stack: independent invariant R+T=1 (guard at RTA.m:183)');
fx.a = 2.0; fx.L = 3; fx.halfnpw = 8;
fx.f = [0.3 0.6 0.5]; fx.d = [1 0.4 0.7 0.5 1];
fx.epsxA = [DA DA DA]; fx.epszA = [DA DA DA];
fx.epsxB = [DB DB DB]; fx.epszB = [DB DB DB];
fx.sigma = zeros(1,4); fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.3; fx.theta_deg = 18; fx.kpar = fx.k0*sin(18*pi/180); fx.lambda0 = 1.3;
F{end+1} = fx;

% Thick metallic layer — S-matrix stability (transfer-matrix would overflow)
fx = base_fx('RTA', 'branch_thick_metal_Smatrix', 'thick lossy metal: S-matrix must stay bounded (no transfer-matrix overflow)');
fx.a = 1.0; fx.L = 1; fx.halfnpw = 10;
fx.f = 0.5; fx.d = [1 6.0 1];
fx.epsxA = DB; fx.epszA = DB; fx.epsxB = -50 + 5i; fx.epszB = -50 + 5i;
fx.sigma = [0 0]; fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.0; fx.theta_deg = 20; fx.kpar = fx.k0*sin(20*pi/180); fx.lambda0 = 1.0;
F{end+1} = fx;

% Wood/Rayleigh anomaly: tune theta so the +1 order kz -> 0 in the superstrate.
% Cutoff: kpar + 2*pi/a = k0*sqrt(epssup)  =>  sin(theta) = 1 - lambda/a.
alam = 1.0; aa = 1.6; sin_cut = 1 - alam/aa;          % order m=+1 grazing
th_cut = asin(max(min(sin_cut,0.999999),-0.999999))*180/pi;
for dth = [0, 0.15, -0.5]     % on cutoff, just above, well below
    fx = base_fx('RTA', sprintf('branch_wood_anomaly_%+0.2f', dth), ...
                 'near Rayleigh/Wood anomaly: an order kz crosses zero -> looser/INFO tier (physics, not port error)');
    fx.a = aa; fx.L = 1; fx.halfnpw = 12;
    fx.f = 0.5; fx.d = [1 0.5 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
    fx.sigma = [0 0]; fx.epssup = 1; fx.epssub = 1;
    fx.k0 = 2*pi/alam; th = th_cut + dth;
    fx.theta_deg = th; fx.kpar = fx.k0*sin(th*pi/180); fx.lambda0 = alam;
    F{end+1} = fx;
end

% Immersion: epssup>1 real (allowed), epssub complex (metal) -> TT = substrate absorbance
fx = base_fx('RTA', 'branch_epssup_gt1_metal_sub', 'epssup=4 (real immersion), epssub metallic complex: TT = substrate absorbance');
fx.a = 1.0; fx.L = 1; fx.halfnpw = 10;
fx.f = 0.5; fx.d = [1 0.4 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
fx.sigma = [0 0]; fx.epssup = 4; fx.epssub = -30 + 3i;
fx.k0 = 2*pi/1.0; fx.theta_deg = 15; fx.kpar = fx.k0*sqrt(4)*sin(15*pi/180); fx.lambda0 = 1.0;
F{end+1} = fx;

% SM branch coverage (complex coefficients) + field branch coverage
fx = base_fx('SM', 'branch_SM_dielectric', 'SM 0-order complex coeffs, lossless dielectric grating');
fx.a = 1.5; fx.L = 2; fx.halfnpw = 8;
fx.f = [0.4 0.6]; fx.d = [1 0.5 0.4 1];
fx.epsxA = [DA DA]; fx.epszA = [DA DA]; fx.epsxB = [DB DB]; fx.epszB = [DB DB];
fx.sigma = zeros(1,3); fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.1; fx.theta_deg = 22; fx.kpar = fx.k0*sin(22*pi/180); fx.lambda0 = 1.1;
F{end+1} = fx;

% SM fixture that IS unitarity-eligible: subwavelength (only 0 order both sides),
% lossless, sigma=0 -> the 0-order 2x2 S must be unitary (INFO invariant target).
fx = base_fx('SM', 'branch_SM_unitary_singleorder', 'SM single-order lossless: 2x2 S unitarity INFO check is valid here');
fx.a = 0.5; fx.L = 1; fx.halfnpw = 8;
fx.f = 0.5; fx.d = [1 0.4 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
fx.sigma = [0 0]; fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.5; fx.theta_deg = 20; fx.kpar = fx.k0*sin(20*pi/180); fx.lambda0 = 1.5;
F{end+1} = fx;

fx = base_fx('field', 'branch_field_maps', 'field_1d_tm: Ex/Ez/Sz maps on sampled grid, lossy grating');
fx.a = 1.5; fx.L = 2; fx.halfnpw = 8;
fx.f = [0.5 0.5]; fx.d = [0.5 0.4 0.6 0.5];
fx.epsxA = [DA DA]; fx.epszA = [DA DA]; fx.epsxB = [DB -20+2i]; fx.epszB = [DB -20+2i];
fx.sigma = zeros(1,3); fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.0; fx.theta_deg = 17; fx.kpar = fx.k0*sin(17*pi/180); fx.lambda0 = 1.0;
fx.nx = 21; fx.nz = [6 10 10 6];
F{end+1} = fx;

% ============================================================================ %
%  C. SYNTHETIC CONDUCTING INTERFACES  (sigma ~= 0)  — v3.0 feature, no demo covers it
% ============================================================================ %
% sigma enters as sigma*376.730/k0. Choose sigma so g = sigma*Z0/k0 spans ~0.05..2.
Z0 = 376.730;
for spec = [struct('sig', 5e-4,       'tag','real_small'); ...
            struct('sig', 5e-3,       'tag','real_mid'); ...
            struct('sig', 2e-2,       'tag','real_large'); ...
            struct('sig', 3e-3+3e-3i, 'tag','complex')].'
    fx = base_fx('RTA', ['sigma_' spec.tag], ...
        sprintf('conducting interface sigma=%s (g=sigma*Z0/k0), v3.0 feature untested by demos', num2str(spec.sig)));
    fx.a = 1.5; fx.L = 2; fx.halfnpw = 8;
    fx.f = [0.5 0.5]; fx.d = [1 0.4 0.5 1];
    fx.epsxA = [DA DA]; fx.epszA = [DA DA]; fx.epsxB = [DB DB]; fx.epszB = [DB DB];
    % nonzero sigma on the middle internal interface (index 2 of L+1=3)
    fx.sigma = [0 spec.sig 0];
    fx.epssup = 1; fx.epssub = 1;
    fx.k0 = 2*pi/1.2; fx.theta_deg = 20; fx.kpar = fx.k0*sin(20*pi/180); fx.lambda0 = 1.2;
    F{end+1} = fx;
end
% sigma on a strate interface + SM readout
fx = base_fx('SM', 'sigma_strate_SM', 'conducting interface at superstrate boundary, SM coeffs');
fx.a = 1.5; fx.L = 1; fx.halfnpw = 6;
fx.f = 0.5; fx.d = [1 0.4 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
fx.sigma = [4e-3 0]; fx.epssup = 1; fx.epssub = 1;
fx.k0 = 2*pi/1.2; fx.theta_deg = 20; fx.kpar = fx.k0*sin(20*pi/180); fx.lambda0 = 1.2;
F{end+1} = fx;

% ============================================================================ %
%  D. OUT-OF-DOMAIN / NEGATIVE FIXTURES  (expected error / ill-defined)
% ============================================================================ %
% epssup complex is docstring-forbidden. EMPIRICALLY CONFIRMED (Phase-3 run_reference,
% Octave): ALL THREE functions SILENTLY return nonphysical results (complex RR/TT).
% RTA does NOT throw — its energy "guard" (RTA.m:183) is a telescoping tautology
% (RR+TT+ΣAA ≡ 1 by construction) that cannot catch a nonphysical input. See
% ORIGINAL_CODE_FINDINGS F2 (revised) + F3 (vacuous guard). The earlier hypothesis
% "RTA throws RTA:EnNotCons" was refuted by the outcome model — kept as a worked
% example of hypothesis-vs-data. Detectable signature: |Im(RR)| large (~0.094 here)
% vs ~0 for every physical fixture; the independent realness invariant (§9) uses this.
epssup_cplx = 1 + 0.3i;
for spec = [struct('fn','RTA',  'exp','nonphysical', 'eid',''); ...
            struct('fn','SM',   'exp','nonphysical', 'eid',''); ...
            struct('fn','field','exp','nonphysical', 'eid','')].'
    fx = base_fx(spec.fn, ['out_of_domain_epssup_complex_' spec.fn], ...
        'epssup complex (docstring-forbidden): all three silently return nonphysical (complex RR/TT). No guard. See ORIGINAL_CODE_FINDINGS F2/F3.');
    fx.expect = spec.exp;
    fx.expect_error = strcmp(spec.exp,'error');
    fx.expect_error_id = spec.eid;
    fx.expect_nonphysical = strcmp(spec.exp,'nonphysical');
    fx.a = 1.5; fx.L = 1; fx.halfnpw = 6;
    fx.f = 0.5; fx.d = [1 0.4 1]; fx.epsxA = DA; fx.epszA = DA; fx.epsxB = DB; fx.epszB = DB;
    fx.sigma = [0 0]; fx.epssup = epssup_cplx; fx.epssub = 1;
    fx.k0 = 2*pi/1.0; fx.theta_deg = 15; fx.kpar = fx.k0*sin(15*pi/180); fx.lambda0 = 1.0;
    if strcmp(spec.fn,'field'), fx.nx = 15; fx.nz = [4 6 4]; end
    F{end+1} = fx;
end

% ============================================================================ %
%  E. RANDOMIZED SWEEP  (N>=200, fixed seed, WIDE/STRESS ranges)
% ============================================================================ %
NRAND = 220;
ladder = RANGES.halfnpw_ladder;
for r = 1:NRAND
    % choose target function: mostly RTA (richest: R/T/A + energy), some SM/field
    u = rand;
    if     u < 0.68, fn = 'RTA';
    elseif u < 0.85, fn = 'SM';
    else             fn = 'field';
    end
    L = RANGES.L_set(randi(numel(RANGES.L_set)));
    lambda0 = pickrng(RANGES.lambda0_um);
    k0 = 2*pi/lambda0;
    a  = pickrng(RANGES.a_um);
    epssup = RANGES.epssup_set(randi(numel(RANGES.epssup_set)));
    theta = pickrng(RANGES.theta_deg);
    kpar = k0*sqrt(epssup)*sin(theta*pi/180);
    halfnpw = ladder(randi(numel(ladder)));

    % substrate: dielectric or metal (complex allowed)
    if rand < 0.5, epssub = mk_diel(); else, epssub = mk_metal(); end

    % internal layers
    epsxA = zeros(1,L); epszA = zeros(1,L); epsxB = zeros(1,L); epszB = zeros(1,L);
    f = zeros(1,L); d_int = zeros(1,L);
    for l = 1:L
        if rand < 0.5, eA = mk_diel(); else, eA = mk_metal(); end
        if rand < 0.5, eB = mk_diel(); else, eB = mk_metal(); end
        epsxA(l) = eA; epsxB(l) = eB;
        % 30% anisotropic: eps_z differs
        if rand < 0.3, epszA(l) = mk_diel(); else, epszA(l) = eA; end
        if rand < 0.3, epszB(l) = mk_metal(); else, epszB(l) = eB; end
        f(l) = pickrng(RANGES.f);
        if rand < 0.15, d_int(l) = pickrng([2 8]); else, d_int(l) = pickrng(RANGES.d_internal_um); end
    end
    d = [1, d_int, 1];

    % conducting interfaces
    sigma = zeros(1, L+1);
    if rand < RANGES.sigma_nonzero_frac
        idx = randi(L+1);
        mag = 10^(log10(RANGES.sigma_mag(1)) + rand*(log10(RANGES.sigma_mag(2))-log10(RANGES.sigma_mag(1))));
        if rand < 0.5, sigma(idx) = mag; else, sigma(idx) = mag*(1 + 1i); end
    end

    fx = base_fx(fn, sprintf('rand_%03d', r), sprintf('randomized sweep (seed %d), wide/stress ranges', SEED));
    fx.a=a; fx.L=L; fx.halfnpw=halfnpw; fx.f=f; fx.d=d;
    fx.epsxA=epsxA; fx.epszA=epszA; fx.epsxB=epsxB; fx.epszB=epszB;
    fx.sigma=sigma; fx.epssup=epssup; fx.epssub=epssub;
    fx.lambda0=lambda0; fx.k0=k0; fx.theta_deg=theta; fx.kpar=kpar;
    if strcmp(fn,'field')
        fx.nx = 15; nz = 2 + randi(6, 1, L+2); fx.nz = nz;
    end
    F{end+1} = fx;
end

% ============================================================================ %
%  WRITE FIXTURES + MANIFEST
% ============================================================================ %
nfx = numel(F);
scenarios = {'S1_OpticalCritical','S2_UniversalLineshapes','S4_SM','S5_TMM_halfnpw0','S6_field'};
cov = struct(); for s=1:numel(scenarios), cov.(scenarios{s}) = {}; end
man_items = {};

% Uniqueness guard: duplicate ids would silently overwrite .mat files.
allids = cellfun(@(x) x.id, F, 'UniformOutput', false);
if numel(unique(allids)) ~= numel(allids)
    error('generate_fixtures:dupID', 'Duplicate fixture ids detected — ids must be unique.');
end

for k = 1:nfx
    fx = F{k};
    % anomaly + propagating-order counts (single authoritative metric in order_metrics)
    [anom, nprop_sup, nprop_sub] = order_metrics(fx.a, fx.epssup, fx.epssub, fx.k0, fx.kpar, fx.halfnpw);
    fx.anomaly_metric = anom;
    fx.n_prop_sup = nprop_sup;
    fx.n_prop_sub = nprop_sub;
    fx.near_anomaly = (anom < RANGES.anomaly_threshold);
    if nprop_sup <= 1, fx.regime = 'single-order'; else, fx.regime = 'multi-order'; end
    % lossless = every permittivity real (dielectric-only, no metal/absorption)
    fx.lossless = all(imag([fx.epssup, fx.epssub, fx.epsxA(:).', fx.epszA(:).', fx.epsxB(:).', fx.epszB(:).']) == 0);
    % Unitarity of the 0-order 2x2 S holds ONLY for exactly one propagating order on
    % BOTH sides, sigma=0, lossless. Outside that the 2x2 block is legitimately
    % non-unitary (energy leaks to nonzero orders) — flag guards compare.py from
    % false-reddening on valid physics.
    fx.unitarity_eligible = fx.lossless && (nprop_sup == 1) && (nprop_sub == 1) && all(fx.sigma(:) == 0);
    fx.seed = SEED;

    fn_mat = fullfile(fixdir, [fx.id '.mat']);
    save_v7(fn_mat, fx);
    sha = sha256_file(fn_mat);

    % coverage tagging (heuristic map to documented scenarios)
    if ~isempty(strfind(fx.id,'S1_Optical')), cov.S1_OpticalCritical{end+1}=fx.id; end
    if ~isempty(strfind(fx.id,'UniversalLineshapes')) && strcmp(fx.func,'RTA'), cov.S2_UniversalLineshapes{end+1}=fx.id; end
    if strcmp(fx.func,'SM'), cov.S4_SM{end+1}=fx.id; end
    if fx.halfnpw==0, cov.S5_TMM_halfnpw0{end+1}=fx.id; end
    if strcmp(fx.func,'field'), cov.S6_field{end+1}=fx.id; end

    item = struct( ...
        'id', fx.id, 'function', fx.func, 'category', fx.category, ...
        'rationale', fx.rationale, 'seed', SEED, ...
        'a_um', fx.a, 'L', fx.L, 'halfnpw', fx.halfnpw, ...
        'k0', fx.k0, 'theta_deg', fx.theta_deg, 'kpar', real(fx.kpar), ...
        'epssup_re', real(fx.epssup), 'epssup_im', imag(fx.epssup), ...
        'n_prop_sup', nprop_sup, 'n_prop_sub', nprop_sub, 'regime', fx.regime, ...
        'anomaly_metric', anom, 'near_anomaly', fx.near_anomaly, ...
        'lossless', fx.lossless, 'unitarity_eligible', fx.unitarity_eligible, ...
        'sigma_nonzero', any(fx.sigma ~= 0), ...
        'expect', fx.expect, 'expect_error', isfield_true(fx,'expect_error'), ...
        'expect_error_id', fx.expect_error_id, 'expect_nonphysical', fx.expect_nonphysical, ...
        'sha256', sha);
    man_items{end+1} = item;
end

% coverage matrix -> console
fprintf('\n================ SCENARIO COVERAGE MATRIX ================\n');
for s = 1:numel(scenarios)
    ids = cov.(scenarios{s});
    fprintf('  %-26s : %3d fixtures\n', scenarios{s}, numel(ids));
    if numel(ids)==0
        fprintf('     *** UNCOVERED -> declare in Scope & Limitations ***\n');
    end
end
fprintf('  %-26s : %3d fixtures (near Wood anomaly, looser/INFO tier)\n', 'near_anomaly', ...
        sum(cellfun(@(it) it.near_anomaly, man_items)));
fprintf('  %-26s : %3d fixtures (sigma~=0 conducting interface)\n', 'sigma_nonzero', ...
        sum(cellfun(@(it) it.sigma_nonzero, man_items)));
fprintf('  %-26s : %3d total\n', 'TOTAL', nfx);
fprintf('=========================================================\n');

% manifest.json  (build scalar struct first; assign struct-array field AFTER to
% avoid the struct() cell-expansion pitfall)
% corpus_version: bump on ANY change to the frozen corpus; log it in CORPUS_CHANGELOG.md.
manifest = struct('schema_version', 1, 'corpus_version', 3, ...
                  'generated', datestr(now,31), ...
                  'seed', SEED, 'engine', engine_string(), ...
                  'n_fixtures', nfx, 'ranges', RANGES, ...
                  'changelog', 'see CORPUS_CHANGELOG.md', ...
                  'pass_logic', 'OR(rel<=tol_rel, abs<=tol_abs) — see report header');
manifest.fixtures = cell2mat_struct(man_items);
jtxt = jsonencode(manifest);
fid = fopen(fullfile(fixdir, 'manifest.json'), 'w');
fprintf(fid, '%s\n', jtxt);
fclose(fid);
fprintf('Wrote %d fixtures + manifest.json to %s\n', nfx, fixdir);

% Helper functions live in separate .m files (cross-engine: base_fx, id_category,
% isfield_true, ternary, engine_string, cell2mat_struct, save_v7, order_metrics,
% sha256_file) — Octave does not resolve script-local functions when run as a script.
