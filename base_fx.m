function fx = base_fx(func, id, rationale)
% Empty fixture struct with all PPML 1-D TM call arguments + bookkeeping fields.
fx = struct();
fx.id = id; fx.func = func; fx.rationale = rationale;
fx.category = id_category(id);
fx.a=[]; fx.L=[]; fx.epssup=[]; fx.epssub=[];
fx.epsxA=[]; fx.epszA=[]; fx.epsxB=[]; fx.epszB=[];
fx.sigma=[]; fx.f=[]; fx.d=[]; fx.halfnpw=[]; fx.k0=[]; fx.kpar=[];
fx.lambda0=[]; fx.theta_deg=[]; fx.nx=[]; fx.nz=[];
% Outcome expectation (first-class): 'ok' | 'error' | 'nonphysical'.
%   'error'       -> run is expected to raise; expect_error_id names the identifier.
%   'nonphysical' -> run returns a value but it is physically invalid (no guard in
%                    the original); port reproduces by default (bug-for-bug).
fx.expect='ok'; fx.expect_error=false; fx.expect_error_id=''; fx.expect_nonphysical=false;
end
