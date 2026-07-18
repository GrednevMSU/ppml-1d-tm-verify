function [anom, nprop_sup, nprop_sub] = order_metrics(a, epssup, epssub, k0, kpar, halfnpw)
% Rayleigh/Wood-anomaly proximity + propagating-order counts on BOTH sides.
%
% SINGLE, AUTHORITATIVE anomaly metric (this file + the manifest are the only
% definition; every report section references it, no second formula anywhere):
%
%   anom = min over diffraction orders m, over each REAL strate, of
%              | Re(eps)*k0^2 - kx_m^2 |  /  ( Re(eps)*k0^2 )
%          i.e. the normalized distance of the longitudinal argument to its
%          cutoff. anom -> 0  <=>  some order's kz -> 0 (Rayleigh/Wood anomaly).
%          We use |argument|, NOT |Re(kz)|: deep-evanescent orders have Re(kz)=0
%          yet are far from cutoff, so |Re(kz)| would false-flag every fixture.
%
%   nprop_sup -> # propagating orders in the (real) superstrate.
%   nprop_sub -> # propagating orders in the substrate if it is real & >0,
%                else -1 (undefined: complex/absorbing substrate).
%
% epssup is real by the RTA/SM contract.

n  = -halfnpw:halfnpw;
kx = (2*pi/a) * n + kpar;

rsup    = real(epssup);
argsup  = rsup*k0^2 - kx.^2;
metrics = abs(argsup) / (rsup*k0^2);

if imag(epssub) == 0 && real(epssub) > 0
    rsub      = real(epssub);
    argsub    = rsub*k0^2 - kx.^2;
    metrics   = [metrics, abs(argsub) / (rsub*k0^2)];
    nprop_sub = sum(argsub > 0);
else
    nprop_sub = -1;   % undefined for complex/absorbing substrate
end

anom      = min(metrics);
nprop_sup = sum(argsup > 0);
end
