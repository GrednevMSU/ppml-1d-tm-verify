function rc = redheffer_cond_1d_tm(a,L,epssup,epssub,epsxA,epszA,epsxB,epszB,sigma,f,d,halfnpw,k0,kpar)
% DIAGNOSTIC (harness-only, non-gating): max over layers of cond(I - S2*S3), the
% Redheffer star-product inversion condition number. Feeds the §10.2 resonance tier.
%
% This is a faithful MIRROR of the S-matrix core of SM_1d_tm.m (lines 76-148): it
% reuses the SAME general/ helpers (sqrt_whittaker, smpropag_fw_cond, smpropag_bw_cond)
% and identical formulas, so it cannot diverge from the original build. It does NOT
% modify or replace the original code and produces NO gated output — only the scalar
% diagnostic. Returns NaN if the build fails (e.g. out-of-domain inputs).

try
    n = -halfnpw:halfnpw;
    npw = size(n,2);
    kx = (2*pi/a)*n + kpar;

    q   = zeros(npw,npw,L+2);
    phi = zeros(npw,npw,L+2);
    A   = zeros(npw,npw,L+2);
    epsx = zeros(npw,npw,L+2);
    etaz = zeros(npw,npw,L+2);

    q(:,:,1)   = diag(sqrt_whittaker(epssup*k0^2 - kx.*kx));
    phi(:,:,1) = eye(npw);
    A(:,:,1)   = diag(k0^2 - kx.*kx/epssup)/q(:,:,1);
    q(:,:,L+2)   = diag(sqrt_whittaker(epssub*k0^2 - kx.*kx));
    phi(:,:,L+2) = eye(npw);
    A(:,:,L+2)   = diag(k0^2 - kx.*kx/epssub)/q(:,:,L+2);

    if L > 0
    for l = 1:L
        F = zeros(npw);
        for i = 1:npw
        for j = 1:npw
            if (i == j)
                F(i,j) = f(l);
            else
                F(i,j) = sin(pi*f(l)*(n(i)-n(j)))/(pi*(n(i)-n(j)));
            end
        end
        end
        epsx(:,:,l+1) = (1/epsxB(l) - 1/epsxA(l))*F + 1/epsxA(l)*eye(npw);
        etaz(:,:,l+1) = (  epszB(l) -   epszA(l))*F +   epszA(l)*eye(npw);
        [phhi,qq] = eig(epsx(:,:,l+1)\(k0^2*eye(npw) - diag(kx)*(etaz(:,:,l+1)\diag(kx))));
        q(:,:,l+1) = diag(sqrt_whittaker(diag(qq)));
        phi(:,:,l+1) = phhi;
        A(:,:,l+1) = (k0^2*eye(npw) - diag(kx)*(etaz(:,:,l+1)\diag(kx)))*phhi/q(:,:,l+1);
    end
    end

    S1 = zeros(npw,npw,L+2); S2 = zeros(npw,npw,L+2);
    S1(:,:,1) = eye(npw);
    for l = 1:L+1
    [S1(:,:,l+1),S2(:,:,l+1)] = smpropag_fw_cond(S1(:,:,l),S2(:,:,l),...
        phi(:,:,l),phi(:,:,l+1),A(:,:,l),A(:,:,l+1),...
        exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l+1)*d(l+1))),sigma(l)*376.730/k0);
    end

    S3 = zeros(npw,npw,L+2); S4 = zeros(npw,npw,L+2);
    S4(:,:,L+2) = eye(npw);
    for ll = 1:L+1
        l = L+3-ll;
    [S3(:,:,l-1),S4(:,:,l-1)] = smpropag_bw_cond(S3(:,:,l),S4(:,:,l),...
        phi(:,:,l),phi(:,:,l-1),A(:,:,l),A(:,:,l-1),...
        exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l-1)*d(l-1))),sigma(l-1)*376.730/k0);
    end

    rc = 0;
    for l = 1:L+2
        rc = max(rc, cond(eye(npw) - S2(:,:,l)*S3(:,:,l)));
    end
catch
    rc = NaN;
end
end
