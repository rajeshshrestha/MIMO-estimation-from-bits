function [R_rec] = relax_2D(R,K)


%---------------- 2D Relax --------------------
y =R;
Fs        =  2*pi;
Nr        =  size(y,1);
Nt        =  size(y,2);
tr         =  (0:Nr-1)';
tt         =  (0:Nt-1)';
Lr         =  ceil( log2(16*Nr) );
Lt         =  ceil( log2(16*Nt) );
dFr        =  Fs/2^Lr;           % hertz
dFt        =  Fs/2^Lt;           % hertz
fr         =  0:dFr:Fs-dFr;
ft         =  0:dFt:Fs-dFt;
build_ar   =  @(fr) exp(1i*tr*fr);
build_at   =  @(ft) exp(1i*tt*ft);

% delta     =  dF;

freqr = zeros(K, 1);
freqt = zeros(K, 1);
amp  = zeros(K, 1);


%% initialization
ytot        =   zeros(Nr,Nt,K);
ytot(:,:,1)   =   y;
X           =   fft2(squeeze(ytot(:,:,1)), 2^Lr,2^Lt);
% [~,I]       =   max(.^2);
% [~,I] = max(abs(X),[],'all', 'linear');
[~,I] = max(abs(X(:)));

[I_row, I_col] = ind2sub(size(X),I) ;

ftemp       =   [fr(I_row),ft(I_col)];
% theta_pg = freqr(1);
% phi_pg = freqt(1);
[theta_pg,phi_pg] = relax_PG(y,Nr,Nt,ftemp(1),ftemp(2));
ar         = build_ar(theta_pg);
at         = build_at(phi_pg);
a           =   ar*at.';
amp(1)      =   trace(a'*squeeze(ytot(:,:,1)))/(Nr*Nt);
freqr(1) = theta_pg;
freqt(1) = phi_pg;
A1 = build_ar(freqr.');
A2 = build_at(freqt.');

%% start iterations
for k_bar = 2:K

    k_bar;

    relerr = 1;
    relerr_tot = 0;
    tt = 0;

    while relerr>1e-3 && tt<10

        err1 = norm(y - A1*diag(amp)*A2.')^2;

        for k2 = [k_bar, 1:k_bar-1]

            A1 = build_ar(freqr.');
            A2 = build_at(freqt.');

            ytot(:,:,k2)  =  y - A1(:, [1:k2-1,k2+1:k_bar])*diag(amp([1:k2-1, k2+1:k_bar]))*A2(:, [1:k2-1,k2+1:k_bar]).';
            X =  fft2(squeeze(ytot(:,:,k2)), 2^Lr,2^Lt);
            [~,I] = max(abs(X(:)));

            [I_row, I_col] = ind2sub(size(X),I) ;
            ftemp       =   [fr(I_row),ft(I_col)];
            [theta_pg,phi_pg] = relax_PG(squeeze(ytot(:,:,k2)),Nr,Nt,ftemp(1),ftemp(2));
            freqr(k2) = theta_pg;
            freqt(k2) = phi_pg;
            ar         = build_ar(theta_pg);
            at         = build_at(phi_pg);
            a           =   ar*at.';
            amp(k2)     =  trace(a'*squeeze(ytot(:,:,k2)))/(Nr*Nt);

        end

        err2 = norm(y - A1*diag(amp)*A2.')^2;
        relerr = abs(err1-err2)/err1;

        tt = tt + 1;

        relerr_tot = [relerr_tot, relerr];

    end

end
R_rec = A1*diag(amp)*A2.';