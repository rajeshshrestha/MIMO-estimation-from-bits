function H_admm=newADMM(z, D, sigma_true, scaler_for_quantization ,Nr, Nt, K)

% Add necessary files to path
mfilepath = fileparts(which(mfilename));

% Find the up and low bin threshold
b_up = z+1;
b_low = z-1;


sigma = max(sigma_true*scaler_for_quantization.scale_ratio,1);

Lambda = zeros(Nr,Nt*2);
pho=str2double(getenv("rho"));
if isnan(pho)
    pho = 1;  % ADMM parameter
end

U_ini_real = zeros(Nr,Nt*2);
U_buff=  U_ini_real(:,1:Nt) + 1i* U_ini_real(:,Nt+1:end);

u_ini = [U_ini_real(:,1:Nt), -U_ini_real(:,Nt+1:end);
    U_ini_real(:,Nt+1:end), U_ini_real(:,1:Nt)];
u_ini = u_ini(:);


num_iter = 100;

for i_ADMM =1:num_iter
    U_buff_pre =U_buff;
    %---------- Relax ---------------------
    Lambda_com = Lambda(:,1:Nt) + 1i* Lambda(:,Nt+1:end) ;
    U_hat = U_buff - Lambda_com/pho;
    R_rec = relax_2D(U_hat,K);

    R_real = [real(R_rec), -imag(R_rec);
              imag(R_rec), real(R_rec)
        ];
    Lambda_real = [ Lambda(:,1:Nt), -Lambda(:,Nt+1:end);
                    Lambda(:,Nt+1:end), Lambda(:,1:Nt)];

    lanm_vec = Lambda_real(:);
    r_vec = R_real(:);

    %----------- EM for updating u -----------------
    y_u = newEM_u(D,sigma, u_ini,lanm_vec,r_vec, pho,b_up,b_low);
    u_ini = y_u;

    U=reshape(y_u,[Nr*2,Nt*2]);
    U_buff = U(1:Nr,1:Nt)+ 1i* U(Nr+1:end,1:Nt);

    %------------- update Lambda  -----------------
    Lambda_real = Lambda_real + pho* (R_real- U);
    Lambda = [Lambda_real(1:Nr,1:Nt), Lambda_real(Nr+1:end, 1:Nt)];

    if norm(U_buff_pre -U_buff,'fro')<1e-3
        break
    end
end

H_admm= scaler_for_quantization.inverse_scale(R_rec);
end