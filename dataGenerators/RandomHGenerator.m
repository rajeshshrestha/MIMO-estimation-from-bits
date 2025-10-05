classdef RandomHGenerator
    properties
        K {mustBeNumeric,mustBePositive}
        Nr {mustBeNumeric,mustBePositive}
        Nt {mustBeNumeric,mustBePositive}
        alpha {mustBeFloat, mustBePositive}
        HSNR {mustBeFloat}

        delta_sd
        lim
        sp_r,sp_t
        HSNR_mag


    end
    methods
        function obj = RandomHGenerator(K, Nr, Nt, alpha, HSNR)
            % Assign values
            obj.K = K;
            obj.Nr = Nr;
            obj.Nt = Nt;
            obj.alpha = alpha;
            obj.HSNR=HSNR;

            % Initialize required variables
            obj.delta_sd = 4*pi/180;
            angle_low = (-pi/2)+obj.delta_sd;
            angle_high = (pi/2)-obj.delta_sd;
            obj.lim=[angle_low,angle_high];
            obj.sp_r = 0:1:obj.Nr-1;
            obj.sp_t = 0:1:obj.Nt-1;
            obj.HSNR_mag=10.^(obj.HSNR/10);


            % Add necessary files to path
            mfilepath = fileparts(which(mfilename));
        end

        function H = generateH(obj)
            % Generate angles
            phi=angle_generation(obj.K, obj.lim, obj.delta_sd);
            theta=angle_generation(obj.K, obj.lim, obj.delta_sd);

            
            beta_low = str2double(getenv("betaLow"));
            beta_high = str2double(getenv("betaHigh"));
            
            if  isnan(beta_low)
                beta_low = 0.5;
            else
                % beta_low
            end
            if  isnan(beta_high)
                beta_high = 1;
            else
                % beta_high
             end


            beta_high = 1;
            beta = (rand(obj.K,1)*(beta_high-beta_low)+beta_low).*exp(1i*rand(obj.K,1)*2*pi)/sqrt(obj.K);


            % Compute the channel H
            theta_true = -2*pi*obj.alpha* sin(theta);
            phi_true = - 2*pi*obj.alpha*  sin(phi);
            Ar = exp(1i*theta_true* obj.sp_r).';
            At = exp(1i*phi_true* obj.sp_t).';
            H = Ar*diag(beta)*At.';
        end

        function H = generateAngularH(obj)
            Ar = dftmtx(obj.Nr);
            At = dftmtx(obj.Nt);

            beta = (rand(obj.K,1)/2+0.5).*exp(1i*rand(obj.K,1)*2*pi)/sqrt(obj.K);
            H_a = zeros(obj.Nr,obj.Nt);
            for k=1:obj.K
                i = randi([1, obj.Nr],1);
                j = randi([1,obj.Nt],1);
                H_a(i,j)=beta(k);
            end
            H = Ar*H_a*At';
        end

        function H_noisy = addNoise(obj, H_true)
            % Add noise
            data_sigma = norm(H_true,'fro')/sqrt(obj.HSNR_mag * numel(H_true));
            Re_H = real(H_true) + data_sigma * randn(size(H_true));
            Im_H = imag(H_true) + data_sigma * randn(size(H_true));
            H_noisy = Re_H + 1i * Im_H;
        end



        function [H_true,H_noisy,h] = getH(obj)
            H_true = obj.generateH();
            H_noisy = obj.addNoise(H_true);
            h = obj.getVech(H_noisy);
        end


    end
    methods(Static)
        function h=getVech(H)
            Re_H = real(H);
            Im_H = imag(H);
            U_true =[Re_H, Im_H];
            h = U_true(:);
        end
    end

end