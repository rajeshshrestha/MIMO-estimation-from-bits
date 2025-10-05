classdef FileHGenerator < handle
    properties
        K {mustBeNumeric,mustBePositive}
        Nr {mustBeNumeric,mustBePositive}
        Nt {mustBeNumeric,mustBePositive}
        HSNR {mustBeFloat}
        file_path

        delta_sd
        lim
        sp_r,sp_t
        HSNR_mag
        H_vals
        total_num_samples
        actual_K

        % Variables to randomly select the H vlaues
        current_index
        fetch_indices



    end
    methods
        function obj = FileHGenerator(K, Nr, Nt, HSNR, file_path)
            % Assign values
            obj.K = K;
            obj.Nr = Nr;
            obj.Nt = Nt;
            obj.HSNR=HSNR;
            obj.file_path = file_path;
            obj.HSNR_mag=10.^(obj.HSNR/10);

            % Load data from file
            data = load(file_path);
            obj.H_vals = data.user_channel_data;
            obj.actual_K = data.actual_num_paths;
            obj.total_num_samples = numel(obj.H_vals);

            rng("shuffle");

            % Randomize the indices
            obj.fetch_indices = randperm(size(obj.H_vals,1));

            % current index for fetching H
            obj.current_index = 1;

            % asser match of Nr and Nt
            assert(size(obj.H_vals,2) == obj.Nr, 'Nr does not match the file data');
            assert(size(obj.H_vals,3) == obj.Nt, 'Nt does not match the file data');


           
        end

        function H = getnextH(obj)
            % Fetch current H
            H= 0;
            while(norm(H, 'fro')==0)
                H = squeeze(obj.H_vals(obj.fetch_indices(obj.current_index),:,:));
                obj.current_index = obj.current_index + 1;
            end
            if obj.current_index > obj.total_num_samples
                obj.current_index = 1;
                obj.fetch_indices = randperm(size(obj.H_vals,1));
            end
        end

        function H_noisy = addNoise(obj, H_true)
            % Add noise
            data_sigma = norm(H_true,'fro')/sqrt(obj.HSNR_mag * numel(H_true));
            Re_H = real(H_true) + data_sigma * randn(size(H_true));
            Im_H = imag(H_true) + data_sigma * randn(size(H_true));
            H_noisy = Re_H + 1i * Im_H;
        end

        function [H_true,H_noisy,h] = getH(obj)
            H_true = obj.getnextH();
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