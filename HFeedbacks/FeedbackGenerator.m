classdef FeedbackGenerator

    properties
        SNR {mustBeNumeric}
        YSNR {mustBeNumeric}
        YSNR_mag {mustBeNumeric}
        bits {mustBeNumeric, mustBePositive}
        BER {mustBeNumeric, mustBeNonnegative, mustBeLessThanOrEqual(BER, 1)}

        SNR_mag
        avoid_dither_noise
        scaler_for_quantization
    end

    methods

        function obj = FeedbackGenerator(SNR, bits, BER, YSNR)

            % Initialize properties
            obj.SNR = SNR;
            obj.bits = bits;
            obj.BER = BER;
            obj.YSNR = YSNR;

            obj.SNR_mag = 10 .^ (obj.SNR / 10);
            obj.YSNR_mag = 10 .^ (obj.YSNR / 10);

            obj.avoid_dither_noise = str2double(getenv("AvoidDitherNoiseAddition"));

            if isnan(obj.avoid_dither_noise)
                obj.avoid_dither_noise = 0;
            elseif (obj.avoid_dither_noise == 0)
                % fprintf("Using normal dithering.\n");
            else
                obj.avoid_dither_noise = 1;
                fprintf("Avoid adding dithering noise.\n");

            end

            obj.scaler_for_quantization = ScalerForQuantization(obj.bits);

        end

        function [z, D, sigma_true] = newCompressQuantize(obj, H_true, nPilots)
            [Nr, Nt] = size(H_true);

            % Pilot Signal
            S = randn(Nt, nPilots);

            % Response of the channel
            Y = H_true * S;

            if ~isnan(obj.YSNR)
                sigma_noise = norm(Y, 'fro') / sqrt(obj.YSNR_mag * numel(Y));
                v = sigma_noise / sqrt(2) * randn(size(Y)) + 1i * sigma_noise / sqrt(2) * randn(size(Y));
                Y = Y + v;
            end

            % Compression
            X_tilde = [real(S); imag(S)];
            Y_tilde = [real(Y); imag(Y)];
            A = randn(size(Y_tilde));
            A = A / sqrt(size(A, 2));
            D = kr(X_tilde, A).';
            actual_compression = diag(A.' * Y_tilde);
            q = actual_compression;

            % Dithering
            sigma_true = norm(q) / sqrt(obj.SNR_mag * length(q));
            v = sigma_true * randn(length(q), 1);

            if (obj.avoid_dither_noise == 1)
                sig = q + v;
            else
                sig = q;
            end

            % Quantization
            partition = -2 ^ obj.bits + 2:2:2 ^ obj.bits - 2;
            codebook = -2 ^ obj.bits + 1:2:2 ^ obj.bits - 1;
            sig = obj.scaler_for_quantization.forward_scale(sig);

            [index, z] = myquantiz(sig, partition, codebook);

            %Introduce Bit error
            if obj.BER ~= 0
                disp('Introducing Bit Error')
                idxs = 1:length(z);
                idxs = idxs(randperm(length(idxs)));
                idxs = idxs(1:round(obj.BER * length(z)));
                z(idxs) = datasample(codebook, length(idxs));
            end

            z = z';
        end

        function [z, A, sigma_true] = compressQuantize(obj, h, R)
            % Compression
            h_len = numel(h);
            A = randn(R, h_len) / sqrt(R);
            q = A * h;

            % Dithering
            sigma_true = norm(q) / sqrt(obj.SNR_mag * length(q));
            v = sigma_true * randn(length(q), 1);

            if (obj.avoid_dither_noise == 1)
                sig = q + v;
            else
                sig = q;
            end

            

            % Quantization
            partition = -2 ^ obj.bits + 2:2:2 ^ obj.bits - 2;
            codebook = -2 ^ obj.bits + 1:2:2 ^ obj.bits - 1;
            sig = obj.scaler_for_quantization.forward_scale(sig);
            [index, z] = myquantiz(sig, partition, codebook);

            %Introduce Bit error
            if obj.BER ~= 0
                disp('Introducing Bit Error')
                idxs = 1:length(z);
                idxs = idxs(randperm(length(idxs)));
                idxs = idxs(1:round(obj.BER * length(z)));
                z(idxs) = datasample(codebook, length(idxs));
            end

            z = z';
        end

        function [H_estimated] = estimateH_OVSF(obj, H_true, Npilots)
            [Nr, Nt] = size(H_true);
            S = generateOVSF(Nt, Npilots);

            % Response of the channel
            Y = H_true * S;

            if ~isnan(obj.YSNR)
                sigma_noise = norm(Y, 'fro') / sqrt(obj.YSNR_mag * numel(Y));
                v = sigma_noise / sqrt(2) * randn(size(Y)) + 1i * sigma_noise / sqrt(2) * randn(size(Y));
                Y = Y + v;
            end

            H_estimated = Y * S';

        end

        function [z, A, sigma_true] = feedbackADMM(obj, H_true, R, Npilotsfactor)
            [Nr, Nt] = size(H_true);
            Npilots = 2 ^ (round(log2(Npilotsfactor * Nt)));

            % Estimate H using OVSF codes
            H_estimated = obj.estimateH_OVSF(H_true, Npilots);

            % Vectorize the channel
            h = obj.getVech(H_estimated);

            % Compress and Quantize the channel
            [z, A, sigma_true] = obj.compressQuantize(h, R);

        end

        function [z, A, sigma_true] = feedbackNewADMM(obj, H_true, R)
            [z, A, sigma_true] = obj.newCompressQuantize(H_true, R);
        end

    end

    methods (Static)

        function h = getVech(H)
            Re_H = real(H);
            Im_H = imag(H);
            U_true = [Re_H, Im_H];
            h = U_true(:);
        end

    end

end
