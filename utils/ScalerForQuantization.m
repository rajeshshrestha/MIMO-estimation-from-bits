classdef ScalerForQuantization < handle

    properties
        scaling_type
        scaling_fn
        inverse_scaling_fn

        nbits
        max_abs_val
        scale_ratio
        mu
        s
    end

    methods

        function obj = ScalerForQuantization(nbits, mu)

            obj.nbits = nbits;

            if nargin < 2
                mu = 255;
            end

            obj.mu = mu;

            % Get the scaling type from the environment variable
            obj.scaling_type = getenv("QuantizationScaleType");

            if isempty(obj.scaling_type)
                obj.scaling_type = "linear";
            end

        end

        function scaled_x = forward_scale(obj, x)
            % Set the max value
            obj.max_abs_val = max(abs(x(:)));

            % Set the scaling function based on the scaling type
            if (obj.scaling_type == "mulaw")
                scaled_x = obj.mulaw_forward(x);
            elseif (obj.scaling_type == "linear")
                scaled_x = x;
            else
                throw(MException("Invalid scale type provided"))
            end

            % Scale the signal so that the quantized range is [0, 2^nbits-1]
            obj.scale_ratio = 2 ^ obj.nbits / obj.max_abs_val;

            scaled_x = scaled_x * obj.scale_ratio;

        end

        function orig_scale_x = inverse_scale(obj, x)

            % Unscale with the scale ratio
            x = x / obj.scale_ratio;

            % Set the scaling function based on the scaling type
            if (obj.scaling_type == "mulaw")
                orig_scale_x = obj.mulaw_inverse(x);
            elseif (obj.scaling_type == "linear")
                orig_scale_x = x;
            else
                throw(MException("Invalid scale type provided"))
            end

        end

        function x = mulaw_forward(obj, y)

            y_normalized = y / obj.max_abs_val;
            % Apply the mu-law scaling
            x = sign(y_normalized) .* log(1 + obj.mu .* abs(y_normalized)) ./ log(1 + obj.mu);

            % Scale the signal to the original range
            x = x * obj.max_abs_val;
        end

        function y = mulaw_inverse(obj, x)
            % Normalize the signal
            x_normalized = x / obj.max_abs_val;

            % Apply the inverse mu-law scaling
            y = sign(x_normalized) .* (1 ./ obj.mu) .* ((1 + obj.mu) .^ abs(x_normalized) - 1);
            % Scale the signal to the original range
            y = y * obj.max_abs_val;
        end

    end

end
