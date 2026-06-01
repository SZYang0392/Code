function str = fortran_ew_d(x, w, d, leading_zero, CapitalizeE)
    % Fortran Ew.d formatted numeric output, with control over leading zero before decimal point
    % Inputs:
    %   x            : numeric value to be formatted
    %   w            : total field width (all characters included, right justified)
    %   d            : number of decimal digits in the mantissa
    %   leading_zero : logical, true -> leading zero before decimal point (0.ddd...), false -> one non-zero digit before decimal point (d.ddd...)
    % Output:
    %   str          : string of length w, positive numbers start with a space, negative numbers start with '-'
    %                   exponent fixed to two digits (with sign), zero-padded

    % Default leading_zero = true
    if nargin < 4
        leading_zero = true;
    end

    % Default CapitalizeE = true
    if nargin < 5
        CapitalizeE = true;
    end

    % Handle sign and absolute value
    if x < 0
        sign_char = '-';
        x_abs = -x;
    else
        sign_char = ' ';
        x_abs = x;
    end
    
    % Special case for zero
    if x_abs == 0
        exponent = 0;
        mantissa = 0;
    else
        if leading_zero
            % Normalize to [0.1, 1): exponent = floor(log10(x_abs)) + 1
            exponent = floor(log10(x_abs)) + 1;
            mantissa = x_abs / (10^exponent);
            % Guard against floating-point errors
            if mantissa >= 1
                mantissa = mantissa / 10;
                exponent = exponent + 1;
            elseif mantissa < 0.1
                mantissa = mantissa * 10;
                exponent = exponent - 1;
            end
        else
            % Normalize to [1, 10): exponent = floor(log10(x_abs))
            exponent = floor(log10(x_abs));
            mantissa = x_abs / (10^exponent);
            % Guard against floating-point errors
            if mantissa >= 10
                mantissa = mantissa / 10;
                exponent = exponent + 1;
            elseif mantissa < 1
                mantissa = mantissa * 10;
                exponent = exponent - 1;
            end
        end
    end
    
    % Construct mantissa string: always use %.*f to get the correct representation
    % For leading_zero=true, mantissa is in [0.1,1) so sprintf produces '0.ddd...'
    % For leading_zero=false, mantissa is in [1,10) so sprintf produces 'd.ddd...'
    mantissa_str = sprintf('%.*f', d, mantissa);
    
    % Construct exponent part: 'E' followed by sign and two digits (zero-padded)
    if CapitalizeE
        exp_str = sprintf('E%+03d', exponent);
    else
        exp_str = sprintf('e%+03d', exponent);
    end
    
    % Base string (sign + mantissa + exponent)
    base_str = [sign_char, mantissa_str, exp_str];
    base_len = length(base_str);   % base length = 1(sign) + (d+2)(mantissa) + 4(exponent) = d+7
    
    % Check if field width is sufficient
    if base_len > w
        error('Width w=%d is insufficient, at least %d characters needed', w, base_len);
    end
    
    % Right justify: pad left with spaces
    padding = repmat(' ', 1, w - base_len);
    str = [padding, base_str];
end