function [B] = Bidiag_Francis_Step(B)

    m = size(B, 1);

     % Compute Givens rotation matrix
    r_m = Givens_rotation([ B(1, 1)^2 - B(m-1, m)^2 - B(m, m)^2
          B(1, 1) * B(1, 2)]);
    % Apply Givens rotation to B(1:2,1:2)
    B(1:2, 1:2) = B(1:2, 1:2) * r_m;

    for z=2:m

        i = z;
        j = z-1;

        t = B(i-1:i, j);
        r_m = Givens_rotation(t);
        
        B(i-1:i, j) = [t' * r_m(:, 1) 0];

        k = j+2;

        if z == m
            k = j+1;
        
        end

        B(i-1:i, j+1:k) = r_m' * B(i-1:i, j+1:k);

        if z ~= m
            
            i = z - 1;
            j = z + 1;
            t = B(i, j-1:j)';
            r_m = Givens_rotation(t);
            B(i+1:i+2, j-1:j) = B(i+1:i+2, j-1:j) * r_m;
            B(i, j-1:j)=[t' * r_m(:, 1), 0];

        end
    end

end

