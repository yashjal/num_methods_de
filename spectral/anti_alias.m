function w_hat = anti_alias(u_hat, v_hat)
% anti-aliasing via oversampling
% computes the anti-aliased Fourier coeffiecients of w=uv
% given Fourier coefficients of u,v
% ref: https://cims.nyu.edu/~donev/Teaching/NMII/Lectures/Aliasing.pdf

    N = size(u_hat,2);
    M = 2*N;
    mid = floor(N/2);
    u_pad = ifft([u_hat(1:mid) zeros(1, M-N) u_hat((mid+1):N)]);
    v_pad = ifft([v_hat(1:mid) zeros(1, M-N) v_hat((mid+1):N)]);
    w_pad = u_pad .* v_pad;
    w_pad_hat = fft(w_pad);
    if mod(N,2) == 0
        w_hat = 2*[w_pad_hat(1:mid) w_pad_hat(M-mid+1:M)];
    else
        w_hat = 2*[w_pad_hat(1:mid+1) w_pad_hat(M-mid+1:M)];
    end
end