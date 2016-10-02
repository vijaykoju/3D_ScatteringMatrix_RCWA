function C = convmat(A,P,Q,R)
% CONVMAT    Rectangular Convolution Matrix
%
% C = convmat(A,P)                 for 1D problems
% C = convmat(A,P,Q)               for 2D problems
% C = convmat(A,P,Q,R)             for 3D problems
%
% This MATLAB function constructs convolution matrices from a real space
% grid.
%

%%% Handle input and output arguments

% Determine size of A
[Nx, Ny, Nz] = size(A);

% Handle number of harmonics for all dimensions
if nargin==2
    Q = 1;
    R = 1;
elseif nargin==3
    R = 1;
end

% Compute indices of spatial harmonics
NH = P*Q*R;                            % total number
p = [-floor(P/2):+floor(P/2)];        % indices along x
q = [-floor(Q/2):+floor(Q/2)];        % indices along y
r = [-floor(R/2):+floor(R/2)];         % indices along z

% Compute Fourier coefficients of A
A = fftshift(fftn(A))/(Nx*Ny*Nz);

% Compute array indices of center harmonic
p0 = 1 + floor(Nx/2);
q0 = 1 + floor(Ny/2);
r0 = 1 + floor(Nz/2);

for rrow = 1:R
for qrow = 1:Q
for prow = 1:P
    row = (rrow-1)*Q*P + (qrow-1)*P + prow;
    for rcol = 1:R
    for qcol = 1:Q
    for pcol = 1:P
        col = (rcol-1)*Q*P + (qcol-1)*P + pcol;
        pfft = p(prow) - p(pcol);
        qfft = q(qrow) - q(qcol);
        rfft = r(rrow) - r(rcol);
        C(row,col) = A(p0+pfft,q0+qfft,r0+rfft);
        %disp([p0+pfft q0+qfft r0+rfft]);
    end
    end
    end
end
end
end

end