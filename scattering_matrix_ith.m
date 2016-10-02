function S = scattering_matrix_ith(ERC,URC,Kx,Ky,k0,L,W0,V0)
    %%% function to compute scattering matrix of ith layer of the device.
    P = [Kx*(ERC\Ky) URC-Kx*(ERC\Kx); Ky*(ERC\Ky)-URC -Ky*(ERC\Kx)];
    Q = [Kx*(URC\Ky) ERC-Kx*(URC\Kx); Ky*(URC\Ky)-ERC -Ky*(URC\Kx)];
    Omega2 = P*Q;
    [W,LAM] = eig(Omega2);
    LAM = sqrt(LAM);
    V = Q*W/LAM;
    
    aa = W\W0;
    bb = V\V0;
    
    A = aa + bb;
    B = aa - bb;
    X = expm(-LAM*k0*L);

    % Components of Scattering matrix for ith layer
    S.s11    = (A-X*B*(A\(X*B)))\(X*B*(A\(X*A))-B);
    S.s12    = ((A-X*B*(A\(X*B)))\X)*(A-B*(A\B));
    S.s22    = S.s11;
    S.s21    = S.s12;
end