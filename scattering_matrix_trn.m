function [St, Wtrn] = scattering_matrix_trn(er2,ur2,Kx,Ky,Kzt,Z,I,II,W0,V0)
    %%% function to compute scattering matrix of transmission region.
    Q = (1/ur2)*[Kx*Ky ur2*er2*I-Kx^2; Ky^2-ur2*er2*I -Ky*Kx];
    Wtrn = [I Z; Z I];
    LAM = [1j*Kzt Z; Z 1j*Kzt];
    Vtrn = Q/LAM;

    aa = W0\Wtrn;
    bb = V0\Vtrn;

    A = aa + bb;
    B = aa - bb;

    St.s11 = B/A;
    St.s12 = 0.5*(A-B*(A\B));
    St.s21 = (2*II)/A;
    St.s22 = -A\B;
end