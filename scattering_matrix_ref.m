function [Sr, Wref] = scattering_matrix_ref(er1,ur1,Kx,Ky,Kzr,Z,I,II,W0,V0)
    %%% function to compute scattering matrix of reflection region.
    Q = (1/ur1)*[Kx*Ky ur1*er1*I-Kx^2; Ky^2-ur1*er1*I -Ky*Kx];
    Wref = [I Z; Z I];
    LAM = [-1j*Kzr Z; Z -1j*Kzr];
    Vref = Q/LAM;
    
    aa = W0\Wref;
    bb = V0\Vref;

    A = aa + bb;
    B = aa - bb;

    Sr.s11 = - A\B;
    Sr.s12 = (2*II)/A;
    Sr.s21 = 0.5*(A-B*(A\B));
    Sr.s22 = B/A;
    
end