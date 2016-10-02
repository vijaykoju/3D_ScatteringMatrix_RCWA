function S = Redheffer_star_product(SA,SB,I)
    % STAR Redheffer Star Product
    %
    % S = Redheffer_star_prodcut(SA,SB)
    %
    % INPUT ARGUMENTS
    % ================
    % SA First Scattering Matrix
    % .s11 S11 scattering parameter
    % .s12 S12 scattering parameter
    % .s21 S21 scattering parameter
    % .s22 S22 scattering parameter
    %
    % SB Second Scattering Matrix
    % .s11 S11 scattering parameter
    % .s12 S12 scattering parameter
    % .s21 S21 scattering parameter
    % .s22 S22 scattering parameter
    %
    % OUTPUT ARGUMENTS
    % ================
    % S Combined Scattering Matrix
    % .s11 S11 scattering parameter
    % .s12 S12 scattering parameter
    % .s21 S21 scattering parameter
    % .s22 S22 scattering parameter
    %
    S.s11    = SA.s11 + SA.s12*((I-SB.s11*SA.s22)\(SB.s11*SA.s21));
    S.s12    = SA.s12*((I-SB.s11*SA.s22)\SB.s12);
    S.s21    = SB.s21*((I-SA.s22*SB.s11)\SA.s21);
    S.s22    = SB.s22 + SB.s21*((I-SA.s22*SB.s11)\(SA.s22*SB.s12));
end