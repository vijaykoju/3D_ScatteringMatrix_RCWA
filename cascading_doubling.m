function SC = cascading_doubling(S,II,ZZ,N)
    % CASC Cascading and Doubling Algorithm
    %
    % SC = cascading_doubling(S,N)
    %
    % This MATLAB function uses an efficient doubling algorithm
    % to cascade N periods.
    % 
    % INPUT ARGUMENTS
    % ================
    % S Scattering Matrix for One Period
    % .S11 S11 scattering parameter
    % .S12 S12 scattering parameter
    % .S21 S21 scattering parameter
    % .S22 S22 scattering parameter
    %
    % N Number of scattering matrices to cascade
    %
    % OUTPUT ARGUMENTS
    % ================
    % SC Overall Scattering Matrix for Cascade
    N_binary = dec2bin(N);
    one = dec2bin(1);
    
    SC.s11 = ZZ;
    SC.s12 = II;
    SC.s21 = II;
    SC.s22 = ZZ;
    
    Sbin = S;
    N_len = length(N_binary);
    for i=1:N_len
        if (N_binary(N_len-i+1) == one)
            SC = Redheffer_star_product(SC,Sbin,II);
        end
        Sbin = Redheffer_star_product(Sbin,Sbin,II);
    end
end