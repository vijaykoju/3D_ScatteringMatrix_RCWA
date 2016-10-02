%RCWA
close all; clc; clear all;
% UNITS
meters = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
degrees = pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE PARAMETERS
lam0 = 2 * centimeters; %free space wavelength
%ginc = [0;0;1]; %unit vector in direction of incident wave
theta = [0:0.1:89] * degrees;             % elevation angle
NN = length(theta);
phi   = 0 * degrees;             % azimuthal angle
pte   = 1;%1/sqrt(2);                       % amplitude of TE polarization
ptm   = 0;%1j/sqrt(2);                       % amplitude of TM polarization

% DEVICE PARAMETERS
ur1 = 1.0; %permeability in reflection region
er1 = 2.0; %permittivity in reflection region
ur2 = 1.0; %permeability in transmission region
er2 = 9.0; %permittivity in transmission region
urd = 1.0; %permeability of device
erd = 6.0; %permittivity of device
Lx = 1.75 * centimeters; %period along x
Ly = 1.5 * centimeters; %period along y
d1 = 0.5 * centimeters; %thickness of layer 1
d2 = 0.3 * centimeters; %thickness of layer 2
w = 0.8*Ly; %length of one side of triangle

% RCWA PARAMETERS
Nx = 512; %number of point along x in real-space grid
Ny = round(Nx*Ly/Lx); %number of point along y in real-space grid
PQ = 3 * [1 1]; %number of spatial harmonics along x and y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: BUILD DEVICE LAYERS ON HIGH RESOLUTION GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROSS SECTIONAL GRID
NL = 2;
dx = Lx/Nx; %grid resolution along x
dy = Ly/Ny; %grid resolution along y
xa = [0:Nx-1]*dx; %x axis array
xa = xa - mean(xa); %center x axis at zero
ya = [0:Ny-1]*dy; %y axis vector
ya = ya - mean(ya); %center y axis at zero
% INITIALIZE LAYERS TO er AND ur
UR = urd * ones(Nx,Ny,NL);
ER = erd * ones(Nx,Ny,NL);
L = [ d1 d2 ];

% BUILD LAYER 1 (TRIANGLE)
h = 0.5*sqrte(3)*w;
ny = round(h/dy);
ny1 = round((Ny - ny)/2);
ny2 = ny1 + ny - 1;
for ny = ny1 : ny2
    f = (ny - ny1)/(ny2 - ny1);
    nx = round(f*w/Lx*Nx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx;
    ER(nx1:nx2,ny,1) = er1;
end
figure(1)
pcolor(1:Ny,1:Nx,ER(:,:,1));shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: COMPUTE CONVOLUTION MATRICES OF EACH LAYER OF DEVICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
URC = zeros(PQ(1)^2,PQ(2)^2,NL);
ERC = zeros(PQ(1)^2,PQ(2)^2,NL);
for i = 1:NL
    URC(:,:,i) = convmat(UR(:,:,i),PQ(1),PQ(2));
    ERC(:,:,i) = convmat(ER(:,:,i),PQ(1),PQ(2));
end


REF(1:NN) = 0;
TRN(1:NN) = 0;

parfor ii=1:NN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 4: COMPUTE WAVE VECTOR EXPANSIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I = eye(PQ(1)^2,PQ(2)^2);
    Z = zeros(PQ(1)^2,PQ(2)^2);

    n1 = sqrte(er1);
    n2 = sqrte(er2);
    k0 = 2*pi/lam0;
    kinc = n1*[sin(theta(ii))*cos(phi); sin(theta(ii))*sin(phi); cos(theta(ii))];


    p = [-floor(PQ(1)/2):+floor(PQ(1)/2)];         % indices along x
    q = [-floor(PQ(2)/2):+floor(PQ(2)/2)];         % indices along y

    Kx = kinc(1) - 2*pi*p/(k0*Lx);
    Ky = kinc(2) - 2*pi*q/(k0*Ly);

    [Ky, Kx] = meshgrid(Ky,Kx);

    Kx = diag(sparse(Kx(:)));
    Ky = diag(sparse(Ky(:)));

    Kzr = -conj(sqrte(ur1*er1*I-Kx^2-Ky^2));
    Kzt = conj(sqrte(ur2*er2*I-Kx^2-Ky^2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 5: COMPUTE EIGEN MODES OF FREE SPACE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Kz = conj(sqrte(I));
    Q = [Kx*Ky I+Ky^2; -(I+Kx^2) -Kx*Ky];
    W0 = [I Z; Z I];
    V0 = -1j*Q;
%     if theta(ii)== 0
%         Kz = conj(sqrt(I-Kx^2-Ky^2));
%         Q = [Kx*Ky I-Kx^2; Ky^2-I -Kx*Ky];
%         W0 = [I Z; Z I];
%         LAM = [1j*Kz Z; Z 1j*Kz];
%         V0 = Q/LAM;
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 6: INITIALIZE GLOBAL SCATTERING MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xx1 yy1] = size(V0);
    Sg.s11 = zeros(xx1,yy1);
    Sg.s12 = eye(xx1,yy1);
    Sg.s21 = eye(xx1,yy1);
    Sg.s22 = zeros(xx1,yy1);
    II = eye(xx1,yy1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 7: MAIN LOOP FOR S-MATRICES THROUGH LAYERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:NL
        S = scattering_matrix_ith(ERC(:,:,i),URC(:,:,i),Kx,Ky,k0,L(i),W0,V0);
        Sg = Redheffer_star_product(Sg,S,II);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 8: COMPUTE RELECTION SIDE CONNECTION S-MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Sr, Wref] = scattering_matrix_ref(er1,ur1,Kx,Ky,Kzr,Z,I,II,W0,V0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 9: COMPUTE TRANSMISSION SIDE CONNECTION S-MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [St, Wtrn] = scattering_matrix_trn(er2,ur2,Kx,Ky,Kzt,Z,I,II,W0,V0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 10: UPDATE GLOBAL SCATTERING MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sgg = Redheffer_star_product(Sr,Sg,II);
    Sgg = Redheffer_star_product(Sgg,St,II);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 11: COMPUTE REFLECTION AND TRANSMITTED FIELDS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta = zeros(PQ(1)^2,1);
    delta(ceil(length(delta)/2),1) = 1;

    n_hat = [0; 0; -1];
    %%% cross product of k_inc and n_hat : k_inc x n_hat
    [k_inc_cross_n_hat, mag_k_inc_cross_n_hat] = vector_cross_product(kinc/k0,n_hat);

    %%% direction of TE polarized light
    if (theta(ii) == 0)
        a_hat_te = [0; 1; 0];
    else
        a_hat_te = [k_inc_cross_n_hat.x/mag_k_inc_cross_n_hat;
                    k_inc_cross_n_hat.y/mag_k_inc_cross_n_hat;
                    k_inc_cross_n_hat.z/mag_k_inc_cross_n_hat];
    end

    %%% cross product of a_hat_te and k_inc : a_hat_te x k_inc
    [a_hat_te_cross_k_inc, mag_a_hat_te_cross_k_inc] = vector_cross_product(a_hat_te,kinc/k0);

    %%% direction of TM polarized light
    a_hat_tm = [a_hat_te_cross_k_inc.x/mag_a_hat_te_cross_k_inc;
                a_hat_te_cross_k_inc.y/mag_a_hat_te_cross_k_inc;
                a_hat_te_cross_k_inc.z/mag_a_hat_te_cross_k_inc];

    %%% Polarization vector
    EP = pte*a_hat_te + ptm*a_hat_tm;

    esrc = zeros(2*PQ(1)^2,1);
    esrc(1:length(esrc)/2,1)=EP(1)*delta;
    esrc(length(esrc)/2+1:end,1)=EP(2)*delta;

    csrc = Wref\esrc;

    cref = Sgg.s11*csrc;
    ctrn = Sgg.s21*csrc;

    eref = Wref*cref;
    etrn = Wtrn*ctrn;

    rx = eref(1:length(eref)/2);
    ry = eref(length(eref)/2+1:end);
    rz = -Kzr\(Kx*rx+Ky*ry);

    tx = etrn(1:length(etrn)/2);
    ty = etrn(length(etrn)/2+1:end);
    tz = -Kzt\(Kx*tx+Ky*ty);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 12: COMPUTE DIFFRACTION EFFICIENCIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r2 = abs(rx).^2+abs(ry).^2+abs(rz).^2;
    R = real(-Kzr/kinc(3))*r2;
    REF(ii) = sum(R);

    t2 = abs(tx).^2+abs(ty).^2+abs(tz).^2;
    T = real((ur1/ur2)*Kzt/kinc(3))*t2;
    TRN(ii) = sum(T);
end
figure(2)
plot(theta/degrees,REF,theta/degrees,TRN);
xlabel('Incident Angle [degree]')
ylabel('Reflectivity')