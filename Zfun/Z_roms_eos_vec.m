function [den,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,z_r)
% Z_roms_eos_vec.m  5/23/2014  Parker MacCready
%
% the equation of state used in ROMS, vectorized to accept and return
% arrays
%
% INPUT (all must be the same size):
% salt = practical salinity
% temp = potential temperature (deg C) [what ROMS uses]
% z_r = z (m) [used as an approximation of pressure] [typically NEGATIVE]
%
% OUTPUT:
% den = in situ density anomaly (kg m-3) [rho - 1000]
% den1 = density anomaly of this water at z=0 (kg m-3) [sigma0]
% alpha = thermal expansion (1/Celsius) coefficient
% beta = saline contraction (1/PSU) coefficient
%   Hence, using rho = den + 1000:
%   alpha = -(1/rho) drho/dtemp
%   beta =  (1/rho) drho/dsalt
%
% we include all code assuming BULK_FLUXES and EOS_TDERIVATIVE are defined

% parameters from ROMS./ROMS./Modules./mod_eoscoef.F

% !svn $Id: mod_eoscoef.F 585 2012-01-03 18:44:28Z arango $
% !================================================== Hernan G. Arango ===
% !  Copyright (c) 2002-2012 The ROMS./TOMS Group                         !
% !    Licensed under a MIT./X style license                              !
% !    See License_ROMS.txt                                              !
% !=======================================================================
% !                                                                      !
% !  Polynomial  expansion  coefficients for the computation of          !
% !  "in situ" density  and other associated quantities via the          !
% !  nonlinear equation of state for seawater  as a function of          !
% !  potential temperature, salinity, and pressure (Jackett and          !
% !  McDougall, 1995). [fixed date]                                      !
% !                                                                      !
% !=======================================================================

A00 = +1.909256e+04;
A01 = +2.098925e+02;
A02 = -3.041638e+00;
A03 = -1.852732e-03;
A04 = -1.361629e-05;
B00 = +1.044077e+02;
B01 = -6.500517e+00;
B02 = +1.553190e-01;
B03 = +2.326469e-04;
D00 = -5.587545e+00;
D01 = +7.390729e-01;
D02 = -1.909078e-02;
E00 = +4.721788e-01;
E01 = +1.028859e-02;
E02 = -2.512549e-04;
E03 = -5.939910e-07;
F00 = -1.571896e-02;
F01 = -2.598241e-04;
F02 = +7.267926e-06;
G00 = +2.042967e-03;
G01 = +1.045941e-05;
G02 = -5.782165e-10;
G03 = +1.296821e-07;
H00 = -2.595994e-07;
H01 = -1.248266e-09;
H02 = -3.508914e-09;
Q00 = +9.99842594e+02;
Q01 = +6.793952e-02;
Q02 = -9.095290e-03;
Q03 = +1.001685e-04;
Q04 = -1.120083e-06;
Q05 = +6.536332e-09;
U00 = +8.24493e-01;
U01 = -4.08990e-03;
U02 = +7.64380e-05;
U03 = -8.24670e-07;
U04 = +5.38750e-09;
V00 = -5.72466e-03;
V01 = +1.02270e-04;
V02 = -1.65460e-06;
W00 = +4.8314e-04;

% from ROMS./ROMS./Nonlinear./rho_eos.F

% !svn $Id: rho_eos.F 585 2012-01-03 18:44:28Z arango $
% !================================================== Hernan G. Arango ===
% !  Copyright (c) 2002-2012 The ROMS./TOMS Group                         !
% !    Licensed under a MIT./X style license                              !
% !    See License_ROMS.txt                                              !
% !=======================================================================
% !                                                                      !
% !  This routine computes  "in situ" density and other associated       !
% !  quantitites as a function of potential temperature,  salinity,      !
% !  and pressure from a polynomial expression (Jackett & McDougall,     !
% !  1992). The polynomial expression was found from fitting to 248      !
% !  values  in the  oceanographic  ranges of  salinity,  potential      !
% !  temperature,  and pressure.  It  assumes no pressure variation      !
% !  along geopotential surfaces, that is, depth (meters; negative)      !
% !  and pressure (dbar; assumed negative here) are interchangeable.     !
% !                                                                      !
% !  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
% !                                                                      !
% !     alpha = 2.1014611551470d-04 (1./Celsius)
% what I got    2.101461155147048e-04 [OK]
% !     beta  = 7.2575037309946d-04 (1./PSU)
% what I got    7.257503730994606e-04 [OK]
% !     gamma = 3.9684764511766d-06 (1./Pa)                               !
% !     den   = 1050.3639165364     (kg./m3)
% what I got    1050.363916536439774 [OK, after adding 1000]
% !     den1  = 1028.2845117925     (kg./m3)
% what I got    1028.284511792455 [OK, after adding 1000]
% !     sound = 1548.8815240223     (m./s)                                !
% !     bulk  = 23786.056026320     (Pa)                                 !
% !                                                                      !
% !  Reference:                                                          !
% !                                                                      !
% !  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
% !    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
% !    and Oceanic Techn., vol. 12, pp. 381-389.                         !
% !                                                                      !
% !=======================================================================

% check values used in "what I got" above
% temp = [3];
% salt = [35.5];
% z_r = [-5000];

% !
% !=======================================================================
% !  Nonlinear equation of state.  Notice that this equation of state
% !  is only valid for potential temperature range of -2C to 40C and
% !  a salinity range of 0 PSU to 42 PSU.
% !=======================================================================
% !
% !
% !  Check temperature and salinity lower values. Assign depth to the
% !  pressure.
% !
Tt=max(-2.0,temp);
Ts=max(0.0,salt);
sqrtTs=sqrt(Ts);
Tp=z_r;
Tpr10=0.1.*Tp;


% !-----------------------------------------------------------------------
% !  Compute density (kg./m3) at standard one atmosphere pressure.
% !-----------------------------------------------------------------------
% !
C_0=Q00+Tt.*(Q01+Tt.*(Q02+Tt.*(Q03+Tt.*(Q04+Tt.*Q05))));
C_1=U00+Tt.*(U01+Tt.*(U02+Tt.*(U03+Tt.*U04)));
C_2=V00+Tt.*(V01+Tt.*V02);

dCdT_0=Q01+Tt.*(2.0.*Q02+Tt.*(3.0.*Q03+Tt.*(4.0.*Q04+Tt.*5.0.*Q05)));
dCdT_1=U01+Tt.*(2.0.*U02+Tt.*(3.0.*U03+Tt.*4.0.*U04));
dCdT_2=V01+Tt.*2.0.*V02;

den1=C_0+Ts.*(C_1+sqrtTs.*C_2+Ts.*W00);

% !
% !  Compute d(den1)./d(S) and d(den1)./d(T) derivatives used in the
% !  computation of thermal expansion and saline contraction
% !  coefficients.
% !
Dden1DS=C_1+1.5.*C_2.*sqrtTs+2.0.*W00.*Ts;
Dden1DT=dCdT_0+Ts.*(dCdT_1+sqrtTs.*dCdT_2);

% !
% !-----------------------------------------------------------------------
% !  Compute secant bulk modulus.
% !-----------------------------------------------------------------------
% !
C_3=A00+Tt.*(A01+Tt.*(A02+Tt.*(A03+Tt.*A04)));
C_4=B00+Tt.*(B01+Tt.*(B02+Tt.*B03));
C_5=D00+Tt.*(D01+Tt.*D02);
C_6=E00+Tt.*(E01+Tt.*(E02+Tt.*E03));
C_7=F00+Tt.*(F01+Tt.*F02);
C_8=G01+Tt.*(G02+Tt.*G03);
C_9=H00+Tt.*(H01+Tt.*H02);

dCdT_3=A01+Tt.*(2.0.*A02+Tt.*(3.0.*A03+Tt.*4.0.*A04));
dCdT_4=B01+Tt.*(2.0.*B02+Tt.*3.0.*B03);
dCdT_5=D01+Tt.*2.0.*D02;
dCdT_6=E01+Tt.*(2.0.*E02+Tt.*3.0.*E03);
dCdT_7=F01+Tt.*2.0.*F02;
dCdT_8=G02+Tt.*2.0.*G03;
dCdT_9=H01+Tt.*2.0.*H02;

bulk0=C_3+Ts.*(C_4+sqrtTs.*C_5);
bulk1=C_6+Ts.*(C_7+sqrtTs.*G00);
bulk2=C_8+Ts.*C_9;
bulk =bulk0-Tp.*(bulk1-Tp.*bulk2);

% this section is because we #define BULK_FLUXES in the ptx_01.h file

% !
% !  Compute d(bulk)./d(S) and d(bulk)./d(T) derivatives used
% !  in the computation of thermal expansion and saline contraction
% !  coefficients.
% !
DbulkDS=C_4+sqrtTs.*1.5.*C_5-Tp.*(C_7+sqrtTs.*1.5.*G00-Tp.*C_9);
DbulkDT=dCdT_3+Ts.*(dCdT_4+sqrtTs.*dCdT_5)- ...
    Tp.*(dCdT_6+Ts.*dCdT_7-Tp.*(dCdT_8+Ts.*dCdT_9));

% !
% !-----------------------------------------------------------------------
% !  Compute local "in situ" density anomaly (kg./m3 - 1000).
% !-----------------------------------------------------------------------
% !
cff=1.0./(bulk+Tpr10);
den=den1.*bulk.*cff;

den=den-1000.0;

% !
% !-----------------------------------------------------------------------
% !  Compute thermal expansion (1./Celsius) and saline contraction
% !  (1./PSU) coefficients.
% !-----------------------------------------------------------------------
% !
Tpr10=0.1.*z_r;
% !
% !  Compute thermal expansion and saline contraction coefficients.
% !
cff=bulk+Tpr10;
cff1=Tpr10.*den1;
cff2=bulk.*cff;
wrk=(den+1000.0).*cff.*cff;
Tcof=-(DbulkDT.*cff1+Dden1DT.*cff2);
Scof= (DbulkDS.*cff1+Dden1DS.*cff2);
cff=1.0./wrk;
alpha=cff.*Tcof;
beta =cff.*Scof;

% redefine den1 to be like sigma0
den1 = den1 - 1000;

