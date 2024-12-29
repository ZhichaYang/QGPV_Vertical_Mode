# QGPV_Vertical_Mode
Caculate the QG Vertical Mode in Ocean

Caculate the quasigeostrophy vertical mode based on the QGPV Sturm-Liouville Problem:

d/dz(f^2/N^2 d(Phi)/dz)+Lambda^2 Phi=0 

f Coriolis frequency; N2 background stratification; Lambda2 eigenvalue; Phi eigenfunction  
Use the a shooting method with a fourth order Runge–Kutta step,
integrating down from the surface.
The eigenvalue is adjusted using Newton’s method

the QGPV_VerticalMode_Shooting.m is the main function,and it needs dynmodes.m to run

2024.12.29 11:49:43
