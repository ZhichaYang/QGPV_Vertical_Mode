function [Lambda,Phi]=QGPV_VerticalMode_Shooting(Nsq, depth, f0, mode_number)
% Caculate the quasigeostrophy vertical mode based on the QGPV Sturm-Liouville Problem:
% d/dz(f^2/N^2 d(Phi)/dz)+Lambda^2 Phi=0 
% f Coriolis frequency; N2 background stratification; Lambda2 eigenvalue; Phi eigenfunction  
% Use the a shooting method with a fourth order Runge–Kutta step,
% integrating down from the surface.
% The eigenvalue is adjusted using Newton’s method
%
% INPUT:
% Nsq: buoyancy frequency(N^2); 
% depth: the Nsq's depth;
% f0: Coriolis frequency;
% mode_number: the nth mode you want to caculate;
%
% OUTPUT:
% Lambda: the eigenvalue
% Phi: nx2 array, the first column is eigenfunction, the second column is
%      its derivatives
% 
% EXAMPLE:
% f0=2*(2*pi/86400)*sind(19.5);
% [Lambda1,pmode1]=QGPV_VerticalMode_Shooting(n2,depth_t(1:end-1),f0,1);
%
% Author: Zhichao Yang 
% send comments and corrections to zhchyang@hhu.edu.cn
% 2024.12.29 11:21:39

% start

[~,~,ce]=dynmodes(Nsq,-depth);

delta=0.5*10^-8;  %Error judgment criteria
surface=1;        %Upper boundary point(Given in advance)
surface_dz=0;               %Upper boundary derivative value
bottom_dz=0;                %Low boundary derivative value
Lambda_ini=f0.^2/ce(mode_number).^2;
Lambda1= Lambda_ini-Lambda_ini/10;               %Eigenvalue guess
Lambda2= Lambda_ini+Lambda_ini/10;               %Eigenvalue guess



Phi1_ini = [surface,surface_dz]; 
Phi1_res = myrungekutta(Nsq, depth, f0, Lambda1, Phi1_ini);
phi_dz_rk1 = Phi1_res(:,2); 
F1=phi_dz_rk1(end)-bottom_dz; %get the error with the low boundary derivative value using fourth order Runge–Kutta

Phi2_ini = [surface,surface_dz]; 
Phi2_res = myrungekutta(Nsq, depth, f0, Lambda2, Phi2_ini);
Phi_dz_rk2 = Phi2_res(:,2); 
F2=Phi_dz_rk2(end)-bottom_dz; %get the error with the low boundary derivative value using fourth order Runge–Kutta

% start to iteration
D=min(abs(F1),abs(F2)); %erroe
nn=0;
while(D>delta)
    nn=nn+1;
%----------adjusted eigenvalue using Newton’s method-----------------

    Lambda3=Lambda2-F2*(Lambda2-Lambda1)/(F2-F1); 
    Phi3_ini = [surface,surface_dz]; 
    Phi3_res = myrungekutta(Nsq, depth, f0, Lambda3, Phi3_ini);  
    Phi_dz_rk3 = Phi3_res(:,2); 
    F3=Phi_dz_rk3(end)-bottom_dz; 
    

    Lambda1=Lambda2;
    F1=F2;
    Phi1_res=Phi2_res;

    Lambda2=Lambda3;
    F2=F3;
    Phi2_res=Phi3_res;
    
    D=min(abs(F1),abs(F2));

if nn>1000
    fprintf('NO!')
    break
end

end



if abs(F1)>abs(F2)
   Phi=Phi1_res;
   Lambda=Lambda1;
else
   Phi=Phi2_res;
   Lambda=Lambda2;
end


end
%%


function res = myrungekutta(Nsq, depth, f0, Lambda, z0)
% numerical integration using fourth order Runge–Kutta method

if depth(1)>0
    depth=[0;depth];
    Nsq=[Nsq(1);Nsq];
end

n = length(depth);

res(1,:)=z0;
for i = 1:n-1
    h = depth(i+1)-depth(i);
    dN_dz1=(Nsq(i)-Nsq(i+1))/(depth(i)-depth(i+1));
    if i<n-1
        dN_dz2=(Nsq(i+1)-Nsq(i+2))/(depth(i+1)-depth(i+2));
    else
        dN_dz2=dN_dz1;
    end

    dN_dz_h2=(dN_dz1+dN_dz2)/2;
    Nsq_h2=(Nsq(i)+Nsq(i+1))/2;

    L1 = res(i,2)/Nsq(i)*dN_dz1-Nsq(i)./f0^2*Lambda*res(i,1);
    L2 = (res(i,2)+h/2*L1)/Nsq_h2*dN_dz_h2-Nsq_h2./f0^2*Lambda*(res(i,1)+h/2*res(i,2));
    L3 = (res(i,2)+h/2*L2)/Nsq_h2*dN_dz_h2-Nsq_h2./f0^2*Lambda*(res(i,1)+h/2*res(i,2)+h^2/4*L1);
    L4 = (res(i,2)+h*L3)/Nsq(i+1)*dN_dz2-Nsq(i+1)./f0^2*Lambda*(res(i,1)+h*res(i,2)+h^2/2*L2);

    res(i+1,1)=res(i,1) + h*res(i,2)+(L1 + L2 + L3)*h^2/6;
    res(i+1,2)=res(i,2) + (L1 + 2*L2 + 2*L3 + L4)*h/6;
end

end


%%
