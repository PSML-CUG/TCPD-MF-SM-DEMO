
%% Function: TCPD-MF three-component polarimetric decomposition method；
% 
% INPUTS
% Elements in the original T3 coherent matrix.
%
% OUTPUTS
% Elements in the surface scattering coherent matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ts11,Ts12,Ts13,Ts22,Ts23,Ts33] = Cal_SurfaceScattering(T11,T12,T13,T22,T23,T33,IncidenceAngle)

T3=[T11 T12 T13;conj(T12) T22 T23;conj(T13) conj(T23) T33];
%%orientation angle estimation by Chen method
phi_lee = 1/4 * atan((-4*real(T3(2,3)))/(-2*T3(2,2)+2*T3(3,3)));
if phi_lee > pi/4
    phi_lee = phi_lee-pi/2;
end
theta1 = real( phi_lee );
R3 = [ 1 0 0;0 cos(2*theta1) sin(2*theta1);0 -sin(2*theta1) cos(2*theta1)];
T3=R3*T3*R3';

% Model-Free 3 component decomposition ( MF3C )
[pd,ps,pv,theta_val,tau_val,dop_fp] = mf3c_T3(T3);

m = real(ps/(ps+pd));


% Volume scattering model
 Tv_FM=[2/4 0 0;0 1/4 0;0 0 1/4];
 Tv = pv*Tv_FM;


 % Ground scattering
 Tp = T3 - Tv;
 pv_FM = Tv(1,1)+Tv(2,2)+Tv(3,3);


    %INVERSION OF GROUND PARAMETERS => 
    %The output 'Tdbl_o' is the estimated coherency matrix for the double-bounce mechanism ; 
    %output 'Tpr_o' is the coherency matrix from estimated parameters from double bounce and direct surface
    
    [error,fs_o,fd_o,beta_o,alpha_o,Tpr_o,Tdbl_o]=inversion_Tg_Xbragg(Tp,real(ps),IncidenceAngle);
    Ts = Tp - Tdbl_o;
    Ts11=Ts(1,1);
    Ts12=Ts(1,2);
    Ts13=Ts(1,3);
    Ts22=Ts(2,2);
    Ts23=Ts(2,3);
    Ts33=Ts(3,3);
    Beta=beta_o;
    Error=error;
 
end



