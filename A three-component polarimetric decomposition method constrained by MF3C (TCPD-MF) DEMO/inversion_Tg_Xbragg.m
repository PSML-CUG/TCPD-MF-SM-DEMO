%%
%% Function: Ground scattering model construction and solution;
%%
% INPUTS
%
% Tp: ground coherency matrix;
% Ps: surface backscattering coefficient;
% theta: incidence angle (radians);
% m: Proportional coefficient.
%
% OUTPUTS
%
% The parameters of the inversion and the coherence matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [error_o,fs_o,fd_o,beta_o,alpha_o,Tpr_o,Tdbl_o]=inversion_Tg_Xbragg(Tp,Ps,theta,m)

t11 = real(Tp(1,1)); t22 = real(Tp(2,2)); t33 = real(Tp(3,3)); t12 = Tp(1,2);
t13 = real(Tp(1,3)); t23 = real(Tp(2,3));
[t11, t12, ~, t22, ~, t33] = Deorientation_Chen(t11, t12, t13, t22, t23, t33);

ts22 = m*t22;
ts33 = m*t33;
sinc4Fai = abs((ts22-ts33)/(ts22+ts33));

%Fai is determined by the Look-Up table method
         if sinc4Fai < 0.0
             Sigma2 = 0.0;
         else
             Sigma2 = log(sinc4Fai)/(-8.0);
         end
         
% Beta values are obtained
[beta,fs_m]=Test_Odd(theta); % fs_m not used

beta = min(beta):0.01:max(beta);
beta_factor = 1./(1 + beta.^2);

error_o=nan;fs_o=nan;fd_o=nan;beta_o=nan;alpha_o=nan;Tpr_o=nan;Tdbl_o=nan;abs_alpha1=nan;abs_alpha2=nan;

MIN_VAL = 1e10;
Tp_aux = Tp;

% They're set to zero for a fair comparison with the reconstructed T matrix
Tp_aux(1,3) = 0; Tp_aux(2,3) = 0; Tp_aux(3,1) = 0; Tp_aux(3,2) = 0;


%% Traverse the values of beta
for ib=1:length(beta)

    % Each parameter is represented by the ground scattering matrix, beta and Ps
    fs = beta_factor(ib)*Ps;
    fd = t22 + t33 - fs*beta(ib)^2;

    cos_term = sqrt((t22-(fs*beta(ib)^2*(1+exp(-8*Sigma2)))/2)/fd);
    phase_dbl = angle( (t12 - fs*beta(ib)*exp(-2*Sigma2))/(fd*cos_term) );
    abs_alpha = abs( (t12 - fs*beta(ib)*exp(-2*Sigma2))/(fd*cos_term) );
    alphar = abs_alphar*exp(1j*phase_dbl);


    % Reconstruction the Tg
    t11r = fd*abs_alphar^2 + fs;
    t22r = fd*cos_term^2 + (fs*beta(ib)^2*(1+(exp(-8*Sigma2))))/2;
    t33r = fd*(1-cos_term^2)+ (fs*beta(ib)^2*(1-(exp(-8*Sigma2))))/2;
    t12r = fd*alphar*cos_term + fs*beta(ib)*exp(-2*Sigma2);

    Tpr = [t11r t12r 0;
           t12r' t22r 0;
           0   0   t33r];
    
    % Calculate the Frobenius norm of the matrix
    error = norm(Tp_aux-Tpr,"fro");  

    if error<MIN_VAL

        fs_o = fs; fd_o = fd;
        beta_o = beta(ib);
        alpha_o = alphar;

        Tpr_o = Tpr;

        Tdbl_o = [fd*abs_alphar^2            fd*alphar*cos_term  0;
                  conj(fd*alphar*cos_term)   fd*cos_term^2      0;
                            0                     0            fd*(1-cos_term^2)];

        error_o = error;

        MIN_VAL = error;

    end

end

%%Test_Odd is an function for generating all possible values for beta with specific incidence angle

function [beta,fs_m]=Test_Odd(theta)

%The default scope for epsoil is [2,41] by Hajnsek 2009
epsoil = linspace(2,41,100);
Rh = ( cos(theta) - sqrt(epsoil-sin(theta)^2) ) ./ ...
         ( cos(theta) + sqrt(epsoil-sin(theta)^2) );
Rv = ( (epsoil-1).*(sin(theta)^2-epsoil*(1+sin(theta)^2)) ) ./ ...
         ( epsoil*cos(theta) + sqrt(epsoil-sin(theta)^2) ).^2;

%For T matrix
beta = (Rh-Rv)./(Rh+Rv) ;
fs_m = 1/2*abs(Rh+Rv).^2;

end
