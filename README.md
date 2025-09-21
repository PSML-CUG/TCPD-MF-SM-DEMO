%   Code_Package Name：A three-component polarimetric decomposition method constrained by MF3C DEMO
%  Author:   Wenxin Xue (wenxin@cug.edu.cn), China University of Geosciences (Wuhan)
%  Tutor: Qinghua Xie (xieqh@cug.edu.cn), China University of Geosciences (Wuhan)
%  Date: Sep. 2025


This TCPD-MF-SM-Demo shows the three-component polarimetric decomposition method constrained by MF3C (TCPD-MF) for soil moisture (SM) retrieval.

% Step 1: Utilize SNAP  for data preprocessing to generate the T3 matrix and local incidence angle (IA). Load the resulting .mat file (Data_T3_IA.mat) into MATLAB.

% Step 2: Input the T3 matrix and incidence angle into the function 'Cal_surfaceScattering'. 
% Within this function:
% Compute the scattering power values of the MF3C  decomposition.
% Calculate the proportionality coefficient m to normalize the volume scattering component.

% Step 3: Establish a revised volume scattering coherency matrix based on the methodology described in the paper. Subtract the volume scattering component to isolate the surface scattering contribution.

% Step 4:  Feed the ground scattering coherency matrix into the function 'inversion_Tg_Xbragg', which models ground scattering as the superposition of surface scattering and double-bounce scattering.
% Iterate over the parameter β (beta) to reconstruct candidate ground scattering matrices.
% Compute the Frobenius norm of the residual matrix (observed vs. simulated values).
% Identify the optimal parameter set by minimizing the Frobenius norm, thereby completing the three-component decomposition.

%%% The programs contained in this package are granted free of charge for research and education purposes only.  Scientific results produced using the software provided shall acknowledge the use of this implementation provided by us.  If you plan to use it for non-scientific purposes, don't hesitate to contact us.  Because the programs are licensed free of charge, there is no warranty for the program, to the extent permitted by applicable law.  except when otherwise stated in writing the copyright holders and/or other parties provide the program "as is" without warranty of any kind, either expressed or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose.  the entire risk as to the quality and performance of the program is with you.  should the program prove defective, you assume the cost of all necessary servicing, repair or correction.  In no event unless required by applicable law or agreed to in writing will any copyright holder, or any other party who may modify and/or redistribute the program, be liable to you for damages, including any general, special,  incidental or consequential damages arising out of the use or inability to use the program (including but not limited to loss of data or data being rendered inaccurate or losses sustained by you or third parties or a failure of the program to operate with any other programs), even if such holder or other party has been advised of the possibility of such damages. 

%NOTE: This is just a demo providing a default initialization.  Training is not at all optimized.  Other initializations, optimization techniques, and training strategies may be of course better suited to achieve improved results in this or other problems.  We just did it in the standard way for illustration purposes and dissemination of these models.
