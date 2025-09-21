clear;clc;

%% Load data:
load('Data_T3_IA', 'data');

[line,col,bs]=size(data);

% restore T3 matrix
T11=double(data(:,:,1));
T12=double(data(:,:,2)+1i*data(:,:,3));
T13=double(data(:,:,4)+1i*data(:,:,5));
T22=double(data(:,:,6));
T23=double(data(:,:,7)+1i*data(:,:,8));
T33=double(data(:,:,9));
IncidenceAngle = double(data(:,:,10));
IA = IncidenceAngle*pi/180;


%% polarimetric decomposition (TCPD-MF)
TS11=zeros(line,col);
TS12=zeros(line,col);
TS22=zeros(line,col);
TS33=zeros(line,col);


for i=1:line
    for j=1:col
        if mask(i,j)==1
            [Ts11,Ts12,Ts13,Ts22,Ts23,Ts33] = Cal_SurfaceScattering(T11(i,j),T12(i,j),T13(i,j),T22(i,j),T23(i,j),T33(i,j),IA(i,j));
            TS11(i,j)=real(Ts11);
            TS12(i,j)=Ts12;
            TS22(i,j)=real(Ts22);
            TS33(i,j)=real(Ts33);
        end
    end
end

