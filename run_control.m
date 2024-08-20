clc;clear;close all

day2sec=24*60*60;

h=500;               % total depth
dt=1800;             % time step
nb=10;               % number of biological tracers
nz=80;               % number of grids
nt=2*365*day2sec/dt; % number of steps for simulations
nr=day2sec/dt;       % number of steps for record freqeuncy

% set-up grid
zr=linspace(-h,0,nz)';
zu=convn(zr,[.5;.5],'valid');
dz=mean(diff(zr));
C=zeros(nz,2,nb);

%% set-up initial conditions
C(:,1,1)=ones(size(zr))*7;    % Silicate
C(:,1,2)=ones(size(zr))*0.1;  % Phytoplankton group A
C(:,1,3)=ones(size(zr))*0.1;  % Zooplankton A
C(:,1,4)=ones(size(zr))*0.1;  % Detritus A

C(:,1,5)=ones(size(zr))*13;   % Nitrate
C(:,1,6)=ones(size(zr))*0.1;  % Phytoplankton group B
C(:,1,7)=ones(size(zr))*0.1;  % Zooplankton B,C
C(:,1,8)=ones(size(zr))*0.1;  % Detritus B,C

C(:,1,9)=ones(size(zr))*0.1;  % Phytoplankton group C
C(:,1,10)=ones(size(zr))*0.9; % Phosphate

B=zeros(nz,floor(nt/nr),nb);
I=zeros(nz,floor(nt/nr));
t=zeros(1,floor(nt/nr));

%% model coefficients
parfrac=0.43;
kappa_w=0.067;   % water light atten. coef.
kapaa_i=1.5;     % ice light atten. coef.
RS=2.0;          % silica:nitrogen ratio 
RP=1/16;         % phosphorous:nitrogen ratio
sedrr=0.75;

% lateral subsurface nutrient supply
UL=4e-7;                             % Lateral nutrient supply rate
no3c=exp(-((zr+175)/70).^2)*1+13;    % nitrate background conc.
sio4c=exp(-((zr+175)/70).^2)*27+7;   % silicate background conc.
po4c=exp(-((zr+175)/70).^2)*0.9+0.9; % phosphate background conc.

% phytoplankton group A
Vm_a=1.5;        % maximum growth rate
kappa_pa=0.0095; % phytoplankton light atten. coef.
Ek_a=2;          % light saturation coef.
kn_a=6.0;        % nitrate half sat. coef.
ks_a=13.0;       % silicate half sat. coef.
kp_a=1.0/16;     % phosphate half sat. coef.
sigma_a=0.1;     % mortality rate

% phytoplankton group B
Vm_b=1.4;        % maximum growth rate
kappa_pb=0.0095; % phytoplankton light atten. coef.
Ek_b=6;          % light saturation coef.
kn_b=0.1;        % nitrate half sat. coef.
kp_b=0.5/16;     % phosphate half sat. coef.
sigma_b=0.1;     % mortality rate

% phytoplankton group C
Vm_c=0.74;       % maximum growth rate
kappa_pc=0.0095; % phytoplankton light atten. coef.
Ek_c=12;         % light saturation coef.
kn_c=0.0;        % nitrate half sat. coef.
kp_c=0.4/16;     % phosphate half sat. coef.
sigma_c=0.1;     % mortality rate

% zooplankton group A
Rm_a=0.52;       % maximum growth rate
lambda_a=0.84;   % phytoplankton sat. coef.
gamma_a=0.1;     % transfer effeciency
xi_a=2.175;      % quadratic mortality rate

% zooplankton group BC
Rm_bc=0.40;      % maximum growth rate
lambda_bc=0.84;  % phytoplankton sat. coef.
gamma_bc=0.1;    % transfer effeciency
xi_bc=2.175;     % quadratic mortality rate

% detritus A
delta_a=0.343;    % remineralization rate
wd_a=-24;         % sinking speed

% detritus BC
delta_bc=1.03;  % remineralization rate
wd_bc=-8;       % sinking speed

%% forcing configuration
load('forcing.mat')
t_frc=[t_frc(1)-mean(diff(t_frc));t_frc; t_frc(end)+mean(diff(t_frc))];
I0=[I0(end);I0;I0(1)];
ice_f=[ice_f(end);ice_f;ice_f(1)];
ice_h=[ice_h(end);ice_h;ice_h(1)];

%% solving governing equations
tic
for j=1:nt
    disp(j/(365*day2sec/dt)*100)
    dy=mod(dt*j/day2sec,365); % day of year
    
    % biogeochem. processes
    SiOH4=C(:,1,1);
    Pa   =C(:,1,2);
    Za   =C(:,1,3);
    Da   =C(:,1,4);
    NO3  =C(:,1,5);
    Pb   =C(:,1,6);
    Zbc  =C(:,1,7);
    Dbc  =C(:,1,8);
    Pc   =C(:,1,9);
    PO4  =C(:,1,10);
    
    % light attenuation
    Ii=zeros(nz,1);
    Ii(nz)=interp1(t_frc,I0,dy);
    icefi=interp1(t_frc,ice_f,dy);
    icehi=interp1(t_frc,ice_h,dy);

    % light attenuation by sea ice
    Ii(nz)=parfrac*Ii(nz)*exp(-kapaa_i*icefi*icehi);
    % light attenuation by water & phytoplanktons
    for i=nz-1:-1:1
        Ii(i)=Ii(i+1)/(1+dz*(kappa_w+...
            kappa_pa*Pa(i)+kappa_pb*Pb(i)+kappa_pc*Pc(i)));
    end
    
    % zooplankton growth rate
    G_a=Rm_a*(1-exp(-lambda_a*Pa));
    G_bc=Rm_bc*(1-exp(-lambda_bc*(Pb+Pc)));

    % phytoplankton growth rate
    U_a=Vm_a*(1-exp(-Ii/Ek_a))...
        .*min([NO3./(NO3+kn_a), SiOH4./(SiOH4+ks_a), PO4./(PO4+kp_a)],[],2);
    U_b=Vm_b*(1-exp(-Ii/Ek_b))...
        .*min([NO3./(NO3+kn_b), PO4./(PO4+kp_b)],[],2);
    U_c=Vm_c*(1-exp(-Ii/Ek_c))...
        .*PO4./(PO4+kp_c);

    % predation probability
    Pib=max(Pb./(Pb+Pc),0);
    Pic=max(Pc./(Pb+Pc),0);
    
    %% biogeochemical processes
    % silicate
    C(:,1,1)=SiOH4+dt/day2sec*(delta_a*Da+gamma_a*G_a.*Za-U_a.*Pa)*RS;
    % phytoplankton A
    C(:,1,2)=Pa+dt/day2sec*(U_a.*Pa-G_a.*Za-sigma_a*Pa);
    % zooplankton A
    C(:,1,3)=Za+dt/day2sec*((1-gamma_a)*G_a.*Za-xi_a*Za.^2);
    % detritus A
    C(:,1,4)=Da+dt/day2sec*(sigma_a*Pa+xi_a*Za.^2-delta_a*Da);
    % nitrate
    C(:,1,5)=NO3+dt/day2sec*(delta_a*Da+gamma_a*G_a.*Za-U_a.*Pa...
                            +delta_bc*Dbc+gamma_bc*G_bc.*Zbc-U_b.*Pb);
    % phytoplankton B
    C(:,1,6)=Pb+dt/day2sec*(U_b.*Pb-Pib.*G_bc.*Zbc-sigma_b*Pb);
    % zooplankton BC
    C(:,1,7)=Zbc+dt/day2sec*((1-gamma_bc)*G_bc.*Zbc-xi_bc*Zbc.^2);
    % detritus BC
    C(:,1,8)=Dbc+dt/day2sec*(sigma_b*Pb+sigma_c*Pc+xi_bc*Zbc.^2-delta_bc*Dbc);
    % phytoplankton C
    C(:,1,9)=Pc+dt/day2sec*(U_c.*Pc-Pic.*G_bc.*Zbc-sigma_c*Pc);
    % phosphate
    C(:,1,10)=PO4+dt/day2sec*(delta_a*Da+gamma_a*G_a.*Za-U_a.*Pa...
                            +delta_bc*Dbc+gamma_bc*G_bc.*Zbc-U_b.*Pb-U_c.*Pc)*RP;

    %% physical processes
    % nutrient supply via lateral transport
    ul=zr*0+UL; ul(zr>-150)=0;
    C(:,2,1)=C(:,1,1)+dt*(-ul.*(C(:,1,1)-sio4c));
    C(:,2,5)=C(:,1,5)+dt*(-ul.*(C(:,1,5)-no3c));
    C(:,2,10)=C(:,1,10)+dt*(-ul.*(C(:,1,10)-po4c));
    C(:,1,[1 5 10])=C(:,2,[1 5 10]);

    % mixing coefficient model
    eps_mix=5e-8;
    drhodz1=-0.0176;
    drhodz2=-0.0048;
    Nf1=sqrt(-10/1025*drhodz1);
    Nf2=sqrt(-10/1025*drhodz2);
    ice_h0=interp1(t_frc,ice_h,dy);
    ice_h1=interp1(t_frc,ice_h,dy+dt/day2sec);
    if ice_h1-ice_h0>0 && ice_h0>1e-4
        Ai=ones(size(zu))*1e-2;
    elseif ice_h1-ice_h0<0 && ice_h0>1e-4
        Ai=ones(size(zu))*1e-5;
    elseif ice_h0<1e-4
        Ai=ones(size(zu))*0.25*eps_mix./(Nf1.^2);
    end
        Ai(zu<=-250)=0.25*eps_mix./(Nf2.^2);

    % vertical mixing
    Af=[0;Ai;0];
    alp=dt/dz^2;
    F=spdiags([-alp*Af(2:end) alp*conv(Af,[1 1],'valid')+1 -alp*Af(1:end-1)],-1:1,nz,nz);
    for i=1:size(C,3)
        C(:,2,i)=F\C(:,1,i);
    end        
    C(:,1,:)=C(:,2,:);
    
    for i=1:nz
        % advection (detritus sinking)
        if i==1      % bottom B.C.
        % sediment (remineralization & denitrification) process
            C(i,2,4)=C(i,1,4)-dt/day2sec/dz*wd_a*(C(i+1,1,4)-C(i,1,4));
            C(i,2,8)=C(i,1,8)-dt/day2sec/dz*wd_bc*(C(i+1,1,8)-C(i,1,8));
            C(i,2,1)=C(i,1,1)-dt/day2sec/dz*wd_a*C(i,1,8)*RS;
            C(i,2,5)=C(i,1,5)-dt/day2sec/dz*(wd_a*C(i,1,4)+wd_bc*C(i,1,8))*sedrr;
            C(i,2,10)=C(i,1,10)-dt/day2sec/dz*(wd_a*C(i,1,4)+wd_bc*C(i,1,8))*RP;
        elseif i==nz % surface B.C.
            C(i,2,4)=C(i,1,4)+dt/day2sec/dz*wd_a*C(i,1,4);
            C(i,2,8)=C(i,1,8)+dt/day2sec/dz*wd_bc*C(i,1,8);
        else         % interior
            C(i,2,4)=C(i,1,4)-dt/day2sec/dz*wd_a*(C(i+1,1,4)-C(i,1,4));
            C(i,2,8)=C(i,1,8)-dt/day2sec/dz*wd_bc*(C(i+1,1,8)-C(i,1,8));
        end
    end
    C(:,1,:)=C(:,2,:);
    
    % blow-up test
    if any(isnan(C(:)))
        error('model blow-up')
    end
    
    
    if mod(j,nr)==0
        B(:,j/nr,:)=C(:,1,:);
        I(:,j/nr)=Ii;
        t(j/nr)=dt*j/day2sec;
    end
end
toc

save('output_control')