clc
clear
close all
global Mav Cav Kav Mnes Cnes Knes typevibration
%% Conctant Parameters
AnalysisType=questdlg('Which Vibration?','Selection','Free Vibraion','Forced Vibration','Free Vibraion')
switch AnalysisType
    case 'Free Vibraion'
     typevibration=0;
    case 'Forced Vibration'
     typevibration=1;  
end
%% aerodynamic properties
M_inf = 4;
U_inf = 300;
rov_inf = 1.288e3;
%% nonlinear energy sink (nes)
Knes=2e10; % 7,8,9,10,11
Mnes=5e-1; %-1,
Cnes=1e8; %5,6,7,8,9
%% plate geometry
a = 0.3;
b = 0.12;
h = 1e-3;
ha = 0.0001;
hs = 0.0001;
%% number of mode shape
nii = 4;
njj = 1;
%
nlayer = 4;
teta =[0,90,0,90];
x1 = 0;
x2 = a;
y1 = 0;
y2 = b;
%
k = 1e4; %
e31 = 1; %
%% mechanical properties
E1 = 138e9;
E2 = 7.2e9;
E = 63e9;
nu = 0.3;
%
nu12 = 0.3;
nu21 = 0.018;
%
rov = 1580;
rovp = 7600; %
%
G12 = 5.5e9;
%
I0 = h*rov;
I1 = (rov*h^3)/12;
%
hbara = k^2*h^2*ha/4 + k*h*ha/2 + ha^3/3;
hbars = k^2*h^2*hs/4 + k*h*hs/2 + hs^3/3;
%% Mode Shape
syms x y
for ii = 1:nii
    for jj = 1:njj  
        W(ii,jj) = sin(ii*pi*x/a)*sin(jj*pi*y/b);   
    end
end
%
for ii = 1:nii
    
    for jj = 1:njj
        
        WT(jj,ii) = sin(ii*pi*x/a)*sin(jj*pi*y/b);
        
    end
end
%
dWx = diff(W,x);
dWTx = diff(WT,x);
dWxx = diff(W,x,2);
dWTxx = diff(WT,x,2);
%
dWy = diff(W,y);
dWTy = diff(WT,y);
dWyy = diff(W,y,2);
dWTyy = diff(WT,y,2);
%
dWxy = diff(diff(W,x),y);
dWTxy = diff(diff(WT,x),y);
%% B & D Matrix
dh = h/nlayer;
z1 = -0.5*h:dh:0.5*h;
nz = numel(z1);
%
Layers=Q_Bar(nlayer,teta,E1,E2,nu12,nu21,G12);
B = zeros(3);
D = zeros(3);
%
for kk = 1:nlayer 
    for ii =1:3
        for jj =1:3
            B(ii,jj) = B(ii,jj) + Layers(kk).Qbar(ii,jj)*((z1(kk+1)^2)-(z1(kk)^2));
            D(ii,jj)=D(ii,jj) + Layers(kk).Qbar(ii,jj)*((z1(kk+1)^3)-(z1(kk)^3));
        end
    end 
end
%
B = (1/2).*B;
D =(1/3).*D;
%
D11 = D(1,1);
D12 = D(1,2);
D22 = D(2,2);
D16 = D(1,3);
D26 = D(2,3);
D66 = D(3,3);
C11 = E/(1-nu^2);
C12 = E*nu/(1-nu^2);
C66 = E/(2*(1+nu));
%% M
C1 = rov*k^3*h^3/12;
C2 = rov*k*h;
C3 = rovp*hbara;
C4 = rovp*ha;
C5 = rovp*hbars;
C6 = rovp*hs;
%
M = C1*int(int(dWx*dWTx,x,0,a),y,0,b) + C1*int(int(dWy*dWTy,x,0,a),y,0,b)+ C2*int(int(W*WT,x,0,a),y,0,b) +...
    C3*int(int(dWx*dWTx,x,x1,x2),y,y1,y2)+ C3*int(int(dWy*dWTy,x,x1,x2),y,y1,y2)+ C4*int(int(W*WT,x,x1,x2),y,y1,y2)+...
    C5*int(int(dWx*dWTx,x,x1,x2),y,y1,y2)+ C5*int(int(dWy*dWTy,x,x1,x2),y,y1,y2)+ C6*int(int(W*WT,x,x1,x2),y,y1,y2);
%
M = double(M);
%% K1
K1 = int(int(D11*dWxx*dWTxx + D12*dWxx*dWTyy + 2*D16*dWxx*dWTxy + D12*dWyy*dWTxx + D22*dWyy*dWTyy + 2*D26*dWyy*dWTxy + 2*D16*dWxy*dWTxx + 2*D26*dWxy*dWTyy + 4*D66*dWxy*dWTxy,x,0,a),y,0,b)+...
    (0.5*C11*hbara)*int(int(dWxx*dWTxx,x,x1,x2),y,y1,y2)+...
    (0.5*C12*hbara)*int(int(dWyy*dWTxx,x,x1,x2),y,y1,y2)+...
    (0.5*C12*hbara)*int(int(dWxx*dWTyy,x,x1,x2),y,y1,y2)+...
    (0.5*C11*hbara)*int(int(dWyy*dWTyy,x,x1,x2),y,y1,y2)+...
    (0.5*C12*hbars)*int(int(dWxx*dWTyy,x,x1,x2),y,y1,y2)+...
    (0.5*C11*hbars)*int(int(dWyy*dWTyy,x,x1,x2),y,y1,y2)+...
    (0.5*C66*4*hbara)*int(int(dWxy*dWTxy,x,x1,x2),y,y1,y2)+...
    (0.5*C66*4*hbars)*int(int(dWxy*dWTxy,x,x1,x2),y,y1,y2);
%
K1 = double(K1);
%% Ka
hhata = k*h/2 + ha;
hhats = k*h/2 + hs;
%
Ka = (0.5*e31*2*hhata)*int(int(dWxx*dWTyy,x,x1,x2),y,y1,y2)+(0.5*e31*2*hhats)*int(int(dWxx*dWTyy,x,x1,x2),y,y1,y2);
Ka = double(Ka);
%% Fdp
Cfdp1 = -rov_inf*U_inf^2/sqrt(M_inf^2 - 1);
Cfdp2 = -(rov_inf*U_inf/sqrt(M_inf^2 - 1))*((M_inf^2 - 2)/(M_inf^2 - 1));
%
Fdp1 = Cfdp1*int(int(dWx*dWTx,x,x1,x2),y,y1,y2);
Fdp1 = double(Fdp1);
%
Fdp2 = Cfdp2*int(int(W*WT,x,x1,x2),y,y1,y2);
Fdp2 = double(Fdp2);
%% Modal Analysis
Wmidd = double(subs(subs(W,x,a/2),b/2));
Mnes = M+Mnes.*Wmidd*Wmidd';  
Cnes = -Fdp2'+Cnes.*Wmidd*Wmidd';
Knes = K1 - Fdp1'+Knes.*Wmidd*Wmidd';
%% structure
Mav = M;
Cav = -Fdp2';
Kav = K1 - Fdp1';
%% ODE Fcn
tend = 10;
Time = 0:0.01:tend;
%% Energy Sink  ode45,ode113, ode23
%
initial_value = [0,0,0,0,0.05,0,0,0];
[t1,y1] = ode45(@OdeFcn,Time,initial_value');
[t2,y2] = ode45(@OdeFcnNES,Time,initial_value');
structure= Wmidd'*y1(:,1:4)';
structureNES= Wmidd'*y2(:,1:4)';
%%
figure(1)
plot(t1,structure)
hold on
plot(t2,structureNES)
xlabel('t(sec)')
ylabel('Wxyt(m)')
xlim([0,1])
legend('no control','energy sink')
axis auto
%%
F=10.*sin(0.5.*Time);
if typevibration==1
    figure(2)
plot(Time,F)
xlabel('t(sec)')
ylabel('disturbance (N)')
axis auto
end



