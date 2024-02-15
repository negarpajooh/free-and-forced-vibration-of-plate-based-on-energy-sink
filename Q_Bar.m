function Layers=Q_Bar(nlayer,Teta,E1,E2,nu12,nu21,G12)
Q = Qmatrix(E1,E2,nu12,nu21,G12);
T = zeros(3);
for i=1:nlayer
    teta = Teta(i);
    tetaR=(pi*teta)/180;
    T(1,1)=cos(tetaR)*cos(tetaR);
    T(1,2)=sin(tetaR)*sin(tetaR);
    T(1,3)=-2*sin(tetaR)*cos(tetaR);
    T(2,1)=T(1,2);
    T(2,2)=T(1,1);
    T(2,3)=-T(1,3);
    T(3,1)=-0.5*T(1,3);
    T(3,2)=-T(3,1);
    T(3,3)=T(1,1)-T(1,2);
    Layers(i).Qbar=(inv(T))*Q*((inv(T))');
end
end