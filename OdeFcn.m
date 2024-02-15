function dx = OdeFcn(t,x)
global Mav Cav Kav  typevibration
dx1 = x(5:8);
dx2 = (-inv(Mav)*Kav)*x(1:4)+(-inv(Mav)*Cav)*x(5:8)+typevibration*10*sin(0.5*t);
dx = [dx1;dx2];
end
