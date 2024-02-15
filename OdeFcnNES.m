function dx = OdeFcnNES(t,x)
global Mnes Cnes Knes typevibration
dx1 = x(5:8);
dx2 = (-inv(Mnes)*Knes)*x(1:4)+(-inv(Mnes)*Cnes)*x(5:8)+typevibration*10*sin(0.5*t);
dx = [dx1;dx2];
end
