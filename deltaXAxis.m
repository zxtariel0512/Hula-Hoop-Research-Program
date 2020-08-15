function[deltaEqu] = deltaXAxis(g,Rh,Rg,Rw,c,deltaF,deltaZ,z_mid,f_mid)

syms x y
% express x-axis wrt f and z (y and x)
Rb = sqrt(Rw^2+Rw^2*x^2/c^2);
xAxis = (2*pi*y)^2*(Rh+Rg-Rb)*sind(atand(-Rw^2*x/(Rb*c^2)))/g; 
% partial diff for x and y
diff_xAxis_z = diff(xAxis,x);
diff_xAxis_f = diff(xAxis,y);

xError_func = sqrt((diff_xAxis_z)^2*deltaZ^2+(diff_xAxis_f)^2*deltaF^2);
deltaEqu = subs(xError_func,{x,y},{z_mid,f_mid});