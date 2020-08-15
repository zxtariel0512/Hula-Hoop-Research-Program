function[deltaSta] = deltaYAxis(alpha,g,Rh,Rg,Rw,c,tangent,deltaF,deltaZ,z_mid,f_mid)

syms x y
k = tangent^2*(-(Rw^2+tangent^2*x^2)^(-3/2)*tangent^2*x^2+(Rw^2+tangent^2*x^2)^(-1/2))/(1+(Rw^2+tangent^2*x^2)^(-1/2)*tangent^2*x)^(3/2);
Rb = sqrt(Rw^2+Rw^2*x^2/c^2);
angle = atand(-Rw^2*x/(Rb*c^2));
yAxis = k*(Rh+Rg-Rb)/(sind(angle)*tand(angle));
diff_yAxis_z = diff(yAxis,x);
diff_yAxis_f = diff(yAxis,y);
yError_func = sqrt((diff_yAxis_z)^2*deltaZ^2+(diff_yAxis_f)^2*deltaF^2);
deltaSta = subs(yError_func,{x,y},{z_mid,f_mid});