function [ y,model ] = Tokamak_boris( y_position,V_origin,num_position,cal_step_long,R_origin,B_origin,q_safefactor,E,type_mfield )
%UNTITLED3 Summary of this function goes here
%   Accoding to Why is Boris algorithm so good?,Hong Qin,Princeton
    


%% main
   y=zeros(6,num_position);
   h=cal_step_long;

   y(1,1)=y_position(1);    % axix-x0 (the initial fix at x direction)
   y(2,1)=y_position(2);    % axix-y0 (the initial fix at y direction)   
   y(3,1)=y_position(3);    % axix-z0 (the initial fix at z direction)
   y(4,1)=V_origin(1);      % vel-x0  (the initial velocity at z direction)
   y(5,1)=V_origin(2);      % vel-y0  (the initial velocity at y direction)
   y(6,1)=V_origin(3);      % vel-z0  (the initial velocity at z direction)
   B(1)=0;                  % mag-x0  (the initial magnetic intensity at x direction)
   B(2)=0;                  % mag-y0  (the initial magnetic intensity at y direction)
   B(3)=0;                  % mag-z0  (the initial magnetic intensity at z direction)
   E_origin=E;              % Electric field intersity  
   Energy(1)=V_origin(1)^2+V_origin(2)^2+V_origin(3)^2;   % kenetic energy

   c=1;                     % Normalized light constant 
   m=1;                     % Normalized mass 
   q=1;                     % safe-factor
   
   
  for n=1:num_position 
      V_minus=[y(4,n) y(5,n) y(6,n)];       %calculate V_minus as the paper
      B=magnetfield_gen_position([y(1,n),y(2,n),y(3,n)],B_origin,q_safefactor,R_origin,type_mfield);        %calculate the magnetic intensity base on the selected position as a approximates results
      B_scalar=sqrt(B(1)^2+B(2)^2+B(3)^2);      % magnetic intensity
      omega=q*B_scalar/m;       % omega as in the paper
      b=B/B_scalar;             % magnetic field direction
      t=omega*h/2*b;            % t as in the paper
      s=omega*h/(1+(omega*h/2)^2)*b;        % s as in the paper
      V_plus=V_minus+cross((V_minus+cross(V_minus,t)),s);       %calculate V_plus
      
      y(4,n+1)=V_plus(1);  % calculate the particle velocity at the new position
      y(5,n+1)=V_plus(2);
      y(6,n+1)=V_plus(3);

      y(1,n+1)=y(4,n+1)*h+y(1,n);       % calculalte the particle position 
      y(2,n+1)=y(5,n+1)*h+y(2,n);
      y(3,n+1)=y(6,n+1)*h+y(3,n);

      
 
      
%% rounge-kutta 4      
%       Vx=y(4,n);
%       Vy=y(5,n);
%       Vz=y(6,n);
%       Bx=B(1);
%       By=B(2);
%       Bz=B(3);

%       y(4,n+1)=(4*Vx-Bx^2*h^2*Vx+2*Bz*h*Vy+Bx*By*h^2*Vy+2*By*h*Vz+Bx*Bz*h^2*Vz)/(4-Bx^2*h^2-By^2*h^2-Bz^2*h^2-Bx*By*Bz*h^3);
%       y(5,n+1)=(4*Vy-By^2*h^2*Vy+2*Bx*h*Vz+By*Bz*h^2*Vz+2*Bz*h*Vx+By*Bz*h^2*Vx)/(4-Bx^2*h^2-By^2*h^2-Bz^2*h^2-Bx*By*Bz*h^3);
%       y(6,n+1)=((1-1/4*h^2*Bz^2)*(Vz+1/2*h*By*Vx)+(1/2*h*Bx+1/4*h^2*By*Bz)*(Vy+1/2*h*Bz*Vx))/((1-1/4*h^2*By^2)*(1-1/4*h^2*Bz^2)-(1/2*h*Bx+1/4*h^2*By*Bz)^2);
      
%       y(4,n+1)=-(-4*c^2*Vx*m^2-4*Bz*c*h*Vy*m*q+4*By*c*h*Vz*m*q-Bx^2*h^2*Vx*q^2+By^2*h^2*Vx*q^2+Bz^2*h^2*Vx*q^2-2*Bx*By*h^2*Vy*q^2-2*Bx*Bz*h^2*Vz*q^2)/(4*c^2*m^2+Bx^2*h^2*q^2+By^2*h^2*q^2+Bz^2*h^2*q^2);
%       y(5,n+1)=-(-4*c^2*Vy*m^2+4*Bz*c*h*Vx*m*q-4*Bx*c*h*Vz*m*q-2*Bx*By*h^2*Vx*q^2+Bx^2*h^2*Vy*q^2-By^2*h^2*Vy*q^2+2*Bz^2*h^2*Vy*q^2-2*By*Bz*h^2*Vz*q^2)/(4*c^2*m^2+Bx^2*h^2*q^2+By^2*h^2*q^2+Bz^2*h^2*q^2);
%       y(6,n+1)=-(-4*c^2*Vz*m^2-4*By*c*h*Vx*m*q+4*Bx*c*h*Vy*m*q+Bx^2*h^2*Vz*q^2+By^2*h^2*Vz*q^2-Bz^2*h^2*Vz*q^2-2*Bx*Bz*h^2*Vx*q^2-2*By*Bz*h^2*Vy*q^2)/(4*c^2*m^2+Bx^2*h^2*q^2+By^2*h^2*q^2+Bz^2*h^2*q^2);


      
           
%for y(1,n+1) axis-x
   %for y(2,n+1) axis-y
   %for y(3,n+1) axis-z
%    k11=y(4,n);
%    k21=y(5,n);
%    k31=y(6,n);
% 
%    B=magnetfield_gen_position([k11*h/3+y(1,n),k21*h/3+y(2,n),k31*h/3+y(3,n)],B_origin,q,R_origin);
%    k12=y(4,n)+(1/3)*h*(B(3)*k21-B(2)*k31);
%    k22=y(5,n)+(1/3)*h*(B(1)*k31-B(3)*k11);
%    k32=y(6,n)+(1/3)*h*(B(2)*k11-B(1)*k21);
% 
%    B=magnetfield_gen_position([k11*h*(-1/3)+k12*h*1+y(1,n),k21*h*(-1/3)+k22*h*1+y(2,n),k31*h*(-1/3)+k32*h*1+y(3,n)],B_origin,q,R_origin);
%    k13=y(4,n)+(2/3)*h*(B(3)*(k21*(-1/2)+k22*3/2)-B(2)*(k31*(-1/2)+k32*3/2));
%    k23=y(5,n)+(2/3)*h*(B(1)*(k31*(-1/2)+k32*3/2)-B(3)*(k11*(-1/2)+k12*3/2));
%    k33=y(6,n)+(2/3)*h*(B(2)*(k11*(-1/2)+k12*3/2)-B(1)*(k21*(-1/2)+k22*3/2));
% 
%    B=magnetfield_gen_position([k11*h-k12*h+k13*h+y(1,n),k21*h-k22*h+k23*h+y(2,n),k31*h-k32*h+k33*h+y(3,n)],B_origin,q,R_origin);
%    k14=y(4,n)+h*(B(3)*(k21-k22+k23)-B(2)*(k31-k32+k33));
%    k24=y(5,n)+h*(B(1)*(k31-k32+k33)-B(3)*(k11-k12+k13));
%    k34=y(6,n)+h*(B(2)*(k11-k12+k13)-B(1)*(k21-k22+k23));
% 
%    y(1,n+1)=y(1,n)+h/8*(k11+3*k12+3*k13+k14);
%    y(2,n+1)=y(2,n)+h/8*(k21+3*k22+3*k23+k24);
%    y(3,n+1)=y(3,n)+h/8*(k31+3*k32+3*k33+k34);
%    y(4,n+1)=1/8*(k11+3*k12+3*k13+k14);
%    y(5,n+1)=1/8*(k21+3*k22+3*k23+k24);
%    y(6,n+1)=1/8*(k31+3*k32+3*k33+k34);

   Energy(n+1)=y(4,n+1)^2+y(5,n+1)^2+y(6,n+1)^2;        % record the energy information
   
  end

   model.Energy=Energy; 


end
