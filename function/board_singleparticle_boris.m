function [ y,model ] = board_singleparticle_boris( y_position,V_origin,num_position,cal_step_long,R_origin,B_origin,q_safefactor,E_field,type_mfield )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% main
    
   y=zeros(6,num_position);
   h=cal_step_long;
   
   y(1,1)=y_position(1);    % axix-x0
   y(2,1)=y_position(2);    % axix-y0    
   y(3,1)=y_position(3);    % axix-z0
   y(4,1)=V_origin(1);      % vel-x0
   y(5,1)=V_origin(2);      % vel-y0
   y(6,1)=V_origin(3);      % vel-z0
   B(1)=0;                  % mag-x0
   B(2)=0;                  % mag-y0
   B(3)=0;                  % mag-z0
   
   E_origin=E_field;              
   Energy(1)=V_origin(1)^2+V_origin(2)^2+V_origin(3)^2;   %energy
   
   c=1;                     % non-sense
   m=1;                     % non_sense
   q=1;                     % 
   
   
   for n=1:num_position
       
      V_minus=[y(4,n) y(5,n) y(6,n)];
      B=magnetfield_gen_position([y(1,n),y(2,n),y(3,n)],B_origin,q_safefactor,R_origin,type_mfield);
      B_scalar=sqrt(B(1)^2+B(2)^2+B(3)^2);
      omega=q*B_scalar/m;
      b=B/B_scalar;
      t=omega*h/2*b;
      s=omega*h/(1+(omega*h/2)^2)*b;
      V_plus=V_minus+cross((V_minus+cross(V_minus,t)),s);
      
      
      y(4,n+1)=V_plus(1);
      y(5,n+1)=V_plus(2);
      y(6,n+1)=V_plus(3);

      y(1,n+1)=y(4,n+1)*h+y(1,n);
      y(2,n+1)=y(5,n+1)*h+y(2,n);
      y(3,n+1)=y(6,n+1)*h+y(3,n);
      

   
      model.Energy(n)=Energy; 
      
    end  
end
     
   

