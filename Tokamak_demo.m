%% describe
% output:
%         y(1)        position_x
%         y(2)        position_y
%         y(3)        position_z
%         y(4)        volecity_x
%         y(5)        volecity_y
%         y(6)        volecity_z
%    
     
clear

%% mian

y_position=[0,5.15,0.6];       % XYZ position
V_origin=[0.001,0.017,0.1];   % V_x,V_y,V_z 
num_position=2e7;         % calculate time
cal_step_long=1e-4;       % step size
R_origin=5;               % R0
B_origin=8;             % B0
q_safefactor=2;                      % safe factor
E=0;                      %Electric field
type_mfield='circle';



% [ y,model ] = Tokamak( y_position,V_origin,num_position,cal_step_long,R_origin,B_origin,q);
[ y,model ] = Tokamak_boris( y_position,V_origin,num_position,cal_step_long,R_origin,B_origin,q_safefactor,E,type_mfield);