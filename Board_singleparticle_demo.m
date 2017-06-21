%% describe
% output:
%         y(1)        position_x
%         y(2)        position_y
%         y(3)        volecity_x
%         y(4)        volecity_y

clear

%% mian

y_position=[-0.5,0,0];       % XYZ position
V_origin=[0.03,0.04,0];   % V_x,V_y,V_z 
num_position=1e6;         % calculate time
cal_step_long=1e-3;       % step size
R_origin=5;               % R0
B_origin=2;             % B0
q_safefactor=2;          %non sense
E_field=0;                      %Electric field
type_mfield='board';

[y,model]=board_singleparticle_boris(y_position,V_origin,num_position,cal_step_long,R_origin,B_origin,q_safefactor,E_field,type_mfield);