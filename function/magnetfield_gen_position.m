function [ B ] = magnetfield_gen_position( y,B_origin,q_safefactor,R_origin,type_mfield)
%generate magnetfield with position.
%  input
%         y            means particle position 
%         B_origin     means B0 (magnetic intensity in the centre of cross section of the toroid)
%         R_origin     means R0 (tokamak toroidal radius)
%   the follow function was base on Pythagorean theorem

switch type_mfield    % determine calculated magnetic field type
    case 'circle' 
        r=sqrt((sqrt(y(1)^2+y(2)^2)-R_origin)^2+y(3)^2);
        q=4-(r/0.8)^2;
        b=R_origin/sqrt(y(1)^2+y(2)^2)*B_origin;


        B_theta=r*b/(R_origin*q);
        B(1)=B_theta*(-y(1)/sqrt(y(1)^2+y(2)^2))*(y(3)/r)+b*(-y(2))/sqrt(y(1)^2+y(2)^2);
        B(2)=B_theta*(-y(2)/sqrt(y(1)^2+y(2)^2))*(y(3)/r)+b*(y(1))/sqrt(y(1)^2+y(2)^2);
        B(3)=B_theta*(sqrt((y(1)^2+y(2)^2))-R_origin)/r;
        
    case 'board'
        B(1)=0;
        B(2)=0.01*B_origin*tanh(y(1)/10);
        B(3)=B_origin;
        
    otherwise 
        B=0;


end
end



