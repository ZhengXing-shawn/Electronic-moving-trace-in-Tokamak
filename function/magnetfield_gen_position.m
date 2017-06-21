function [ B ] = magnetfield_gen_position( y,B_origin,q_safefactor,R_origin,type_mfield)
%generate magnetfield with position.
%  input
%         y   position 
%         B_origin     origin 
%         R_origin     




%% main
% B_orgin=150;
% q=2;
% R_orgin=5;



%% version 1
% r=sqrt((sqrt(y(1)^2+y(2)^2)-R_origin)^2+y(3)^2);
% B_theta=B_origin*r/(q*R_origin);
% B(1)=B_theta*(-y(1)/sqrt(y(1)^2+y(2)^2))*y(3)/r+B_origin*(-y(2))/sqrt(y(1)^2+y(2)^2);
% B(2)=B_theta*(-y(2)/sqrt(y(1)^2+y(2)^2))*y(3)/r+B_origin*(y(1))/sqrt(y(1)^2+y(2)^2);
% B(3)=B_theta*(sqrt((y(1)^2+y(2)^2))-R_origin)/r;

switch type_mfield
    case 'circle' 
        r=sqrt((sqrt(y(1)^2+y(2)^2)-R_origin)^2+y(3)^2);
        q=4-(r/0.8)^2;
        b=R_origin/sqrt(y(1)^2+y(2)^2)*B_origin;
%         B_theta=B_origin/(0.3*q*R_origin);
%         B(1)=B_theta*(-y(1)/sqrt(y(1)^2+y(2)^2))*y(3)+b*(-y(2))/sqrt(y(1)^2+y(2)^2);
%         B(2)=B_theta*(-y(2)/sqrt(y(1)^2+y(2)^2))*y(3)+b*(y(1))/sqrt(y(1)^2+y(2)^2);
%         B(3)=B_theta*(sqrt((y(1)^2+y(2)^2))-R_origin);

        B_theta=r*b/(R_origin*q);
        B(1)=B_theta*(-y(1)/sqrt(y(1)^2+y(2)^2))*(y(3)/r)+b*(-y(2))/sqrt(y(1)^2+y(2)^2);
        B(2)=B_theta*(-y(2)/sqrt(y(1)^2+y(2)^2))*(y(3)/r)+b*(y(1))/sqrt(y(1)^2+y(2)^2);
        B(3)=B_theta*(sqrt((y(1)^2+y(2)^2))-R_origin)/r;
        
    case 'board'
        B(1)=0;
        B(2)=0.01*B_origin*tanh(y(1)/10);
        B(3)=B_origin;
%         B_origin+0.1*(y(2)/R_origin)^2;

    otherwise 
        B=0;


end
end



