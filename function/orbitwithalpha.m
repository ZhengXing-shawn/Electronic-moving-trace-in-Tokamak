% calculate the centre of electrinic moving trace in tokamak at non-othogonal 
% reference to White R B. The theory of toroidally confined plasmas[M].World Scientific Publishing Co Inc, 2013.
% all the follows functions was base at white books,in which more detail can be found. 


function [dy] = orbitwithalpha( t,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        
    global q q1 q2 q3 mu psiw;
    
    dy=zeros(4,1);
    zeta=y(1); theta=y(2); psip=y(3); rhopara=y(4);
    psin=psip/psiw; % normalization psip
    q=q1+q2*psin+q3*psin^2;
    
    % field and other parameters
    psi=psip*(q1+q2/2*psin+q3/3*psin^2);
    r=sqrt(2*psi); % r/R0
    R=1+r*cos(theta);
    g=1.0; % g: poloidal current outside psi
    gppsi=0; % p(prime): derivative of psi
    I=r^2/q; % I: toroidal current inside psi
    Ippsi=2/q^3*(q^2-psin*(q2+2*q3*psin)*(q1+q2/2*psin+q3/3*psin^2));
    Bt=g/R;
    Bp=r/(q*R); % Bp=r/(q*R)
    B=sqrt(Bt^2+Bp^2);
    Bppsip=(r*(1-r^2*R*(q2+2*q3*psin)/psiw/q^2)-g^2*q^2*cos(theta))/(B*R^3*r*q); % dB/dpsip
    Bptheta=B*r*sin(theta)/(R); % dB/dtheta
    % Bpzeta=0;
    zoom=0.0005;
    
%    alpha21
    alpha=zoom*3*exp(-(psip-0.02)^2/0.008374^2)*sin(2*theta-zeta);
    alphappsip=zoom*85562.9*exp(-(psip-0.02)^2/0.008374^2)*sin(-2*theta+zeta)*(psip-0.02);
    alphaptheta=zoom*6*exp(-(psip-0.02)^2/0.008374^2)*cos(-2*theta+zeta);
    alphapzeta=-zoom*2*exp(-(psip-0.02)^2/0.008374^2)*cos(-2*theta+zeta);
    %alpha11
%     alpha=zoom*3*exp(-(psip-0.015)^2/0.008374^2)*sin(theta-zeta);
%     alphappsip=zoom*85562.9*exp(-(psip-0.015)^2/0.008374^2)*sin(-theta+zeta)*(psip-0.015);
%     alphaptheta=zoom*6*exp(-(psip-0.015)^2/0.008374^2)*cos(theta+zeta);
%     alphapzeta=-zoom*2*exp(-(psip-0.015)^2/0.008374^2)*cos(theta+zeta);
    
    
    D=g*q+I+rhopara*(g*Ippsi-I*gppsi);

%     E=rhopara^2*B^2/2+mu*B;

    dy(1)=rhopara*B^2*(q+(rhopara+alpha)*Ippsi+I*alphappsip)/D-(mu+rhopara^2*B)*I*Bppsip/D;
    dy(2)=rhopara*B^2*(1-(rhopara+alpha)*gppsi-g*alphappsip)/D+(mu+rhopara^2*B)*g*Bppsip/D;
    dy(3)=-g*(mu+rhopara^2*B)*Bptheta/D+rhopara*B^2*g*alphaptheta/D-rhopara*B^2*I*alphapzeta/D;
    dy(4)=-(mu+rhopara^2*B)*((1-(rhopara+alpha)*gppsi-g*alphappsip)*Bptheta)/D+(I*alphapzeta-g*alphaptheta)/D*((mu+rhopara^2*B)*Bppsip);
    
end

