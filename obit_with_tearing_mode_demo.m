% calculate the centre of electrinic moving trace in tokamak at non-othogonal 
% reference to White R B. The theory of toroidally confined plasmas[M].World Scientific Publishing Co Inc, 2013.
% all the follows functions was base at white books,in which more detail can be found. 


close all; clear; clc;
   
    global q1 q2 q3 mu psiw;
    psiw=0.043636; % psiw=0.043636, a=0.40
    q1=1.0; q2=1.0; q3=1.0;
    
    a=sqrt(2*psiw*(q1+q2/2+q3/3)); % minor radius
    
    % initial values
    psip0=0.6*psiw; theta0=0*pi/4; zeta0=0;
    psin0=psip0/psiw;
    r0=sqrt(2*psip0*(q1+q2/2*psin0+q3/3*psin0^2));
    q=q1+q2*psin0+q3*psin0^2;

    % calculate initial parameters
    R=1+r0*cos(theta0);
    g=1.0; Bt=g/R;
    Bp=r0/(q*R);
    B=sqrt(Bt^2+Bp^2);
    lambda0=0.9; % pitch angle = mu*B/E
    % v ~ rhoi/Omega ~ (10^-3)R0/Omega, then E0 ~ v^2 ~ 10^-6
    E0=1/5e4; % sensitive, E0> e.g., 1/2e6..., banana
    mu=lambda0*E0;
    drc=-1; % v0= +/- give different orbits
    v0=sqrt(2*E0);
	rhopara0=drc*v0*sqrt(1-lambda0*B)/B;
    
    % solve
    y0=[zeta0, theta0, psip0, rhopara0];
    options=odeset('RelTol',1e-10,'AbsTol',[1e-10 1e-10 1e-11 1e-10],'MaxStep',0.1);
%     [t,y] = ode45(@orbit,0:200,y0,options);
    tend=200/abs(rhopara0); dt=tend/2e4;
    [t,y] = ode45(@orbitwithalpha2,0:dt:tend,y0);
    
    % plot
    figure; set(gcf,'DefaultAxesFontSize',15);
    zeta=y(:,1); theta=y(:,2); psip=y(:,3); rhopara=y(:,4);
    psi=psip.*(q1+q2/2*psip./psiw+q3/3*(psip./psiw).^2);
    r=sqrt(2*psi); R=1+r.*cos(theta); x2=r.*cos(theta); y2=r.*sin(theta);
    x3=R.*cos(zeta); y3=R.*sin(zeta); z3=r.*sin(theta);
    subplot(221);plot(t,zeta,'r-','LineWidth',2);
    axis tight; xlabel('t/\Omega_c'); ylabel('\zeta'); 
    title(['\zeta-t',', E=',num2str(E0),', \Lambda=',num2str(lambda0)]);
    subplot(222);plot(t,rhopara,'r-','LineWidth',2); axis tight; 
    xlabel('t/\Omega_c'); ylabel('\rho_{||}'); 
    title(['\rho_{||}-t',', r0=',num2str(r0),',q=',num2str(q)]);
    subplot(223); plot(x2,y2,'LineWidth',2); axis equal;
    xlabel('x'); ylabel('y'); 
    title(['poloidal projection, a=',num2str(a)]); hold on;
    plot(a.*cos(0:pi/20:2*pi),a.*sin(0:pi/20:2*pi),'r--');
    subplot(224); plot3(x3,y3,z3,'g',x3(1),y3(1),z3(1),'r*');    
    title(['direction=',num2str(drc)]);
    axis equal; % title('3D plot'); view(2);
    
    print(gcf,'-dpng',['E=',num2str(E0),',r=',num2str(r0),',Lambda=',...
    num2str(lambda0),',a=',num2str(a),...
    ',drc=',num2str(drc),'.png']);