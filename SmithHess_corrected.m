function mainPanelSingle

% written by Ebenezer P Gnanamanickam while teaching AAE 490D Aerodynamic
% design at Purdue University in the Spring of 2012. The code is entirely
% self containing and it calculates the Lift coefficient and Pressure
% distribution around a NACA 4 digit airfoil.
% It follows almost identically the Smith Hess Method as outlined in
% Moran's Aerodynamics text (the exact name skips my memory).
% You can use it for any arbitrary airfoil by changing the subroutine NACA4
% to input a bunch of x,z points. The points start from the trailing edge
% go through the bottom surface of the airfoil and then over the top and
% finish at the trailing edge again. It should be non-dimensionalized with
% respect to the chord. So it should start with the point (1,0) and end
% with the point (1,0)
% The inputs are right below. They are the 
% 1) airfoil number (af_no) if using NACA4 digit series
% 2) number of panels (no_pls) 
% 3) angle of attack (alpha)
clc;
clear all;
close all;

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
af_no=1408;
no_pls=100;
alpha=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,z]=NACA4(af_no,no_pls,0);
% If you want to input points from a dat file, please comment out the line
% above and use the following format to import points
% a=load('n0012.dat');
% x=a(:,1);
% z=a(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sth,cth]=setUpPanel1(x,z);
V_inf=1;

figure;
plot(x,z,'o-r')
axis equal
grid;
xlabel('x/c')
ylabel('z/c')

[A,b,beta,r]=calCoeffs1(x,z,...
    sth,cth,alpha,V_inf);
% 

xx=inv(A)*b';
q=xx(1:end-1);
gamma=xx(end)
[Vt]=calVt(sth,cth,r,beta,V_inf,alpha,q,gamma);

figure;
plot(x(1:round(length(x)/2)),1-(Vt(1:round(length(x)/2))./V_inf).^2,'r',...
    x(round(length(x)/2)+1:end),1-(Vt(round(length(x)/2)+1:end)./V_inf).^2,'b');
grid;
legend('Bottom surface','Top Surface','Location','SouthWest')
set(gca,'YDir','reverse');
xlabel('x/c')
ylabel('C_p')

C_l=4*(gamma)/V_inf

return
% % % theta*180/pi
% % % C_l/(alpha*180/pi)

function [x,z]=NACA4(af_no,no_pls,doPlot)

x_c=linspace(0,pi,round(no_pls/2));
x_c=(1-cos(x_c))./2;

t=(af_no-floor(af_no/100)*100)/100; % maximum thickness in tenth's of chord
m=floor(af_no/1000)/100; %maximum thickness in tenth's of chord
p=(floor(2412/100)-floor(2412/1000)*10)/10; %position of maximum thickness in tenth's of chord
zu=(t/0.2).*(.2969.*sqrt(x_c)-.126.*x_c-.35160.*x_c.^2+...
    .2843.*x_c.^3-.1036.*x_c.^4);

for i=1:1:length(x_c)
    if x_c(i)<=p
        z_c(i)=m*(1/p)^2.*(2.*p.*x_c(i)-x_c(i).^2);
    else
        z_c(i)=m*(1/(p-1))^2.*((1-2*p)+2.*p.*x_c(i)-x_c(i).^2);
    end
end
z_u=zu+z_c;
z_l=z_c-zu;

% P=polyfit(z_u(end-round(no_pls*.1):end),x_c(end-round(no_pls*.1):end),2);
% z_u(end+1)=0;
% x_c(end+1)=polyval(P,z_u(end));
% z_l(end+1)=0;

if doPlot==1
    h6=figure;
    plot(x_c,z_u,'o-r',x_c,z_l,'o-b',x_c,z_c,'-k');
    grid;
    axis equal;
end

z=[fliplr(z_l(2:end)) z_u(1:end-1)]';
x=[fliplr(x_c(2:end)) x_c(1:end-1)]';
z(1)=0;

return;

function [A,b,beta,r]=calCoeffs1(x,z,...
    sth,cth,alpha,V_inf)

A=zeros(length(x)+1,length(x)+1);
b=zeros(1,length(x)+1);

[r,beta]=calRBeta(x,z);

alpha=alpha/180*pi;

for i=1:1:length(x)
    sum=0;
    for j=1:1:length(x)
        if j==length(x)
            mm=1;
        else
            mm=j+1;
        end
        
        sd=sinDiff(sth,cth,i,j);
        cd=cosDiff(sth,cth,i,j);
        
        A(i,j)=sd*log(r(i,mm)/r(i,j))+cd*beta(i,j);
        sum=sum+cd*log(r(i,mm)/r(i,j))-sd*beta(i,j);
        
      
        A(i,j)=A(i,j)/(2*pi);
    end
    A(i,length(x)+1)=sum/(2*pi);
    
    s1=(cos(alpha)*sth(i)-sin(alpha)*cth(i));
    b(i)=V_inf*s1;
end

tt=[1 length(x)];

for j=1:1:length(x)
    if j==length(x)
        mm=1;
    else
        mm=j+1;
    end
    sum=0;
    for k=1:1:length(tt)
        sd=sinDiff(sth,cth,tt(k),j);
        cd=cosDiff(sth,cth,tt(k),j);
        
       sum=sum+(sd*beta(tt(k),j)-...
            cd*log(r(tt(k),mm)/r(tt(k),j)));
    end
    A(length(x)+1,j)=sum/(2*pi);
end


sum2=0;
for k=1:1:length(tt)
    sum1=0;
    for j=1:1:length(x)
        if j==length(x)
            mm=1;
        else
            mm=j+1;
        end
        sd=sinDiff(sth,cth,tt(k),j);
        cd=cosDiff(sth,cth,tt(k),j);
        
        sum1=sum1+(sd*log(r(tt(k),mm)/r(tt(k),j))+...
            cd*beta(tt(k),j));
    end
    sum2=sum2+sum1;
end

A(length(x)+1,length(x)+1)=sum2/(2*pi);


c1=(sin(alpha)*sth(1)+cos(alpha)*cth(1));
c2=(sin(alpha)*sth(length(x))+cos(alpha)*cth(length(x)));
b(length(x)+1)=-V_inf*c1-V_inf*c2;
return;

function [sth,cth,len]=setUpPanel1(x,z)

N=length(x);

len=zeros(1,N);
sth=len;
cth=len;

for i=1:1:N
   if i==length(x)
       mm=1;
   else
       mm=i+1;
   end
   len(i)=sqrt((x(mm)-x(i))^2+(z(mm)-z(i))^2);
   sth(i)=(z(mm)-z(i))/len(i);
   cth(i)=(x(mm)-x(i))/len(i);
   
end

return;

function [Vt]=calVt(sth,cth,r,beta,V_inf,alpha,q,gamma)

alpha=alpha/180*pi;


N=size(beta,1);
Vt=zeros(1,N);
for i=1:1:N
    sum1=0;
    sum2=0;
    for j=1:1:N
        if j==N
            mm=1;
        else
            mm=j+1;
        end
        
        sd=sinDiff(sth,cth,i,j);
        cd=cosDiff(sth,cth,i,j);
               
        sum1=sum1+q(j)/(2*pi)*...
            (sd*beta(i,j)-cd*log(r(i,mm)/r(i,j)));
        sum2=sum2+gamma/(2*pi)*...
            (sd*log(r(i,mm)/r(i,j))+cd*beta(i,j));
    end
    c1=sin(alpha)*sth(i)+cos(alpha)*cth(i);
    Vt(i)=V_inf*c1+sum1+sum2;
end

return;

function [r,beta]=calRBeta(x,z)

r=zeros(length(x),length(x));
beta=r;
for i=1:1:length(x)
    if i==length(x)
        xcntr=(x(1)+x(i))/2;
        zcntr=(z(1)+z(i))/2;
    else
        xcntr=(x(i+1)+x(i))/2;
        zcntr=(z(i+1)+z(i))/2;
    end
    
    for j=1:1:length(x)
        r(i,j)=sqrt((xcntr-x(j))^2+(zcntr-z(j))^2);
        
        if i~=j
            if j==length(x)
                k=1;
            else
                k=j+1;
            end
            temp1=(zcntr-z(k))*(xcntr-x(j))-...
                (xcntr-x(k))*(zcntr-z(j));
            temp2=(xcntr-x(k))*(xcntr-x(j))+...
                (zcntr-z(k))*(zcntr-z(j));
            beta(i,j)=atan2(temp1,temp2);
        else
            beta(i,j)=pi;
        end
        
    end
end

function [sinD]=sinDiff(sin_theta,cos_theta,i,j)

sinD=(cos_theta(j)*sin_theta(i)-sin_theta(j)*cos_theta(i));

return;

function [cosD]=cosDiff(sin_theta,cos_theta,i,j)

cosD=(sin_theta(j)*sin_theta(i)+cos_theta(j)*cos_theta(i));

return;