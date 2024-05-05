%This function is to graph confidence regions and confidence intervals
%for z0, based on a SINGLE SAA solution.



function graph2d_oneRep(zN,Q,r,simCI,tilde_zN,indCI,z0)
%Inputs: 
%zN: the SAA solution, also the center of the confidence region

%Q and r: the matrix and the constant that defines the confidence region
%{z|(z-zN)^T Q (z-zN) <= r}

%simCI: the half lengths of the simultaneous confidence intervals around
%zN. The simultaneous confidence intervals are given by 
%[zN(i)-simCI(i),zN(i)+simCI(i)]_i
%These intervals form the minimal bounding box of the confidence region

%tilde_zN: the center of the individual confidence intervals

%indCI: the half lengths of the individual confidence intervals for 
%z0. The individual confidence intervals are given by 
%[tilde_zN(i)-indCI(i),tilde_zN(i)+indCI(i)]_i

%optional input: z0: the true soluiton, used to check how the confidence intervals work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graph confidence region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%new codes that graph the ellipse {z|(z-zN)^T Q (z-zN) <= r} by using cos and sin


%bdryZ: 1000 points that satisfy z-zN)^T Q (z-zN) = r

bdryZ=zeros(2,1000);

tmp=zeros(2,1);

for i=0:1:999
    %angle: the angle that changes from almost 0 to almost 2pi.
    angle=i/1000*2*pi;
    tmp(2)= sqrt( r * Q(1,1)/ (Q(2,2)* Q(1,1)-Q(1,2)^2 )) * sin(angle);
    tmp(1)= - Q(1,2)/Q(1,1)*tmp(2)+ sqrt(r/Q(1,1))* cos(angle);
    bdryZ(:,i+1)=zN+tmp;
end
hold off;
plot(bdryZ(1,:),bdryZ(2,:));



% %old codes that graph the confidence region by contours of the quadratic
% %function
% %I use the contour of the quadratic function (z-zN)^T Q (z-zN) at value r to plot the
% %boundary of the confidence region. 
% minX=zN(1)-simCI(1);
% maxX=zN(1)+simCI(1);
% minY=zN(2)-simCI(2);
% maxY=zN(2)+simCI(2);
% 
% scale=1.1;
% 
% lobddX=round(zN(1)-scale*simCI(1),2,'significant');
% upbddX=round(zN(1)+scale*simCI(1),2,'significant');
% lobddY=round(zN(2)-scale*simCI(2),2,'significant');
% upbddY=round(zN(2)+scale*simCI(2),2,'significant');
% 
% numX=201;
% numY=201;
% 
% gapX=(upbddX-lobddX)/(numX-1);
% gapY=(upbddY-lobddY)/(numY-1);
% [X,Y] = meshgrid(lobddX:gapX:upbddX,lobddY:gapY:upbddY);
% 
% %X and Y are both matrices with numY rows and numX columns
% %X(1:numY,1)=lobddX, x(1:numY,numX)=upbddX
% %Y(1,1:numX)=lobddY, Y(numY,1:numX)=upbddY
% 
% Z=zeros(numY,numX);
% 
% for col=1:numX
%     for row=1:numY
%         xdiff=X(row,col)-zN(1);
%         ydiff=Y(row,col)-zN(2);
%         tmp=Q(1,1)*xdiff^2 +2* Q(1,2)* xdiff*ydiff+Q(2,2)*ydiff^2;
%         Z(row,col)=round(tmp,2,'significant');
%     end
%     
% end
% hold off;
% r_tmp=round(r,2,'significant');
% 
% [C,h]=contour(X,Y,Z,[r_tmp,r_tmp],'LineStyle','-','Color','k','LineWidth',1);
% %clabel(C,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot zN and z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on;
plot(zN(1),zN(2),'kx');
text(zN(1)+0.005,zN(2),'$z_N$','Interpreter','latex','FontSize',15);

if nargin==7
    plot(z0(1),z0(2),'ko','markerfacecolor','k');
    text(z0(1)+0.005,z0(2),'$z_0$','Interpreter','latex','FontSize',15);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graph simultaneous confidence intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minX=zN(1)-simCI(1);
maxX=zN(1)+simCI(1);
minY=zN(2)-simCI(2);
maxY=zN(2)+simCI(2);

line([minX,maxX],[maxY,maxY],'LineStyle','--','Color','k','LineWidth',1);
line([minX,maxX],[minY,minY],'LineStyle','--','Color','k','LineWidth',1);
line([minX,minX],[minY,maxY],'LineStyle','--','Color','k','LineWidth',1);
line([maxX,maxX],[minY,maxY],'LineStyle','--','Color','k','LineWidth',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graph individual confidence intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

line([tilde_zN(1)-indCI(1),tilde_zN(1)+indCI(1)],[tilde_zN(2),tilde_zN(2)],'LineStyle','-','Color','k','LineWidth',1);
line([tilde_zN(1),tilde_zN(1)],[tilde_zN(2)-indCI(2),tilde_zN(2)+indCI(2)],'LineStyle','-','Color','k','LineWidth',1);

plot(tilde_zN(1),tilde_zN(2),'kx');
text(tilde_zN(1)+0.005,tilde_zN(2),'$\tilde{z}_N$','Interpreter','latex','FontSize',15);

%add markers to end points of individual confidence intervals
if indCI(1)>0
    line([tilde_zN(1)-indCI(1),tilde_zN(1)-indCI(1)],[tilde_zN(2)-0.003,tilde_zN(2)+0.003],'LineStyle','-','Color','k');
    line([tilde_zN(1)+indCI(1),tilde_zN(1)+indCI(1)],[tilde_zN(2)-0.003,tilde_zN(2)+0.003],'LineStyle','-','Color','k');
end

if indCI(2)>0
    line([tilde_zN(1)-0.003,tilde_zN(1)+0.003],[tilde_zN(2)-indCI(2),tilde_zN(2)-indCI(2)],'LineStyle','-','Color','k');
    line([tilde_zN(1)-0.003,tilde_zN(1)+0.003],[tilde_zN(2)+indCI(2),tilde_zN(2)+indCI(2)],'LineStyle','-','Color','k');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


set(gca,'DataAspectRatio',[1 1 1],'fontsize',13);



hold on;


end