%%% tilt in y direction for 23.5*pi/180
tic;
clc
clear all
close all
t=0.001*pi;  
tilt=24*pi/180;
eux = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
euy = [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
euz = [cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];
euy_tilt = [cos(tilt) 0 sin(tilt);0 1 0;-sin(tilt) 0 cos(tilt)];
% eum = [cos(t) sin(t) 1;-sin(t) cos(t) 1;0 0 1];
radius = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b1 = 0.3;
b2 = 25;
step = 0.005;
bflies = input('enter no of butterflies:  ');
iter_max = input('enter no of iterations: ');
iter = 1;
swarm = zeros(bflies,6,iter);
for i=1:bflies
a=2*pi*rand;
r=sqrt(rand);
rc = 2;
swarm(i,1,iter)=(rc*r)*cos(a);
swarm(i,2,iter)=(rc*r)*sin(a);
    if mod(i,2)==0
        swarm(i,3,iter)=sqrt(4-(swarm(i,1,iter).^2)-(swarm(i,2,iter).^2));
    else
        swarm(i,3,iter)=-sqrt(4-(swarm(i,1,iter).^2)-(swarm(i,2,iter).^2));
    end
end
swarm(:,5,iter) = swarm(:,1,iter);
swarm(:,6,iter) = swarm(:,2,iter);
swarm(:,7,iter) = swarm(:,3,iter);
UV = ones(bflies,1).*10;
swarm(:,4,iter)= UV;
%%% main loop
for iter=2:iter_max
%     if iter==2
%         swarm(:,1,iter)=swarm(:,1,iter-1);
%         swarm(:,2,iter)=swarm(:,2,iter-1);
%         swarm(:,3,iter)=swarm(:,3,iter-1);
%     else
%         swarm(:,1,iter)=swarm(:,5,iter-1);
%         swarm(:,2,iter)=swarm(:,6,iter-1);
%         swarm(:,3,iter)=swarm(:,7,iter-1);
%     end
    %%%%% modified %%%%%%%%%%%%%%
    
        swarm(:,1,iter)=swarm(:,5,iter-1);
        swarm(:,2,iter)=swarm(:,6,iter-1);
        swarm(:,3,iter)=swarm(:,7,iter-1);
        for m=1:bflies
            s = euz*[swarm(m,1,iter);swarm(m,2,iter);swarm(m,3,iter)];
            swarm(m,1,iter) = s(1);
            swarm(m,2,iter) = s(2);
            swarm(m,3,iter) = s(3);
        end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:bflies
        for j=1:bflies
%             dist(i,j)=sqrt(sum((swarm(i,1:3,iter)-swarm(j,1:3,iter)).^2));
              k = [swarm(i,1,iter),swarm(i,1,iter),swarm(i,1,iter)];
              l = [swarm(j,1,iter),swarm(j,1,iter),swarm(j,1,iter)];
              dist(i,j) = distanc(k,l);
              heart(i,j) = dist(i,j);
        end
    end
    for i=1:bflies
        for j=1:bflies
            if i==j
                dist(i,j)=0;
            else
                dist(i,j)=1./dist(i,j);
            end
        end
    end
    dist = dist';
    g = sum(dist);
    g=g';
    dist = dist';
    for i=1:bflies
        for j=1:bflies
            r(i,j) = UV(i)*dist(i,j)./g(i);
        end
    end
    [B,I] = sort(r,'descend');
    z = swarm(:,3,iter);
    z=z';
    for i=1:bflies
        for j=1:bflies
            xy(i,j) = z(I(i,j));
        end
    end
for i=1:bflies
    flag = 0;
    for j=1:bflies
        if xy(i,j)>z(i)
            index(i) = I(i,j);
            flag = 1;
            break
        end
    end
    if flag==0
        index(i) = 0;
    end
end
t=0.001*pi; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eux = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
euy = [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
euz = [cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];


for i=1:bflies
       if index(i)~=0
            nom = sqrt(sum((swarm(index(i),1:3,iter)-swarm(i,1:3,iter)).^2));
            x_u = swarm(i,1,iter)+step*(swarm(index(i),1,iter)-swarm(i,1,iter))/nom;
            y_u = swarm(i,2,iter)+step*(swarm(index(i),2,iter)-swarm(i,2,iter))/nom;
            z_u = swarm(i,3,iter)+step*(swarm(index(i),3,iter)-swarm(i,3,iter))/nom;
            [az,el,r] = cart2sph(x_u,y_u,z_u);
            r = radius;
            [x_u,y_u,z_u] = sph2cart(az,el,r);
%             if iter~=2
%                rev_trans = inv(euy_tilt)*[x_u;y_u;z_u];
%                x_u = rev_trans(1);
%                y_u = rev_trans(2);
%                z_u = rev_trans(3);
%             end
%             h = euy_tilt*euz*[x_u;y_u;z_u]; %%%%%%%%%%% GAME CHANGER-1 %%%%%%%%%%%%%%
%             h = euz*[x_u;y_u;z_u];
%             swarm(i,5,iter)=h(1);
%             swarm(i,6,iter)=h(2);
%             swarm(i,7,iter)=h(3);
       else 
           nom=0;
%            if iter~=2
%                rev_trans = inv(euy_tilt)*[x_u;y_u;z_u];
%                x_u = rev_trans(1);
%                y_u = rev_trans(2);
%                z_u = rev_trans(3);
%             end
%             h = euy_tilt*euz*[swarm(i,1,iter);swarm(i,2,iter);swarm(i,3,iter)];%%%%%%%%%% GAME CHANGER-2 %%%%%%%%
%             h = euz*[swarm(i,1,iter);swarm(i,2,iter);swarm(i,3,iter)];        
%             swarm(i,1,iter) = h(1);
%             swarm(i,2,iter) = h(2);
%             swarm(i,3,iter) = h(3);
            swarm(i,5,iter)=swarm(i,1,iter);
            swarm(i,6,iter)=swarm(i,2,iter);
            swarm(i,7,iter)=swarm(i,3,iter);
       end
end
for i=1:bflies
    UV(i)=b1*UV(i)+b2*swarm(i,7,iter);
    swarm(i,4,iter)=UV(i);
end
fprintf('itaration...  %d\n',iter);
end
rho = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tetha = linspace(0,2*pi,30);
phi = linspace(-pi/2,pi/2,30);
[tetha phi] = meshgrid(tetha,phi);
[X Y Z] = sph2cart(tetha,phi,rho);
cdata = imread('image_file.jpg');
% t=0.001*pi; %%%%%%%%%%%%%%%%%%%%%%%%%%%
% eux = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
% euy = [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
% euz = [cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];
for i=1:iter_max
clf
for r=1:length(X)
for c=1:length(X)     
%      h = euy_tilt*euz*[X(r,c);Y(r,c);Z(r,c)];%%%%%%%%%%%%%%%%%%% GAME CHANGER-3 %%%%%%%%%%%%
     h = euz*[X(r,c);Y(r,c);Z(r,c)];
     X(r,c) = h(1);
     Y(r,c) = h(2);
     Z(r,c) = h(3);

end
end
% globe = surf(X,Y,Z);
% set(globe,'FaceColor','texturemap','Cdata',cdata);
    plot3(X,Y,Z);
   zlim([-3 3]);
   ylim([-3 3]);
   zlim([-3 3]);
   grid on
   hold on

% plot3(swarm(:,1,i),swarm(:,2,i),swarm(:,3,i),'r*');
plot3(swarm(:,5,i),swarm(:,6,i),swarm(:,7,i),'r*');
   title(strcat('Iteration: ',num2str(i)))
   axis([-3 3 -3 3]);
   pause(0.00001);
   
%    for r=1:length(X)
%     for c=1:length(X)
%      h = inv(euy_tilt)*[X(r,c);Y(r,c);Z(r,c)];
%      X(r,c) = h(1);
%      Y(r,c) = h(2);
%      Z(r,c) = h(3);
%     end
%    end

end
% figure(2)
% for i=1:iter_max
%     clf
%     t=0.001*pi; %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     eux = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
%     euy = [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
%     euz = [cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];

%     for r=1:length(X)
%     for c=1:length(X)     
%         h = euz*[X(r,c);Y(r,c);Z(r,c)];
%         X(r,c) = h(1);
%         Y(r,c) = h(2);
%         Z(r,c) = h(3);
%    end
%    end
%    contour(X,Y,Z,15)
%    grid on
%    hold on
%    plot3(swarm(:,1,i),swarm(:,2,i),swarm(:,3,i),'r>');
%    title(strcat('Iteration: ',num2str(i)));
%    axis([-3 3 -3 3 -3 3]);
%    xlabel('x-axis');
%    ylabel('y-axis');
%    zlabel('z-axis');
%    pause(0.00001)
% end
% figure(2)
% for i=1:bflies
%    contour(X,Y,Z,15)
%    hold on
%    x=swarm(i,1,:);
%    y=swarm(i,2,:);
%    z=swarm(i,3,:);
%    plot3(x(:),y(:),z(:),'m-',swarm(i,1,iter_max),swarm(i,2,iter_max),swarm(i,3,iter_max),'r*',swarm(i,1,1),swarm(i,2,1),swarm(i,3,1),'r*');
%    axis([-3 3 -3 3 -3 3]);
% end
% grid on
toc;
        
plot3((swarm(:,5,10)),(swarm(:,5,10)),(swarm(:,5,10)),'r*')

