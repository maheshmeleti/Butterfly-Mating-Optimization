clc
clear all
close all
bflies = input('enter no of bfies: ');
iterations = input('enter no of iterations: ');


%%% rotation functions %%%%
t = 0.001*pi;
tilt = 24*pi/180;
eux = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
euy = [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
euz = [cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];
euy_tilt = [cos(tilt) 0 sin(tilt);0 1 0;-sin(tilt) 0 cos(tilt)];
%%%


radius = 3;
b1 = 0.3;
b2 = 25;
step = 0.009;
iter = 1;
swarm = zeros(bflies,5,iter);

for i=1:bflies
a=2*pi*rand;
r=sqrt(rand);
rc = radius;
swarm(i,1,iter)=(rc*r)*cos(a);
swarm(i,2,iter)=(rc*r)*sin(a);
    if mod(i,2)==0
        swarm(i,3,iter)=sqrt((radius^2)-(swarm(i,1,iter).^2)-(swarm(i,2,iter).^2));
    else
        swarm(i,3,iter)=-sqrt((radius^2)-(swarm(i,1,iter).^2)-(swarm(i,2,iter).^2));
    end
end


UV = ones(bflies,1).*10;
swarm(:,4,iter) = UV;

for iter = 2:iterations
    
    for c=1:bflies
        if iter==2
            k = euy_tilt*euz*[swarm(c,1,iter-1);swarm(c,2,iter-1);swarm(c,3,iter-1)];
        else
            T =  inv(euy_tilt)*[swarm(c,1,iter-1);swarm(c,2,iter-1);swarm(c,3,iter-1)];
            swarm(c,1,iter-1) = T(1);
            swarm(c,2,iter-1) = T(2);
            swarm(c,3,iter-1) = T(3);
            k = euy_tilt*euz*[swarm(c,1,iter-1);swarm(c,2,iter-1);swarm(c,3,iter-1)];
        end
        
        swarm(c,1,iter-1) = k(1);
        swarm(c,2,iter-1) = k(2);
        swarm(c,3,iter-1) = k(3);
    end
    swarm(:,1,iter) = swarm(:,1,iter-1);
    swarm(:,2,iter) = swarm(:,2,iter-1);
    swarm(:,3,iter) = swarm(:,3,iter-1);

    for i=1:bflies
        for j=i+1:bflies
%             fprintf('\n%d %d',i,j);
            dist(i,j) = sqrt(sum((swarm(i,1:3,iter)-swarm(j,1:3,iter)).^2));
            dist(j,i) = dist(i,j);
              
        end
    end
    
    for i=1:bflies
        for j=1:bflies
            if i==j
                dist(i,j)=0;
            else
                dist(i,j) = 1./dist(i,j);
            end
        end
    end
    
    dist = dist';
    g = sum(dist);
    g = g';
    dist = dist';
    
    for i=1:bflies
        for j=1:bflies
            r(i,j) = UV(i)*dist(i,j)/g(i);
        end
    end
    
    [B,I] = sort(r,'descend');
    z = swarm(:,3,iter);
    z = z';
    
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
    
    for i=1:bflies
        if index(i)~=0
            nom = sqrt(sum((swarm(index(i),1:3,iter)-swarm(i,1:3,iter)).^2));
            swarm(i,1,iter) = swarm(i,1,iter)+step*(swarm(index(i),1,iter)-swarm(i,1,iter))/nom;
            swarm(i,2,iter) = swarm(i,2,iter)+step*(swarm(index(i),2,iter)-swarm(i,2,iter))/nom;
            swarm(i,3,iter) = swarm(i,3,iter)+step*(swarm(index(i),3,iter)-swarm(i,3,iter))/nom;
            [az,el,r] = cart2sph(swarm(i,1,iter),swarm(i,2,iter),swarm(i,3,iter));
            r = radius;
            [swarm(i,1,iter),swarm(i,2,iter),swarm(i,3,iter)] = sph2cart(az,el,r);
    
        end
%         if ((swarm(i,1,iter)^2)+(swarm(i,2,iter)^2)+(swarm(i,3,iter)^2))~=radius^2
%             fprintf('BFLY %d error %d \n',i,(swarm(i,1,iter).^2)+(swarm(i,2,iter).^2)+(swarm(i,3,iter).^2));
%         end
     end
    %%%(swarm(5,1,10)^2)+(swarm(5,2,10)^2)+(swarm(5,3,10)^2)
    for i=1:bflies
    UV(i)=b1*UV(i)+b2*swarm(i,3,iter);
    swarm(i,4,iter)=UV(i);
    end
    
    fprintf('iteration %d...\n',iter);
end
% figure(1);
%%%%%%%%%%%% sphere %%%%%%%%%%%%
cdata = imread('image_file.jpg');
tetha = linspace(0,2*pi,50);
phi =  linspace(-pi./2,pi./2,50);
[tetha phi] = meshgrid(tetha,phi);
[x,y,z] = sph2cart(tetha,phi,radius);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = input('enter 1: ');
if k==1;
for iter=1:iterations
    clf
    for i=1:length(x)
        for j=1:length(x)
            k = euy_tilt*euz*[x(i,j);y(i,j);z(i,j)];
            x(i,j) = k(1);
            y(i,j) = k(2);
            z(i,j) = k(3);
        end
    end
%     plot3(x,y,z);
    globe = surf(x,y,-z);
    set(globe,'FaceColor','texturemap','Cdata',cdata);
    hold on
    plot3(swarm(:,1,iter),swarm(:,2,iter),swarm(:,3,iter),'r*');
    title(strcat('Iteration: ',num2str(iter)))
    grid on
    pause(0.00001);
    for i=1:length(x)
        for j=1:length(x)
            k = inv(euy_tilt)*[x(i,j);y(i,j);z(i,j)];
            x(i,j) = k(1);
            y(i,j) = k(2);
            z(i,j) = k(3);
        end
    end
end
end