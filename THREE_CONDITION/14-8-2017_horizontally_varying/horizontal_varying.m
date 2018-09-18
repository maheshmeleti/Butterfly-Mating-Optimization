clc
clear all
close all
b1=0.3;
b2=25;
step=0.007;
bflies=100;                                                                       
iter_max=1000; %  Change the value of iter_max so that we can observe the variation.
iter=1;
access = linspace(-3,3,iter_max);
swarm=zeros(bflies,6,iter); 
% [X Y] = meshgrid(-3:0.1:3,-3:0.1:3);
% Z=fitness(X,Y,tetha); 
 for i=1:bflies
swarm(i,1,iter)=3*(rand-0.5)/0.5;
swarm(i,2,iter)=3*(rand-0.5)/0.5;
 end
z = fitness(swarm(:,1,iter),swarm(:,2,iter),access(1)); %Fitness calculating from 3-peaks function
swarm(:,3,iter)=z;
UV = ones(bflies,1)*10;
swarm(:,4,iter)=UV;
% disp(swarm);
% fprintf('\n');
x=swarm(1:bflies,1,1);
y=swarm(1:bflies,2,1);
z=fitness(x,y,access(1));
for iter=2:iter_max
% tetha = tetha+(5*pi/(2.*iter_max));
%  fprintf('tetha %d\n',tetha);
     if iter==2
         swarm(:,1,iter)=swarm(:,1,iter-1);
         swarm(:,2,iter)=swarm(:,2,iter-1);
   else
       swarm(:,1,iter)=swarm(:,5,iter-1);
       swarm(:,2,iter)=swarm(:,6,iter-1);
     end
     z = fitness(swarm(:,1,iter),swarm(:,2,iter),access(iter));
 %Fitness calculating from 3-peaks function
swarm(:,3,iter)=z;
% randomly initialize bee with their positions
for i=1:bflies
        for j=i+1:bflies
            dist(i,j)= sqrt(sum((swarm(i,1:2,iter)-swarm(j,1:2,iter)).^2));%Distance between each bfly to other bflies
            dist(j,i)=dist(i,j);
        end
end
   for i=1:bflies
    for j=1:bflies
        if i==j
            dista(i,j)=0;
        else
dista(i,j)=1./dist(i,j); % Finding the inverse UV for the sake of calculating the UV distribution
        end
    end
   end
% fprintf('\ndistanace\n');
% disp(dista);
    g=sum(dista);
    
    for i=1:bflies
        j=1:bflies;
        r(i,j)=(dista(:,i)*UV(i,1))./(g(1,i)); % Measuring UV distribution
    end
    [B,I]=sort(r,'descend');% Sorting the distributed UV values
%     fprintf('\n value of B-\n');
%     disp(B);
%     fprintf('\n value of I-\n');
%     disp(I);

for i=1:bflies
    for j=1:bflies
        xy(i,j)=z(I(i,j)); %% comparing with the fitness of other bflies
    end
end
% fprintf('xy\n');
%     disp(xy);
% fprintf('\n');
   for i=1:bflies           %Bfly number
       flag=0;
       for j=1:bflies       %row values in bfly column
        if (xy(j,i)> z(i))
            %m(i,1)=z(i);
            index(i,1)=I(j,i);
            flag=1;
            break;
        %else
                %m(i,1)=0;
        end
       end
       if(flag==0)
           index(i,1)=0; % This gives the l-mate of the each bfly in the row wise order.
       end
       b(i)=index(i);
       lmate=transpose(b);
   end
%    fprintf('\n value of l-mate-\n');
%     disp(lmate);
%     fprintf('\n');
   for i=1:bflies
       if lmate(i)
%            fprintf('\n%d',lmate(i));
      nom=(sqrt((swarm(lmate(i),1,iter)-swarm(i,1,iter))^2 + (swarm(lmate(i),2,iter)-swarm(i,2,iter))^2));
   swarm(i,5,iter)=swarm(i,1,iter)+step*((swarm(lmate(i),1,iter)-swarm(i,1,iter)))/nom;
   swarm(i,6,iter)=swarm(i,2,iter)+step*((swarm(lmate(i),2,iter)-swarm(i,2,iter)))/nom;
     
       else 
           nom=0;
            swarm(i,5,iter)=swarm(i,1,iter);
            swarm(i,6,iter)=swarm(i,2,iter);
       end
   end
  
%   %UV updation
  for i=1:bflies
         UV(i)=b1*UV(i)+b2*fitness(swarm(i,5,iter),swarm(i,6,iter),access(iter));
%          swarm(i,7,iter)=UV(i);
           swarm(i,4,iter)=UV(i);
%     if round(UV(i))==round(UV(lmate(i)))
%      swarm(i,7,iter:iter_max)=UV(i);
%     end
  end
  fprintf('iteration %d\n',iter);
%   disp(swarm);
end
% [X Y] = meshgrid(-3:0.1:3,-3:0.1:3);
% Z=fitness(X,Y,tetha);
%  figure(1)
figure('Renderer','zbuffer');
% zlim([-20 20]);
% tetha = 0;
nframes = iter_max;
mo(1:iter_max)=struct('cdata',[],...
    'colormap',[]);
set(gca,'nextplot','replacechildren');
for i=1:iter_max
    fprintf('-iterationfor plotting %d\n-----',i);
   clf
%    zlim([-20 20]);
   [X Y] = meshgrid(-3:0.1:3,-3:0.1:3);
   Z=fitness(X,Y,access(i));
   surf(X,Y,Z);
   zlim([-8 8]);
   grid on
   hold on
%    plot3(swarm(:,1,i),swarm(:,2,i),fitness(swarm(:,1,i),swarm(:,2,i),tetha),'r*');
plot3(swarm(:,1,i),swarm(:,2,i),swarm(:,3,i),'b*');
   title(strcat('Iteration: ',num2str(i)))
   axis([-3 3 -3 3]);
   mov(i) = getframe(gcf);
   pause(0.00001);
end
movie2avi(mov,'mahesh.avi','compression','None');
% figure(2)
% for i=1:iter_max
%    clf
%    contour(X,Y,Z,15);
%    grid on
%    hold on
%    plot(swarm(:,1,i),swarm(:,2,i),'r*');
%    title(strcat('Iteration: ',num2str(i)));
%    axis([-3 3 -3 3]);
%    pause(0.00001)
% end
% figure(3)
% for i=1:bflies
%    contour(X,Y,Z,15);
%    hold on
%    x=swarm(i,1,:);
%    y=swarm(i,2,:);
%    plot(x(:),y(:),'-b',swarm(i,1,1),swarm(i,2,1),'r*');
%    axis([-3 3 -3 3]);
%    %axis tight
% end
% grid on
  
