%% 3 Class Fuzzy K-means part 1
clear all;clc;close all
% Define the classes
class=[-7.82,-4.58,-3.97;-6.68,3.16,2.71;4.36,-2.19,2.09;6.72,0.88,2.80;-8.64,3.06,3.50;-6.87,0.57,-5.45;4.47,-2.62,5.76;6.73,-2.01,4.18;-7.71,2.34,-6.33;-6.91,-0.49,-5.68;6.18,2.81,5.82;6.72,-0.93,-4.04;-6.25,-0.26,0.56;-6.94,-1.22,1.13;8.09,0.20,2.25;6.81,0.17,-4.15;-5.19,4.24,4.04;-6.38,-1.74,1.43;4.08,1.30,5.33;6.27,0.93,-2.78];

% intitalize membership matrix of random values
U=[];
for i=1:20
    U(i,:)=rand(1,3);    
end

% Define FUZZIFIER parameter
q=2;
% Define matrix m
m=[0,0,0;1,1,1;-1,0,2];
% Initialize counter n and learning rate E
n=1;
E=0.001;
% store m for comparison
prev_m=m;
while n<100  
    cluster1=[];cluster2=[];cluster3=[];
    n=n+1;
    % Normalize matrix U
    for i=1:length(U)
        U(i,:)=U(i,:)./sum(U(i,:));
    end
    % Update Matrix m
    for i=1:3
        m(i,1:3)=sum((U(:,i).^q).*class)./sum(U(:,i).^2);
    end
    % Euclidean distance to find the cluster 
    dist1=[];dist2=[];dist3=[];
    for i=1:length(class)
        dist1=[dist1;sqrt(sum((class(i,:)-m(1,:)).^2))];
        dist2=[dist2;sqrt(sum((class(i,:)-m(2,:)).^2))];
        dist3=[dist3;sqrt(sum((class(i,:)-m(3,:)).^2))];
        if dist1(i)<dist2(i) & dist1(i)<dist3(i);
            cluster1=[cluster1;class(i,:)];
        elseif dist2(i)<dist1(i) & dist2(i)<dist3(i);
            cluster2=[cluster2;class(i,:)];
        else dist3(i)<dist1(i) & dist3(i)<dist2(i);
            cluster3=[cluster3;class(i,:)];
        end   
    end
    % Update membership matrix U
    for i=1:20
        U(i,1)=(inv(dist1(i)^2))/(inv(dist1(i))^2+inv(dist2(i))^2+inv(dist3(i))^2);
        U(i,2)=(inv(dist2(i)^2))/(inv(dist1(i))^2+inv(dist2(i))^2+inv(dist3(i))^2);
        U(i,3)=(inv(dist3(i)^2))/(inv(dist1(i))^2+inv(dist2(i))^2+inv(dist3(i))^2);
    end
    % Check if value of m has changed
    if (prev_m(1,:)-m(1,:)<E) & (prev_m(2,:)-m(2,:)<E) & (prev_m(3,:)-m(3,:)<E);
        break;
    end
    prev_m=m;    
end

% Display Number of iterations required
X=['No of iterations = ',num2str(n)];
disp(X) 

% Plot
figure,scatter3(class(:,1),class(:,2),class(:,3))
figure,scatter3(cluster1(:,1),cluster1(:,2),cluster1(:,3),'MarkerFaceColor',[0 0 1])
hold on
scatter3(cluster2(:,1),cluster2(:,2),cluster2(:,3),'MarkerFaceColor',[1 0 0])
scatter3(cluster3(:,1),cluster3(:,2),cluster3(:,3),'MarkerFaceColor',[1 1 0])
scatter3(m(1,1),m(1,2),m(1,3),'d','MarkerFaceColor',[0 0 1])
scatter3(m(2,1),m(2,2),m(2,3),'d','filled','MarkerFaceColor',[1 0 0])
scatter3(m(3,1),m(3,2),m(3,3),'d','MarkerFaceColor',[1 1 0])
hold off
legend('cluster1','cluster2','cluster3','mean1','mean2','mean3');

%% 3 Class Fuzzy K-means part 1
clear all;clc;close all

% Define the classes
class=[-7.82,-4.58,-3.97;-6.68,3.16,2.71;4.36,-2.19,2.09;6.72,0.88,2.80;-8.64,3.06,3.50;-6.87,0.57,-5.45;4.47,-2.62,5.76;6.73,-2.01,4.18;-7.71,2.34,-6.33;-6.91,-0.49,-5.68;6.18,2.81,5.82;6.72,-0.93,-4.04;-6.25,-0.26,0.56;-6.94,-1.22,1.13;8.09,0.20,2.25;6.81,0.17,-4.15;-5.19,4.24,4.04;-6.38,-1.74,1.43;4.08,1.30,5.33;6.27,0.93,-2.78];

% intitalize membership matrix of random values
U=[];
for i=1:20
    U(i,:)=rand(1,3);    
end

% Define FUZZIFIER parameter
q=2;
% Define matrix m
m=[-0.1,0,0.1;0,-0.1,0.1;-0.1,-0.1,0.1];
% Initialize counter n and learning rate E
n=1;
E=0.001;
% store m for comparison
prev_m=m;
while n<100  
    cluster1=[];cluster2=[];cluster3=[];
    n=n+1;
    % Normalize matrix U
    for i=1:length(U)
        U(i,:)=U(i,:)./sum(U(i,:));
    end
    % Update Matrix m  
    for i=1:3
        m(i,1:3)=sum((U(:,i).^q).*class)./sum(U(:,i).^2);
    end
    % Euclidean distance to find the cluster    
    dist1=[];dist2=[];dist3=[];
    for i=1:length(class)
        dist1=[dist1;sqrt(sum((class(i,:)-m(1,:)).^2))];
        dist2=[dist2;sqrt(sum((class(i,:)-m(2,:)).^2))];
        dist3=[dist3;sqrt(sum((class(i,:)-m(3,:)).^2))];
        if dist1(i)<dist2(i) & dist1(i)<dist3(i);
            cluster1=[cluster1;class(i,:)];
        elseif dist2(i)<dist1(i) & dist2(i)<dist3(i);
            cluster2=[cluster2;class(i,:)];
        else dist3(i)<dist1(i) & dist3(i)<dist2(i);
            cluster3=[cluster3;class(i,:)];
        end   
    end
    % Update membership matrix U
    for i=1:20
        U(i,1)=(inv(dist1(i)^2))/(inv(dist1(i))^2+inv(dist2(i))^2+inv(dist3(i))^2);
        U(i,2)=(inv(dist2(i)^2))/(inv(dist1(i))^2+inv(dist2(i))^2+inv(dist3(i))^2);
        U(i,3)=(inv(dist3(i)^2))/(inv(dist1(i))^2+inv(dist2(i))^2+inv(dist3(i))^2);
    end
    % Check if value of m has changed
    if (prev_m(1,:)-m(1,:)<E) & (prev_m(2,:)-m(2,:)<E) & (prev_m(3,:)-m(3,:)<E);
        break;
    end
    prev_m=m;    
end

% Display Number of iterations required
X=['No of iterations = ',num2str(n)];
disp(X) 

%Plot
figure,scatter3(class(:,1),class(:,2),class(:,3))
figure,scatter3(cluster1(:,1),cluster1(:,2),cluster1(:,3),'MarkerFaceColor',[0 0 1])
hold on
scatter3(cluster2(:,1),cluster2(:,2),cluster2(:,3),'MarkerFaceColor',[1 0 0])
scatter3(cluster3(:,1),cluster3(:,2),cluster3(:,3),'MarkerFaceColor',[1 1 0])
scatter3(m(1,1),m(1,2),m(1,3),'d','MarkerFaceColor',[0 0 1])
scatter3(m(2,1),m(2,2),m(2,3),'d','filled','MarkerFaceColor',[1 0 0])
scatter3(m(3,1),m(3,2),m(3,3),'d','MarkerFaceColor',[1 1 0])
hold off
legend('cluster1','cluster2','cluster3','mean1','mean2','mean3');
