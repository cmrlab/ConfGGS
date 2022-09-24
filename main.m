%% Main function for ConfGGS algorithm towards optimization
function iteration=ConfGGS(dim,sigma)
tic
clc
N=50;
max_it=200;
ElitistCheck=1; Rpower=1;
min_flag=1; % 1: minimization, 0: maximization

[iteration,Fbest,Lbest,BestChart,MeanChart]=GGSA(N,max_it,ElitistCheck,min_flag,Rpower,dim,sigma);

% plot(BestChart,'--k');
% xlabel('\fontsize{12}\bf Iteration');
% ylabel('\fontsize{12}\bf Evaluated value'); 
% legend('\fontsize{10}\bf GGS',1);
toc
positions=csvread('Output.csv');
n=size(positions,1); % num. Of Points
cod1=csvread('cod_bond.csv'); % read the cod_bond
cod2=csvread('cod_angle.csv'); % read the cod_angle
cod3=csvread('cod_dihedral.csv'); % read the cod_dihedral
m1=size(cod1,1);
m2=size(cod2,1);
m3=size(cod3,1);
cosfi=zeros(1,m3);

for i=1:m1
      dist(i)=findDist(positions(cod1(i,1),:),positions(cod1(i,2),:));
end

dist'
for i=1:m2
      angle(i)=findAngle(positions(cod2(i,1),:),positions(cod2(i,2),:),positions(cod2(i,3),:));
end

angle'
for i=1:m3
     cosfi(i)=findDihedral(positions(cod3(i,1),:),positions(cod3(i,2),:),positions(cod3(i,3),:),positions(cod3(i,4),:));
%      acosfi(i)=acos(cosfi(i))*180/3.141;
     
end

cosfi'

csvwrite('Output_bond.csv',dist');
csvwrite('Output_angle.csv',angle');
csvwrite('Output_dihedral.csv',cosfi');
end

%% first call from main
function [iteration,Fbest,Lbest,BestChart,MeanChart]=GGSA(N,max_it,ElitistCheck,min_flag,Rpower,dim,sigma)
Rnorm=2;
%get allowable range and dimension of the test function.
%[low,up,dim]=test_functions_range(F_index);
low=-8;
up=8;

%random initialization for agents.
X=initialization(dim,N,up,low);
%create the best so far chart and average fitnesses chart.
BestChart=[];
MeanChart=[];

V=zeros(N,dim);

for iteration=1:max_it
    %Checking allowable range.
    X=space_bound(X,up,low);
    
    %Evaluation of agents.
    fitness=evaluateF(X);
    
    if min_flag==1
        [best best_X]=min(fitness); %minimization.
    else
        [best best_X]=max(fitness); %maximization.
    end
    
    if iteration==1
        Fbest=best;
        Lbest=X(best_X,:);
    end
    if min_flag==1
        if best<Fbest  %minimization.
            Fbest=best;
            Lbest=X(best_X,:);
        end
    else
        if best>Fbest  %maximization
            Fbest=best;
            Lbest=X(best_X,:);
        end
    end
    
    BestChart=[BestChart Fbest];
    s=sprintf('%10.3f',Fbest);
    fprintf('%d\t%10.3f\n',iteration,Fbest);

    MeanChart=[MeanChart mean(fitness)];
    
    %Calculation of M.
    [M]=massCalculation(fitness,min_flag);
    
    %Calculation of Gravitational constant.
    G=Gconstant(iteration,max_it);
    
    %Calculation of accelaration in gravitational field. eq.7-10,21.
    a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it,sigma);
    if(Fbest<0.009)
        break; 
    end
    %Agent movement.
    [X,V]=move(X,a,V);
    chromosomes(iteration,1:30)=Lbest;
end %iteration
ppp=30;
ref=csvread('ref.csv');
[r c]=size(ref);
for i=1:r
    reff(1,(i-1)*3+1:i*3)=ref(i,1:3);
end
for i=1:(iteration-1)
    [a, b, c, d,e]=finderror(chromosomes(i,1:ppp),3,10,reff);
    chromosomes(i,1+ppp)=a;
    chromosomes(i,2+ppp)=b;
    chromosomes(i,3+ppp)=c;
    chromosomes(i,4+ppp)=d; 
    chromosomes(i,5+ppp)=e;
end
csvwrite('Best_Coordinates.csv',chromosomes);
optimal_Coordinate=chromosomes(iteration-1,1:ppp);
for i=1:r
    output(i,1:3)=optimal_Coordinate(1,(i-1)*3+1:i*3);
end
csvwrite('Output.csv',output);
% plot(chromosomes(:,1+ppp)); % Plot Total Energy
% plot(chromosomes(:,2+ppp)); % Plot Bonded Energy
% plot(chromosomes(:,3+ppp)); % Plot Non-Bonded Energy
% plot(chromosomes(:,4+ppp)); % Plot Dihedral Energy   
plot(chromosomes(:,5+ppp)); % Plot Fitness
end

%% This function initializes the position of the agents in the search space, randomly.
function [X]=initialization(dim,N,up,down)

% if size(up,2)==1
    X=rand(N,dim).*(up-down)+down;
     PPP1=csvread('Input.csv');
    [r c]=size(PPP1);
    for i=1:r
        PPP(1,(i-1)*3+1:i*3)=PPP1(i,1:3);
    end
   
    X(N,:)=PPP';

end

%% This function checks the search space boundaries for agents.
function  X=space_bound(X,up,low)
[N,dim]=size(X);
for i=1:N
    %     %%Agents that go out of the search space, are reinitialized randomly .
    Tp=X(i,:)>up;
    Tm=X(i,:)<low;
    X(i,:)=(X(i,:).*(~(Tp+Tm)))+((rand(1,dim).*(up-low)+low).*(Tp+Tm));
end
end                                    

%% This function Evaluates the agents.
function   fitness=evaluateF(X)
[N,dim]=size(X);
ref=csvread('ref.csv');
[r c]=size(ref);
for i=1:r
    reff(1,(i-1)*3+1:i*3)=ref(i,1:3);
end
for i=1:N
    [a,b,c,d,e]=finderror(X(i,:),3,10,reff);
    fitness(i)=e;
end
end

%% This function calculates the mass of each agent. eq.14-20
function [M]=massCalculation(fit,min_flag)
%we can also change the masscalculation equation based on our requirement
Fmax=max(fit);
Fmin=min(fit);
Fmean=mean(fit);
[~, N]=size(fit);

if Fmax==Fmin
    M=ones(N,1);
else
    
    if min_flag==1 %for minimization
        best=Fmin;
        worst=Fmax;
    else %for maximization
        best=Fmax;
        worst=Fmin;
    end
    
    M=(fit-worst)./(best-worst);
    
end

M=M./sum(M);
end

%% This function calculates Gravitational constant.
function G=Gconstant(iteration,max_it)
%%%here, make your own function of 'G'
alfa=20;
G0=100;
G=G0*exp(-alfa*iteration/max_it); %Exponentially decreasing over iterations

end

%% Calculation of Gradient G-field
function a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it,sigma)

[N,dim]=size(X);
final_per=2; %In the last iteration, only 2 percent of agents apply force to the others.

%total force calculation
if ElitistCheck==1
    kbest=final_per+(1-iteration/max_it)*(100-final_per); %kbest in eq. 21.
    kbest=round(N*kbest/100);
else
    kbest=N;
end
[~, ds]=sort(M,'descend');

for i=1:N
    E(i,:)=zeros(1,dim);
    for ii=1:kbest
        j=ds(ii);
        if j~=i
            R=norm(X(i,:)-X(j,:),Rnorm); %Euclidian distanse.
            for k=1:dim
                E(i,k)=E(i,k)+rand*(M(j))*((X(j,k)-X(i,k))/(R^Rpower+eps));
                %note that Mp(i)/Mi(i)=1
            end
        end
    end
    min1=min(min(gradient(E(i,:))));
    max1=max(max(gradient(E(i,:))));
    E(i,:)=E(i,:)+sigma*(max1-min1);
end

%%compute acceleration
a=E.*G; %note that Mp(i)/Mi(i)=1
end

%% This function updates the velocity and position of agents.
function [X,V]=move(X,a,V)
%movement.
[N,dim]=size(X);
V=rand(N,dim).*V + a;
X=X + V;
end