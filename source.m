% Building Data
clc
clear
%load Brain1dose.mat;
%load DTruncatedWorkspace.mat;
%load matlab.mat;
load LP1Workspace.mat
%load LP3dose.mat
[oar, nbeam]= DefineStructure(cst,dij);

%clear cst
%clear ct
%clear dij
%clear pln
%clear stf

disp("Data loaded!");
%% 
prolongoar= IncrResInterior(oar,ct);

%%
% Preprocessing data
%oar=prolongoar;
scal=1;
theta=10^(-0.5);
tol=[10^(-5),10^(-5),10^(-5),10^(-5)];
method=4;
sparseoar = oar;
tic
n1=0;
for o = 1:4
    fprintf('OR%d\n',o);
    fprintf('Nonzeros before sparsification %d\n',nnz(oar(o).doseinfluence));
    [tr, tc, tv]=find(max(max(oar(o).doseinfluence)));
    fprintf('Maximum: %d\n',tv(1));
    [tr, tc, tv]=find(min(oar(o).doseinfluence(oar(o).doseinfluence>0)));
    fprintf('Minimum nonzero: %d\n',tv);
    [row, col, val]=find(oar(o).doseinfluence); 
    [m,n]=size(oar(o).doseinfluence);
    if method==1
        maxC=zeros(n,1);
        for j =1:size(val)
            if val(j)>maxC(col(j))
                maxC(col(j))=val(j);
            end
        end
        SS = 10;
        jStep = round(size(val,1) / SS);
        
        for i =1:size(val)
            if mod(i, jStep) == 0
                fprintf('- processing %d percent...\n', SS * i / jStep);
            end
            if val(i)<theta * maxC(col(i))
                 n1=n1+val(i);
                 val(i)=0;
            end
        end
    elseif method==2
        maxC=zeros(n,1);
        for j =1:size(val)
            if val(j)>maxC(col(j))
                maxC(col(j))=val(j);
            end
        end
        SS = 10;
        jStep = round(size(val,1) / SS);
        
        for i =1:size(val)
            if mod(i, jStep) == 0
                fprintf('- processing %d percent...\n', SS * i / jStep);
            end
            if val(i)<theta * maxC(col(i)) || maxC(col(i))<=tol(o)
                 n1=n1+val(i);
                 val(i)=0;
            end
        end
      elseif method==3
        maxC=zeros(m,1);
        for j =1:size(val)
            if val(j)>maxC(row(j))
                maxC(row(j))=val(j);
            end
        end
        SS = 10;
        jStep = round(size(val,1) / SS);
        
        for i =1:size(val)
            if mod(i, jStep) == 0
                fprintf('- processing %d percent...\n', SS * i / jStep);
            end
            if val(i)<theta * maxC(row(i)) %|| maxC(row(i))<=tol(o)
                 n1=n1+val(i);
                 val(i)=0;
            end
        end
      elseif method==4
        maxC=max(val);
        SS = 10;
        jStep = round(size(val,1) / SS);
        
        for i =1:size(val)
            if mod(i, jStep) == 0
                fprintf('- processing %d percent...\n', SS * i / jStep);
            end
            if val(i)<theta * maxC 
                 n1=n1+val(i);
                 val(i)=0;
            end
        end
      elseif method==5
        SS = 10;
        jStep = round(size(val,1) / SS);
        
        for i =1:size(val)
            if mod(i, jStep) == 0
                fprintf('- processing %d percent...\n', SS * i / jStep);
            end
            if val(i)<theta  
                 n1=n1+val(i);
                 val(i)=0;
            end
        end
      


    end
    
    sparseoar(o).doseinfluence=sparse(row, col, scal*val,m,n);
    fprintf('Nonzeros after sparsification: %d\n',nnz(sparseoar(o).doseinfluence));
    fprintf('Sparsification ratio: %d\n',nnz(sparseoar(o).doseinfluence)/nnz(oar(o).doseinfluence));
    [tr, tc, tv]=find(max(max(sparseoar(o).doseinfluence)));
    fprintf('Maximum: %d\n',tv(1));
    [tr, tc, tv]=find(min(sparseoar(o).doseinfluence(sparseoar(o).doseinfluence>0)));
    fprintf('Minimum nonzero: %d\n',tv);
end
fprintf('Norm one of all modifications: %d\n',n1);
disp("Data preprocessed!")
toc
%%
% Assume A is your sparse matrix
A = sparseoar(4).doseinfluence;
% Find the row indices, column indices, and values of the non-zero elements
[i, j, v] = find(A);
i = int32(i);
j = int32(j);
    
    % Combine data into a matrix
data = [double(i-1), double(j-1), v];

% Combine these into a single matrix where the first column is row index,
% the second column is column index, and the third column is the value


% Write the data to a CSV file
%csvwrite('LP1Target.csv', data);
writematrix(data, 'Lungs1.csv');
%%
names = {'Lungs', 'Heart', 'Esophagus', 'Target'};
for t = 1:4
    A = sparseoar(t).doseinfluence;
    [i, j, v] = find(A); % Extract non-zero elements
    
    % Ensure indices are integers
    i = int32(i);
    j = int32(j);
    
    % Combine data into a matrix
    data = [double(i-1), double(j-1), v]; % Convert indices to double for CSV compatibility
    
    % Write data to CSV file with exact precision
    writematrix(data, sprintf('%s1.csv', names{t}));
end

%%
% Nonlinear Model
A1=[sum(sparseoar(1).doseinfluence,1) sparse(1,2*size(sparseoar(4).voilist,1))];
A2=[sum(sparseoar(2).doseinfluence,1) sparse(1,2*size(sparseoar(4).voilist,1))];
A3=[sum(sparseoar(3).doseinfluence,1) sparse(1,2*size(sparseoar(4).voilist,1))];
A4=[sparseoar(4).doseinfluence -speye(size(sparseoar(4).voilist,1)) sparse(size(sparseoar(4).voilist,1),size(sparseoar(4).voilist,1))];
A5=[-sparseoar(4).doseinfluence sparse(size(sparseoar(4).voilist,1),size(sparseoar(4).voilist,1)) -speye(size(sparseoar(4).voilist,1))];
A=sparse([A1;A2;A3;A4;A5]);
b1=[(min(sparseoar(1).meandose, sparseoar(1).firstmom))*size(sparseoar(1).voilist,1)];
b2=[(min(sparseoar(2).meandose, sparseoar(2).firstmom))*size(sparseoar(2).voilist,1)];
b3=[(sparseoar(3).meandose)*size(sparseoar(3).voilist,1)];
b4= sparseoar(4).maxdose*ones(size(sparseoar(4).voilist,1),1);
b5=-sparseoar(4).mindose*ones(size(sparseoar(4).voilist,1),1);
b=[b1;b2;b3;b4;b5];
lb = zeros(nbeam+(2* size(sparseoar(4).voilist,1)),1);
c= [sparse(1,nbeam) ones(1, size(sparseoar(4).voilist,1)) ones(1, size(sparseoar(4).voilist,1))];
x0 = zeros(nbeam+(2* size(sparseoar(4).voilist,1)),1);
ub = [];
Aeq = [];
beq = [];
fun = @(x) myfun(x,c);
noncon = @(x) secondMomentum(x,sparseoar,A1,A2);
options = optimoptions('fmincon','Display','iter','Algorithm','active-set','SpecifyObjectiveGradient',true);%
disp("Model built for fmincon!")
%%
disp("fmincon")
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,noncon,options);

%% 
% LP model
oar1=sparseoar;
%clear sparseoar;
tic

cvx_begin
    cvx_solver gurobi_2;
    %cvx_solver_settings('timelimit', 1800)
    %cvx_solver_settings('GURO_PAR_DUMP', 1)
    %cvx_solver_settings('GURO_PAR_BARDENSETHRESH', 20)
    %cvx_solver_settings('BarIterLimit', 13)
    cvx_solver_settings('BarConvTol', 1e-4)
    %cvx_solver_settings('LogFile', 'C:\Users\Bahar\Desktop\New folder (3)\Log Files\Scenario 4.log')
    cvx_solver_settings( 'Crossover', 0 )
    cvx_solver_settings( 'Method', 2 )
    %cvx_solver_settings( 'dumpfile', 'soco2' )
    %cvx_solver_settings( 'GURO_PAR_DUMP', 1 )
    variable x(nbeam) nonnegative
    variable s(size(oar1(4).doseinfluence,1)) nonnegative
    variable su(size(oar1(4).doseinfluence,1)) nonnegative
    variable aux1(size(oar1(1).doseinfluence,1)) nonnegative
    variable aux2(size(oar1(2).doseinfluence,1)) nonnegative
    variable aux3(size(oar1(3).doseinfluence,1)) nonnegative
    variable aux4(size(oar1(4).doseinfluence,1)) nonnegative
    minimize sum(s) + sum(su)+0.00001* sum(x)
    
      
    subject to
    
        aux1 == oar1(1).doseinfluence*x;
        aux2 == oar1(2).doseinfluence*x;
        aux3 == oar1(3).doseinfluence*x;
        aux4 == oar1(4).doseinfluence*x;
        
        sum(aux1)/size(aux1,1) <= min(oar1(1).meandose, oar1(1).firstmom)
        sum(aux2)/size(aux2,1) <= min(oar1(2).meandose, oar1(2).firstmom)
        sum(aux3)/size(aux3,1) <= min(oar1(3).meandose, oar1(3).firstmom)
        %sum(aux3)/size(aux3,1) <= oar1(3).meandose
        
        aux4 - su <= oar1(4).maxdose
        aux4 + s >= oar1(4).mindose
        
        %norm(aux1,2)/size(aux1,1) <= sqrt(oar1(1).secondmom)
        %norm(aux2,2)/size(aux2,1) <= sqrt(oar1(2).secondmom)
        %sum((sparse(oar(1).doseinfluence*x)).^2)/size(oar(1).voilist,1) <= (oar(1).secondmom)
        %sum((sparse(oar(2).doseinfluence*x)).^2)/size(oar(2).voilist,1) <= (oar(2).secondmom)
cvx_end

%%
% CVaR Model
scal=1;
M = 72*scal;
OR1D=40*scal;
OR1V=0.1;
OR2D=50;
OR2V=0.01;
OR3D=30*scal;
OR3V=0.2;
OR4Max=72*scal;
OR4Min=70*scal;
eps1=45;
eps2=55;
eps3=35;
oar1=sparseoar;
%clear sparseoar;
tic

cvx_begin
    cvx_solver gurobi_2;
    %cvx_solver_settings('timelimit', 1800)
    %cvx_solver_settings('GURO_PAR_DUMP', 1)
    %cvx_solver_settings('GURO_PAR_BARDENSETHRESH', 20)
    %cvx_solver_settings('BarIterLimit', 13)
    cvx_solver_settings('BarConvTol', 1e-4)
    %cvx_solver_settings('LogFile', 'C:\Users\Bahar\Desktop\New folder (3)\Log Files\Scenario 4.log')
    cvx_solver_settings( 'Crossover', 0 )
    cvx_solver_settings( 'Method', 2 )
    %cvx_solver_settings( 'dumpfile', 'soco2' )
    %cvx_solver_settings( 'GURO_PAR_DUMP', 1 )
    variable x(nbeam) nonnegative
    variable s(size(oar1(4).doseinfluence,1)) nonnegative
    variable su(size(oar1(4).doseinfluence,1)) nonnegative
    variable aux1(size(oar1(1).doseinfluence,1)) nonnegative
    variable aux2(size(oar1(2).doseinfluence,1)) nonnegative
    variable aux3(size(oar1(3).doseinfluence,1)) nonnegative
    variable aux4(size(oar1(4).doseinfluence,1)) nonnegative
    variable d1p(size(oar1(1).doseinfluence,1)) nonnegative
    variable d2p(size(oar1(2).doseinfluence,1)) nonnegative
    variable d3p(size(oar1(3).doseinfluence,1)) nonnegative
    minimize sum(s) + sum(su)
    
      
    subject to
    
        aux1 == oar1(1).doseinfluence*x;
        aux2 == oar1(2).doseinfluence*x;
        aux3 == oar1(3).doseinfluence*x;
        aux4 == oar1(4).doseinfluence*x;
        
        aux1 - d1p <= OR1D
        aux2 - d2p <= OR2D
        aux3 - d3p <= OR3D

        sum(d1p)/(size(d1p,1)*OR1V) <= eps1-OR1D
        sum(d2p)/(size(d2p,1)*OR2V) <= eps2-OR2D
        sum(d3p)/(size(d3p,1)*OR3V) <= eps3-OR3D
        
        
        aux4 - su <= oar1(4).maxdose
        aux4 + s >= oar1(4).mindose
       
cvx_end




%% 
% SOCO model
oar1=sparseoar;
%clear sparseoar;
tic

cvx_begin
    cvx_solver gurobi;
    %cvx_solver_settings('timelimit', 1800)
    %cvx_solver_settings('GURO_PAR_DUMP', 1)
    %cvx_solver_settings('GURO_PAR_BARDENSETHRESH', 20)
    %cvx_solver_settings('BarIterLimit', 13)
    cvx_solver_settings('BarConvTol', 1e-4)
    %cvx_solver_settings('LogFile', 'C:\Users\Bahar\Desktop\New folder (3)\Log Files\Scenario 4.log')
    cvx_solver_settings( 'Crossover', 0 )
    cvx_solver_settings( 'Method', 2 )
    %cvx_solver_settings( 'dumpfile', 'soco2' )
    %cvx_solver_settings( 'GURO_PAR_DUMP', 1 )
    variable x(nbeam) nonnegative
    variable s(size(oar1(4).doseinfluence,1)) nonnegative
    variable su(size(oar1(4).doseinfluence,1)) nonnegative
    variable aux1(size(oar1(1).doseinfluence,1)) nonnegative
    variable aux2(size(oar1(2).doseinfluence,1)) nonnegative
    variable aux3(size(oar1(3).doseinfluence,1)) nonnegative
    variable aux4(size(oar1(4).doseinfluence,1)) nonnegative
    minimize sum(s) + sum(su)+0.00001* sum(x)
    
      
    subject to
    
        aux1 == oar1(1).doseinfluence*x;
        aux2 == oar1(2).doseinfluence*x;
        aux3 == oar1(3).doseinfluence*x;
        aux4 == oar1(4).doseinfluence*x;
        
        sum(aux1)/size(aux1,1) <= min(oar1(1).meandose, oar1(1).firstmom)
        sum(aux2)/size(aux2,1) <= min(oar1(2).meandose, oar1(2).firstmom)
        sum(aux3)/size(aux3,1) <= min(oar1(3).meandose, oar1(3).firstmom)
        %sum(aux3)/size(aux3,1) <= oar1(3).meandose
        
        aux4 - su <= oar1(4).maxdose
        aux4 + s >= oar1(4).mindose
        
        norm(aux1,2)/size(aux1,1) <= sqrt(oar1(1).secondmom)
        norm(aux2,2)/size(aux2,1) <= sqrt(oar1(2).secondmom)
        %sum((sparse(oar(1).doseinfluence*x)).^2)/size(oar(1).voilist,1) <= (oar(1).secondmom)
        %sum((sparse(oar(2).doseinfluence*x)).^2)/size(oar(2).voilist,1) <= (oar(2).secondmom)
cvx_end

%%
xax = linspace(1,100,100);
figure
plot(xax, dvhl(1,:),'--g',xax, dvhl(2,:),'--b',xax, dvhl(3,:),'--y',xax, dvhl(4,:),'--r');
%legend('Lungs', 'Heart', 'Esoph', 'PTV')
%legend('BrainStem', 'SpinalCord', 'Eyes', 'PTV')
hold on
plot(xax, dvhB(1,:),':g',xax, dvhB(2,:),':b',xax, dvhB(3,:),':y',xax, dvhS(4,:),':r');
plot(xax, dvhB(1,:)+dvhS(1,:),'g',xax, dvhB(2,:)+dvhS(2,:),'b',xax, dvhB(3,:)+dvhS(3,:),'y',xax, dvhB(4,:),'r');
plot(OR1D,OR1V,'^g')
plot(OR2D,OR2V,'^b')
plot(OR3D,OR3V,'^y')
hold off


