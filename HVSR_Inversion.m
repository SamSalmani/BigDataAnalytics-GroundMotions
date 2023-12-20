clc; clear; close all;

% Creating parallel computing toolbox to use 10 cores
poolobj = gcp('nocreate');
seedchange = 'No';
nparpro=10;
if ~isempty(poolobj) && strcmp(seedchange,'Yes')
    delete(poolobj);
    parpool(nparpro);
elseif isempty(poolobj)
    parpool(nparpro);
end

%% Getting directory of HVSR data derived using HVSR code (freq, hvsr-mean, hvsr-std)
d = dir(pwd+'*.txt'); % reading all hvsr files with txt format created by HVSR code using ground motion data
mkdir result      % creating a directory to save the results for each station

%% Create optimization options
options = optimoptions('quadprog','Display','off');

%% Start Performing HVSR Inversion
for ss = 1:size(d,1) 
    clearvars -except ss d options
    
    data_hvsr = readtable([pwd,d(ss).name]);  %reading HVSR file
    freqhvsrreal = data_hvsr.Var1;     % range of frequncy in experimental hvsr
    yhvsrreal0 = data_hvsr.Var2;       % hvsr values in experimental hvsr
    freq = logspace(log10(0.4),log10(50),500);         % defining a frequency row vector of 500 logarithmically spaced points 
    yhvsrreal = interp1(freqhvsrreal,yhvsrreal0,freq,'linear','extrap'); % linear extrapolation over eexperimental hvsr values over "freq" as query points
    yhvsrreal = yhvsrreal';
        
    PTable = [0:1:5 7.5:2.5:30 35:5:50 60:10:100 125:25:350 logspace(log10(375),log10(5000),30)]; % Defining Soil layers 
    PTable = PTable'; soilThick = max(PTable);
    
    yreal  = yhvsrreal; % assigning real hvsr values to yreal 
    Rfac   = 0.05*yhvsrreal; %)*ones(size(yhvsrreal));  %% noise created based on real values to ba added to predctions
    
    Gamma  = diag(Rfac).^2;   % diagonal matrix of noise values 
    
    % deriving thickness of each soil layer
    for i = 1:length(PTable)-1
        hTable(i,1) = PTable(i+1)-PTable(i);
    end

    hTable = [hTable;0]; % adding 0 as the thickness of halfspace
    nb = length(hTable); % number of layers
    nparam = 2*nb; % number of parameters (Vs and Vp for each layer)
    nParticle = 3*nparam; % creating intial ensemble of N particles 

    % defining A as coefficient matrix that should satisfy Au<g (g includes defined constraints)
    % where uT = [VS1, VSr, VP1, . . . , VPr] is a vector of size 2nb
    % A is a coefficient matrix of size 3nb×2nb and g is an upperbound vector of size 3nb
    Auu1 = zeros(nb-1,nparam); % 64*130
    Auu2 = zeros(nb  ,nparam); % 65*130
    Auu3 = zeros(nb-1,nparam); % 64*130

    % updating A matrices based on constraints
    % Auu1 for constructing Vs(i) <= Vs(i+1)
    for i = 1:nb-1
        Auu1(i,i)   = 1;
        Auu1(i,i+1) =-1;
    end

    % Auu2 for constructing Vs(i) <= 1.5*Vp(i)
    for i = 1:nb
        Auu2(i,i)    = 1.5;
        Auu2(i,nb+i) = -1;
    end

    % Auu3 for constructing Vp(i) <= Vp(i+1)
    for i = 1:nb-1
        Auu3(i,nb+i)    = 1;
        Auu3(i,nb+i+1) = -1;
    end

    % Adding 3 rows to be updated afterward to add following constraints: 
    Auu = [Auu1;Auu2;Auu3;zeros(3,nparam)];
    
    ind = size(Auu,1); 
    Auu(ind-2,1 ) =-1; % constraint: -Vs1<-50
    Auu(ind-1,nb) = 1; % constraint: Vs(end)<5000
    Auu(ind-0,2*nb) = 1; % constraint: Vp(end)<10000
    gu = [zeros((3*nb)-2,1);-50/200;5000/200;10000/200]; % upperbound values vector
    
    % Defining Vs values for each depth (row=i)
    for i =1:length(PTable)
        ParticleVs(i,:) = 2.5+12.5*(PTable(i)/soilThick)^0.5*rand(1,nParticle);
    end
    ParticleVs(1,:) = 2.5+12.5*(PTable(2)/soilThick)^0.5*rand(1,nParticle);  % updating Vs value of the first row based on the bottom side of the first layer
    
    uhat = [ParticleVs;2*ParticleVs]; %joining Vs and Vp values(considering Vp=2*Vs)-uhat=[Vs1,Vs2,...,Vsr,Vp1, Vp2,..Vpr]
    
    count = 0;
    for j = 1:nParticle 
        theta = uhat(:,j); %reading the values of Vs and Vp for each particle(130*1)
        V = find(Auu*theta>gu); %finding the index of rows which violates the rule : Au<g
        if isempty(V) %if empty,nothing happens
        else % use quadprog optimization to update uhat to satisfy the constraint
            count = count+1;
            M    =  eye(nparam);
            f    =  -theta;
            uhat(:,j) = quadprog(M,f,Auu,gu,[],[],[],[],[],options);
        end
    end

    uvec(:,1) = mean(uhat,2); %calculating the mean of all columns of uhat
   
    iteration = 0; maxiter = 300;
    uhattemp = uhat;
    
    s = 1; 
    Gammau = diag((0.05*uvec).^2); %considering a diagonal noise matrix on uvec
    
    mypath0 = [pwd,'/'];
    dir='Est'; mypath=[mypath0,dir]; mkdir(mypath);
    
    fig1 = figure(1); fig2 = figure(2); fig3 = figure(3); %fig4 = figure(4);
    
    while iteration <  maxiter
        iteration = iteration+1;
        %% run in parallel
        for ii = 1:nParticle
            mypath=[mypath0,'Est/S',num2str(ii-1)];
            writematrix(uhat(:,ii),[mypath,'param.txt']); %creating txt file for each column of uhat(Vs,Vp values)
        end
        parfor ii = 1:nParticle
            what(:,ii) = hvsr([mypath0,'Est/S',num2str(ii-1)],hTable,freq); %calculating theoritical hvsr values for each set of Vs and Vps (uhat)
        end
        
        for ii = 1:nParticle
            mypath=[mypath0,'Est/S',num2str(ii-1)];
            what(:,ii) = csvread([mypath,'hvsr.out']); %reading created hvsr files and joining them
        end
        uhatold = uhat;
        
        wnp1 = mean(what,2); %calculating mean of hvsr of all particles 
        unp1 = mean(uhat,2); %calculating mean of uhat of all particles 
        Cuwnp1 = zeros(nparam,length(yreal)); %covariance matrix of uhat and what
        Cwwnp1 = zeros(length(yreal)); %covariance matrix of what
        for i = 1:nParticle
            Cuwnp1   = Cuwnp1   + (uhat(:,i)-unp1)*(what(:,i)-wnp1)'/(nParticle-1);
            Cwwnp1   = Cwwnp1   + (what(:,i)-wnp1)*(what(:,i)-wnp1)'/(nParticle-1);
        end
        
        Ik = eye(size(Cwwnp1)); %Identity matrix-size 
        K = Cuwnp1*((Cwwnp1+Gamma)\Ik); % Adding noise to Cww
        
        temp = mvnrnd(zeros(length(yreal),1),Gamma,nParticle); % returns a matrix "temp" of "nParticle"
        % random vectors chosen from the same multivariate normal distribution, 
        % with mean vector  "zeros" and covariance matrix Gamma
        
        yperturb = temp'; 
        
        for i = 1:nParticle
            dy = yreal-what(:,i)+s*yperturb(:,i); %calculating the difference between experimental and theoritical
            % hvsr and adding noise to them
            uhattemp(:,i) = uhat(:,i) + K*dy; % ???? which part of the paper?
        end
        
        CC = zeros(size(Auu,1) ,nParticle); 
        Bu = zeros(nparam      ,nParticle); 
        Bw = zeros(length(yreal),nParticle); 
        for i = 1:nParticle
            Bu(:,i) = uhat(:,i) - unp1;
            Bw(:,i) = what(:,i) - wnp1;
        end
        for i = 1:nParticle
            V = find(Auu*uhattemp(:,i)>gu);
            if isempty(V)
                uhat(:,i) = uhattemp(:,i);
            else
                M    =  1/nParticle/nParticle*Bw'*(Gamma\Bw)+eye(nParticle)/nParticle;
                f    =  1/nParticle*Bw'*(Gamma\(what(:,i)-yreal-s*yperturb(:,i)));
                A    = Auu*Bu/nParticle;
                ub   = gu-Auu*uhat(:,i);
                bhat = quadprog(M,f,A,ub,[],[],[],[],[],options);
                uhat(:,i) = uhat(:,i)+1/nParticle*Bu*bhat;
            end
        end
        

        
        temp = mvnrnd(zeros(nparam,1),Gammau,nParticle);
        
        uperturb = temp';
        
        if iteration < maxiter 
            uhat = uhat+uperturb;
            for j = 1:nParticle
                while ~isempty(find(Auu*uhat(:,j)>gu,1))
                    M     =  eye(nparam);
                    f     =  -uhat(:,j);
                    uhat(:,j) = quadprog(M,f,Auu,gu,[],[],[],[],[],options);
                end
            end
        end
        unp1   = mean(uhat,2); 
        uvec(:,iteration+1) = unp1;
        
        writematrix(unp1,[mypath0,'param.txt'])
        yhvsr = hvsr(mypath0,hTable,freq);
        yhvsr = csvread([mypath0,'hvsr.out']);
        set(0,'CurrentFigure',fig1);clf; hold all;
        plot(freq,yhvsrreal','r','linewidth',1.5);
        plot(freq,yhvsr,'b-','linewidth',1.5);
        set(gca,'xscale','log');
        xlim([0.1 100]);
        grid on;
        drawnow;
        
        set(0,'CurrentFigure',fig2);clf;
        hold all;
        for i=1:nParticle
            stairs([uhat(1,i);uhat(1:nb,i)]*200,-[PTable;soilThick+10],'color',[0.5 0.5 0.5]);
            stairs([uhat(nb+1,i);uhat(nb+1:2*nb,i)]*200,-[PTable;soilThick+10],'color',[0.5 0.5 0.5]);
        end
        stairs([unp1(1);unp1(1:nb)]*200,-[PTable;soilThick+10],'r','linewidth',1.5);
        stairs([unp1(nb+1);unp1(nb+1:2*nb)]*200,-[PTable;soilThick+10],'r','linewidth',1.5);
        drawnow;
        
        set(0,'CurrentFigure',fig3);clf; hold all;
        for i =1:nparam
            plot(1:iteration+1,uvec(i,:),'r-x');grid on;
        end
        drawnow
    end
     plot(freq,yhvsrreal','r','linewidth',1.5); hold all;
     plot(freq,yhvsr,'b-','linewidth',1.5);
     set(gca,'xscale','log'); xlim([0.1 100]);
     grid on;xlabel('Freq'); ylabel("HVSR");title(strcat('NP-',d(ss).name(4:end-10)));drawnow; 
     fig_name = strcat('FIG_',d(ss).name(1:end));
     saveas(fig1, [[pwd+ string(fig_name)+".png"]);
     mkdir([[pwd,d(ss).name]);
     saveas(fig1, [[pwd,d(ss).name,'/hvsr.png']);
     saveas(fig3, [[pwd,d(ss).name,'/iteration.png']);
    csvwrite([pwd,d(ss).name,'/uvec.txt'],uvec);
    csvwrite([pwd,d(ss).name,'/unp1.txt'],unp1*200);
    csvwrite([[pwd,d(ss).name,'/uhat.txt'],uhat*200);
    csvwrite([[pwd,d(ss).name,'/depth.txt'],PTable);
end
