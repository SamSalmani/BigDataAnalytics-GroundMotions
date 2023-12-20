clc;clear; close all;

% parallel computing
poolobj = gcp('nocreate');
seedchange = 'No';
nparpro=10;
if ~isempty(poolobj) && strcmp(seedchange,'Yes')
    delete(poolobj);
    parpool(nparpro);
elseif isempty(poolobj)
    parpool(nparpro);
end

% reading ground motion data
d = dir(pwd);
dirlist = d([d.isdir]);
dirlist = dirlist(~ismember({dirlist.name}, {'.','..'}));
N = size(dirlist,1);
b = 40; g = 9.81;

for i = 1: N
    stname{i,1} = dirlist(i).name;
end

for i = 1:N
    clear hvsr95 dt dataE dataN dataZ tE tN tZ aE aN aZ 
    disp(i)
    folder = [dirlist(i).folder,'/',stname{i}];
    d = dir(folder);
    folder_list = d([d.isdir]);
    folder_list = folder_list(~ismember({folder_list.name}, {'.','..'}));
      
    M = size(folder_list,1);
    for j = 1:M
        clear dataZ dataE dataN fileE fileZ fileN
        dir_to_search = [folder,'/',folder_list(j).name];   
        
        if exist ([dir_to_search,'/HHE.sac' ])
            fileE = dir([dir_to_search,'/HHE.sac']);
            dataE = rdsac([fileE.folder,'/',fileE.name]);
        elseif exist ([dir_to_search,'/HN1.sac' ])
            fileE = dir([dir_to_search,'/HN1.sac']);
            dataE = rdsac([fileE.folder,'/',fileE.name]);
        elseif exist ([dir_to_search,'/HNE.sac' ])
            fileE = dir([dir_to_search,'/*E.sac']);
            dataE = rdsac([fileE.folder,'/',fileE.name]);
        end
        dataE = rdsac([fileE.folder,'/',fileE.name]);
       
        if exist ([dir_to_search,'/HHN.sac' ]) 
            fileN = dir([dir_to_search,'/HHN.sac']);
            dataN = rdsac([fileN.folder,'/',fileN.name]);
        elseif exist ([dir_to_search,'/HN2.sac' ])
            fileN = dir([dir_to_search,'/HN2.sac']);
            dataN = rdsac([fileN.folder,'/',fileN.name]);
        elseif exist ([dir_to_search,'/HNN.sac' ])
            fileN = dir([dir_to_search,'/*N.sac']);
            dataN = rdsac([fileN.folder,'/',fileN.name]);
        end
        dataN = rdsac([fileN.folder,'/',fileN.name]);
        
        if exist ([dir_to_search,'/HNZ.sac' ]) 
            fileZ = dir([dir_to_search,'/HNZ.sac']);
            dataZ = rdsac([fileZ.folder,'/',fileZ.name]);
        elseif exist ([dir_to_search,'/HHZ.sac' ])
                fileZ = dir([dir_to_search,'/*Z.sac']);
                dataZ = rdsac([fileZ.folder,'/',fileZ.name]);
        elseif exist ([dir_to_search,'/HN3.sac' ])
                fileZ = dir([dir_to_search,'/HN3.sac']);
                dataZ = rdsac([fileZ.folder,'/',fileZ.name]);
        end
        dataZ = rdsac([fileZ.folder,'/',fileZ.name]);
        
        
        station = [dataE.HEADER.STLA,dataE.HEADER.STLO];
        
        if exist('dataE', "var")
            tE = dataE.t; tN = dataN.t; tZ = dataZ.t;
            dtE = dataE.HEADER.DELTA; dtN = dataN.HEADER.DELTA; dtZ = dataZ.HEADER.DELTA;
            dtlist(j,:) = [dtE dtN dtZ];
            event(j,:) = [dataE.HEADER.EVLA,dataE.HEADER.EVLO,dataE.HEADER.MAG,dataE.HEADER.DIST, ];
        end 
    end
    dt = min(min(dtlist));

    count = 0;
    for j = 1:M
        
        dir_to_search = [folder,'/',folder_list(j).name];   
        if exist ([dir_to_search,'/HHE.sac' ])
            fileE = dir([dir_to_search,'/HHE.sac']);
            dataE = rdsac([fileE.folder,'/',fileE.name]);
        elseif exist ([dir_to_search,'/HN1.sac' ])
            fileE = dir([dir_to_search,'/HN1.sac']);
            dataE = rdsac([fileE.folder,'/',fileE.name]);
        elseif exist ([dir_to_search,'/HNE.sac' ])
            fileE = dir([dir_to_search,'/*E.sac']);
            dataE = rdsac([fileE.folder,'/',fileE.name]);
        end
        %dataE = rdsac([fileE.folder,'/',fileE.name]);
       
        if exist ([dir_to_search,'/HHN.sac' ])
            fileN = dir([dir_to_search,'/HHN.sac']);
            dataN = rdsac([fileN.folder,'/',fileN.name]);
        elseif exist ([dir_to_search,'/HN2.sac' ])
            fileN = dir([dir_to_search,'/HN2.sac']);
            dataN = rdsac([fileN.folder,'/',fileN.name]);
        elseif exist ([dir_to_search,'/HNN.sac' ])
            fileN = dir([dir_to_search,'/*N.sac']);
            dataN = rdsac([fileN.folder,'/',fileN.name]);
        end
        dataN = rdsac([fileN.folder,'/',fileN.name]);
        
        if exist ([dir_to_search,'/HNZ.sac' ])
            fileZ = dir([dir_to_search,'/*Z.sac']);
            dataZ = rdsac([fileZ.folder,'/',fileZ.name]);
        elseif exist ([dir_to_search,'/HHZ.sac' ])
            fileZ = dir([dir_to_search,'/*Z.sac']);
            dataZ = rdsac([fileZ.folder,'/',fileZ.name]);
        elseif exist ([dir_to_search,'/HN3.sac' ])
            fileZ = dir([dir_to_search,'/HN3.sac']);
            dataZ = rdsac([fileZ.folder,'/',fileZ.name]);
        end
        dataZ = rdsac([fileZ.folder,'/',fileZ.name]);
        
        
        if exist('dataE', "var")
            tE = dataE.t; tN = dataN.t; tZ = dataZ.t;
            dtE = dataE.HEADER.DELTA; dtN = dataN.HEADER.DELTA; dtZ = dataZ.HEADER.DELTA;
            nE =length(tE); nN =length(tN); nZ =length(tZ);
        end 
        
        
        if dtE~=dt
            aE = resample(detrend(dataE.d),dtE/dt,1);
            tE = linspace(tE(1),tE(end),size(aE,1)); tE=tE';
        else
            aE = detrend(dataE.d);
        end
        if dtN~=dt
            aN = resample(detrend(dataN.d),dtN/dt,1);
            tN = linspace(tN(1),tN(end),size(aN,1)); tN=tN';
        else
            aN = detrend(dataN.d);
        end
        if dtZ~=dt
            aZ = resample(detrend(dataZ.d),dtZ/dt,1);
            tZ = linspace(tZ(1),tZ(end),size(aZ,1)); tZ=tZ';
        else
            aZ = detrend(dataZ.d);
        end
        
        IE = cumtrapz(tE,aE.^2);IE = IE/IE(end);
        IN = cumtrapz(tN,aN.^2);IN = IN/IN(end);
        IZ = cumtrapz(tZ,aZ.^2);IZ = IZ/IZ(end);
        [IE,iii] = unique(IE); tE = tE(iii);
        [IN,iii] = unique(IN); tN = tN(iii);
        [IZ,iii] = unique(IZ); tZ = tZ(iii);
                
        hhh = figure(1);clf;hold all; plot(tE,IE);plot(tN,IN);plot(tZ,IZ);drawnow; 
        lowerboundE = find(IE >= 0.05,1);
        upperboundE = find(IE <= 0.95); upperboundE=upperboundE(end);
        lowerboundN = find(IN >= 0.05,1);
        upperboundN = find(IN <= 0.95); upperboundN=upperboundN(end);
        lowerboundZ = find(IZ >= 0.05,1);
        upperboundZ = find(IZ <= 0.95); upperboundZ=upperboundZ(end);

        lowerbound = max([lowerboundE, lowerboundN, lowerboundZ]);
        upperbound = min([upperboundE, upperboundN,upperboundZ]);

       
        PGA(j,:) = [max(abs(aN)) max(abs(aE)) max(abs(aZ))];
        event(j,:) = [dataE.HEADER.EVLA,dataE.HEADER.EVLO,dataE.HEADER.MAG,dataE.HEADER.DIST];
       
        if PGA(j,1)>0.00005*g && PGA(j,2)>=0.00005*g %&& lenm95-len095+1>0
            
            L = tukeywin( upperbound - lowerbound +1, 0.05);
            aE95 = L.*aE( lowerbound : upperbound );
            aN95 = L.*aN( lowerbound : upperbound );
            aZ95 = L.*aZ( lowerbound : upperbound );
        
            
            fs = 1/dt; nfft = 2^15;
            
            nfreq = round(100/(fs/nfft));
            
            freq = linspace(0,1,nfft)*fs;
            
            aN95fft = fft(aN95,nfft);
            aE95fft = fft(aE95,nfft);
            aZ95fft = fft(aZ95,nfft);
            
            h195 = (aN95fft).^2;%/trapz(freq,abs(aN95fft).^2);
            h295 = (aE95fft).^2;%/trapz(freq,abs(aE95fft).^2);
            v95  = (aZ95fft);%/trapz(freq,abs(aZ95fft).^2);
            count = count+1; 

            hvsr95(:,count)  = 1/sqrt(2)*fasterKonnoOhmachi(abs(sqrt((h195(1:nfreq)+h295(1:nfreq)))./(v95(1:nfreq)+1e-6)),freq(1:nfreq),25);
            hvsr95(:,count)  = 1/sqrt(2)*(abs(sqrt((h195(1:nfreq)+h295(1:nfreq)))./(v95(1:nfreq)+1e-6)));
            
        end
    end
    hvsr_mean295 = exp(mean(log(hvsr95),2));
    hvsr_mean295 = mean((hvsr95),2);
    hvsr_std295 = std((hvsr95),1,2);


% PLOTTING
     mypath0 = [pwd,'/'];
     dir='Est'; mypath=[mypath0,dir]; mkdir(mypath);
     mkdir(fullfile(mypath0,dirlist(i).name));
     mypath2=[pwd,'/',dirlist(i).name];

     h = figure(2);figinit(8,6,h);clf;
     scatter(event(:,4),event(:,3))
     ylabel('Magnitude'); xlabel('Distance'); ylim([3 5]); xlim([0 200]);
     print(h,'-dpdf',[dataE.HEADER.KNETWK,'_',dataE.HEADER.KSTNM,'_Mag_Dist.pdf']);

    csvwrite([dataE.HEADER.KNETWK,'_',dataE.HEADER.KSTNM,'_HVSR95_2power15_with_knno75_1.txt'],[freq(1:nfreq)' hvsr_mean295(1:nfreq,1) hvsr_std295(1:nfreq,1) hvsr95]);
    close all;
end
