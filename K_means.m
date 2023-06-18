clc, clear;
%Vstupy
m = readtable('data.csv'); %Soubor vstupnich dat
metricsType = 1;  %Volba metriky 1-čtverec euklidovske, 2-euklidovska, 3-manhattanska
k = 6; %Pocet clusteru
centroids_manual = false; %Manualni zadavani centroidu, pokud false, tak automaticky
steps = true; %Krokovani vypoctu
c_manual = [5.6312 5.6375; %Manualni zadavani souradnic centroidu
            5.0163 3.4408;
            7.5062 3.0125];
%-------------------------------------------------------------------------------------
m = table2array(m);

%Vytvoreni matice centroidu
if centroids_manual %Manualni volba
    c = c_manual;
else %Automatick8 volba
    max_limit = max(m);
    min_limit = min(m);
    c = zeros(k,length(m(1,:)));
    for i = 1:k
        for k = 1:length(m(1,:))
            c(i,k) = (max_limit(k) - min_limit(k)).*rand(1,1) +  min_limit(k);
        end
    end
end

D=zeros(length(m(:,1)), length(c(1,:))); %Matice D, matice vzalenosti prvku od centroidu
G=zeros(length(m(:,1)), length(c(1,:))); %Matice G, Group matrix, matice prislusnosti prvku do clusteru
Gold=zeros(length(m(:,1)), length(c(1,:))); %%Matice G z minule iterace
Gnew=zeros(length(m(:,1)), length(c(1,:)))+1; %%Matice G soucasna
numIterations = 0; %Pocet iteraci

%Zobrazeni neprirazenych prvku mnoziny a nahodne vygenerovanych centroidu
if steps
    for i = 1:k
        scatter(m(:,1),m(:,2),'filled','k');
        hold on
    end
    scatter(c(:,1),c(:,2),200,'x', 'k');
    title 'Cluster Assignments and Centroids'
    disp"Stisknete enter pro dalsi krok!";
    xlabel('x1'), ylabel('x2');
    pause;
end

AlreadyInsomegroup = false;
while ~isequal(Gold,Gnew); %cyklus se bude opakovat, dokud se prvnky neustali v urcitem clusteru
    Gold=G;
    for i=1:length(c(:,1))
        for j=1:length(m(:,1))
            
            % volba metriky
            switch metricsType
                case 1 %Ctverec euklidovske metriky
                    D(j,i)= ((m(j,1)-c(i,1))^2+(m(j,2)-c(i,2))^2);
                case 2 %Euklidovska metrika
                    D(j,i)= sqrt((m(j,1)-c(i,1))^2+(m(j,2)-c(i,2))^2);
                case 3 %Manhattanská metrika
                    D(j,i)= (m(j,1)-c(i,1))+(m(j,2)-c(i,2));
            end
        end
    end
    
    %Generovani matice G
    for i= 1:length(D(:,1))
        minimum = min(D(i,:));
        for k= 1:length(D(1,:))
            if D(i,k) == minimum & ~AlreadyInsomegroup
                G(i,k)=1;
                AlreadyInsomegroup = true;
            else
                G(i,k)=0;
            end
        end
        AlreadyInsomegroup = false;
    end
    
    %Vypocet centroidu
    ctemp = c.* 0;
    divtemp = 0;
    for  i= 1 : length(G(1,:))
        % i je kazdy jeden cluster, takze sloupec
        for j =1 : length(G(:,1))
            % j je kazdy jeden bod tedy radek G
            ctemp(i,1) = ctemp(i,1) + G(j,i)*m(j,1);
            ctemp(i,2) = ctemp(i,2) + G(j,i)*m(j,2);
            divtemp = divtemp + G(j,i);
        end
        if divtemp ~= 0
            ctemp(i,1) = ctemp(i,1)./divtemp;
            ctemp(i,2) = ctemp(i,2)./divtemp;
        else
            ctemp(i,1) = ctemp(i,1);
            ctemp(i,2) = ctemp(i,2); 
        end
        divtemp = 0;
    end
    c = ctemp;
    ctemp = c .* 0;
    Gnew = G;
    
    numIterations = numIterations +1;
    
    %Krokovani
    if steps
        clusterPlot(m,G,c);
        disp"Stisknete enter pro dalsi krok!";
        pause;
        clc;
    end
end
close all;
G=Gnew;
clusterPlot(m,G,c);
writematrix(c,"Centroids.txt");
