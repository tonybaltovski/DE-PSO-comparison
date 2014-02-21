%{
Copyright (c) 2012, Tony Baltovski 
All rights reserved. 

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met: 

 * Redistributions of source code must retain the above copyright notice, 
   this list of conditions and the following disclaimer. 
 * Redistributions in binary form must reproduce the above copyright 
   notice, this list of conditions and the following disclaimer in the 
   documentation and/or other materials provided with the distribution. 
 * Neither the name of  nor the names of its contributors may be used to 
   endorse or promote products derived from this software without specific 
   prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE. 

%}

%Comparison of DE and PSO using 8 benchmark functions

%Check to see if the pool is open to allow for parallel processing
if matlabpool('size') == 0 
    matlabpool open 4  %Opens 4 Matlab workers 
end
tic  %Starts a timer
K = 1000;  %Multiplier for NFC
runsToComplete = 1;   %Total number of runs need to collect stats
DEavgRes = zeros(8,2);    %Memory allocation for the results for DE
DEstdRes = zeros(8,2);    %Memory allocation for the results for DE
PSOavgRes = zeros(8,2);    %Memory allocation for the results for PSO
PSOstdRes = zeros(8,2); 
DEavgData50 = zeros(8,runsToComplete);    %Memory allocation for the results for DE
PSOavgData50 = zeros(8,runsToComplete);    %Memory allocation for the results for PSO
DEavgData100 = zeros(8,runsToComplete);    %Memory allocation for the results for DE
PSOavgData100 = zeros(8,runsToComplete);    %Memory allocation for the results for PSO


for N = 1:2
    if N == 1
        D = 50;
        %Memory allocation for the results for DE & PSO
        DE50NFC = zeros(8,K*D);
		PSO50NFC = zeros(8,K*D);
    else
        D = 100;
        %Memory allocation for the results for DE & PSO
        DE100NFC = zeros(8,K*D);
		PSO100NFC = zeros(8,K*D);
    end
    
    N_p = D;
    
   
    parfor f = 1:8
        
		%Common variables
        [H,L] = fitInt(f);  %Returns the bounds based on the benchmark function
		NFC_MAX = K*D;   %Maximum number of function calls
        
		
		%DE variables
		F = 0.5;    %Mutation Constant
		C_r = 0.9;    %Crossover rate
		X = zeros(N_p,D);   %Population memory allocation
		Xprime = zeros(N_p,D);  %Next generation population memory allocation
        V = zeros(N_p,D); 	%Noise Vector memory allocation
        U = zeros(N_p,D); 	%Trial Vector memory allocation
        DEbestFitOfRun = zeros(1,NFC_MAX);  %memory allocation for the last run to store all FCs
        DEbestFitOfRuns = zeros(1,runsToComplete);  %Temp  memory allocation
        
		%PSO variables
		c1 = 2;  %Learning Factor
		c2 = 2;  %Learning Factor
		P = zeros(N_p,D);  %Population memory allocation
        gBest = zeros(1,D); %best particle ever for whole popluation memory allocation
        pBest = zeros(N_p,D); %best particle ever for each particle memory allocation
        PSObestFitOfRuns = zeros(1,runsToComplete); %memory allocation for the last run to store all FCs
        PSObestFitOfRun = zeros(1,NFC_MAX);
        
		for loop = 1:runsToComplete
            %Initial Population to be shared by DE & PSO
            intPop = L + (H-L).*rand(N_p,D);
			
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%                     DE                                  %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            X = intPop;
            NFC = 0;
            
			while (NFC <=(NFC_MAX-D))
			% A while loop to terminate when the maximum NFC is reached
            % -D is needed since the following for loop will go pass
            % NFC_MAX since the counter is within it
                                
				for i =1:N_p
                    %A for loop to perform the crossover and mutation

                    %Generate 3 Unique Random Numbers
                    a = randi(N_p);
                    while (a == i)
                        a = randi(N_p);
                    end
                    b = randi(N_p);
                    while (a == b || b == i)
                        b = randi(N_p);
                    end
                    c = randi(N_p);
                    while  (a == c || b == c || i == c)
                        c = randi(N_p);
                    end

                    %Mutation Operator
                    V(i,:) = (X(a,:) + F .* (X(c,:) - X(b,:)));

                    %Crossover Operator 
                    for j = 1:D
                        if rand(1) < C_r
                            U(i,j) = V(i,j);
                        else
                            U(i,j) = X(i,j);
                        end

                    end  

                    %Evaluation for best fitness value of current generation
                    if fitEval(f,U(i,:),H,L,D,N_p) <= fitEval(f,X(i,:),H,L,D,N_p)
                    %Checks if the trial vector at k has a better fitness value
                    %If it does, the trial vector at k is stored in the new generation at k
                        Xprime(i,:) = U(i,:);
                    else
                        Xprime(i,:) = X(i,:);
                    end
                    NFC = NFC+1;  % Increases the NFC counter 
                
                    if NFC == 1
                    %Stores the first fitness value
                        bestFit = fitEval(f,Xprime(i,:),H,L,D,N_p);
                    end

                    
                    if bestFit > fitEval(f,Xprime(i,:),H,L,D,N_p)
                    %Checks for a better fitness value is found
                        bestFit = fitEval(f,Xprime(i,:),H,L,D,N_p);
                    end
                    
                    if loop == runsToComplete
                    %Stores the fitness values for each function call
                    %Only for the last run
                        DEbestFitOfRun(NFC) = bestFit;
                    end
                end
                   
                                      
                X = Xprime;   %Stores the new generation as the current generation
                X(X < L) = L;  %Limits new population
                X(X > H) = H;
                DEbestFitOfRuns(1,loop) = bestFit; %Stores the best fitness value found
            end
            
            %Stores the last run of each functions fitness value per FC           
            if N == 1
                    DE50NFC(f,:) = DEbestFitOfRun(1,:); 
            elseif N == 2
                    DE100NFC(f,:) = DEbestFitOfRun(1,:);
            end

           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%                           PSO                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            NFC = 0;	%Resets the NFC counter
            P = intPop;  %Take the initial population
            pBest = P;  %places values in the pBest
            v = zeros(N_p,D);  % velocity memory allocation

            while (NFC <=NFC_MAX-D)
            % A while loop to terminate when the maximum NFC is reached    
                for i=1:N_p   
                        NFC = NFC + 1;  %Increases NFC counter
                        
                        if NFC == 1
                            %Stores the first fitness value
                                gBest = P(1,:);
                        end

                        %Evaluation for best fitness value of current
                        %particle compared to its gBest record and than to
                        %gBest
                        if fitEval(f,pBest(i,:),H,L,D,N_p) > fitEval(f,P(i,:),H,L,D,N_p)
                                pBest(i,:) = P(i,:);

                        elseif fitEval(f,gBest(1,:),H,L,D,N_p) > fitEval(f,pBest(i,:),H,L,D,N_p)
                                gBest(1,:) = pBest(i,:);
                        end

                        if loop == runsToComplete
                        %Stores the fitness values for each function call
                        %Only for the last run
                            PSObestFitOfRun(NFC) = fitEval(f,gBest,H,L,D,N_p);
                        end
       
                        %Velocity calculation
                        vtemp = v(i,:)+c1.*rand(1).*(pBest(i,:)-P(i,:)) + c2.*rand(1).*(gBest-P(i,:));
                        limit = max(abs(P(i,:)));  % Determines the limit for the velocity
                        for j=1:D
                            %Limits the value in each element based on sign
                            if abs(vtemp(1,j)) > limit
                                if vtemp(1,j) < 0
                                    vtemp(1,j) = -limit;
                                elseif vtemp(1,j) > 0
                                    vtemp(1,j) = limit;
                                end
                            end
                        end 
                        v(i,:) = vtemp; %stores the calculation into the velocity matrix
                 end

                    P=P+v;  %Updates the particles
                    P(P < L) = L;  %Limits the particles
                    P(P > H) = H;
                    PSObestFitOfRuns(1,loop) = fitEval(f,gBest,H,L,D,N_p);  %Stores the best value so far

            end
            
            %Stores the last run of each functions fitness value per FC  
            if N == 1
                PSO50NFC(f,:) = PSObestFitOfRun(1,:); 
            elseif N == 2
                PSO100NFC(f,:) = PSObestFitOfRun(1,:);
            end
            
        end
        
        %Stores the best fitness value for every function's run
        if N == 1
                DEavgData50(f,:) = DEbestFitOfRuns(1,:);   
                PSOavgData50(f,:) = PSObestFitOfRuns(1,:);    
        elseif N == 2
                DEavgData100(f,:) = DEbestFitOfRuns(1,:);    
                PSOavgData100(f,:) = PSObestFitOfRuns(1,:);
        end
        
        %Calculates the mean and standard deviation
        DEavgRes(f,N) = mean(DEbestFitOfRuns); 
        DEstdRes(f,N) = std(DEbestFitOfRuns); 
        
        PSOavgRes(f,N) = mean(PSObestFitOfRuns); 
        PSOstdRes(f,N) = std(PSObestFitOfRuns);
              
           
	end
	
		     
end

toc  %returns the elapsed time in seconds
matlabpool close  %close the matlab pool
