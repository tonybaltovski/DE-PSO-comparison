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

function [fitVal] = fitEval(f,A,H,L,D,N_p)
%Performs fitness evaluations based on the fitness fucntion [f], the vector
%submitted [A], the bounds [H,L], dimensionality [D] and the population
%size [N_p] and returns the value [fitVal]. **Note all these parameters are
%needed for the following functions.

A(A < L) = L;  %Limits the input values based on the bounds
A(A > H) = H;
fitVal = 0;  %variable initalization

    switch(f)   
        case 1 % 1st De Jong
            fitVal = sum(A.^2);  
        case 2 % Axis Parallel Hyer-Ellipsoid
            for j=1:D
                fitVal = fitVal + j*A(1,j)^2;
            end
        case 3 % Schwefel's Problem 1.2
            fitVal = sum((cumsum(A)).^2);
        case 4 % Rosenbrock's Vally
            for j=1:(D-1)
                fitVal = (100*(A(1,(j+1))-A(1,j)^2)^2+(1-A(1,j))^2) + fitVal;
            end
        case 5 % Rastrigin's Function
            fitVal = 10*D + sum((A.^2)-10.*cos((2*pi).*A));
        case 6 % Griewangk's Function
            product = 1;
            for j=1:D
                k = (A(1,j)/(sqrt(j)));
                product = cos(k) * product;
            end
            fitVal = 1 + ((1/4000)*sum(A.^2))-product;
        case 7 %Sum of Different Power
            for j=1:D
                fitVal = abs(A(1,j))^(j+1) + fitVal;
            end
        case 8 % Ackley's Problom
            sq = A(1,:).^2;
            p2 = 2.*pi.*A(1,:);
            c = sum(cos(p2));
            s = sum(sq);
            fitVal = 20 + exp(1) - (20*(exp(-0.2*(sqrt((s/D)))))) - (exp((c/D)));
    end
end

        


