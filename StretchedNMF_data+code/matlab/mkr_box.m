function [x, s, t, A1, A2, objnew, k, maxfloor, solves, avg_inact, I1,I2, floor, time1, tol, sil, flag] = mkr_box( Q, d, a, b, cA, A1, A2, I1, I2, floor,tol,avg_inact,solves,maxfloor,sil)


% straightforward extension from the one-sided case
% notation has to be adapted to be consistent


% call: [x,alpha,Aopt,time,iter,objvalue] = hr_performance(Q,q,a,b)
% default starting active set is empty
% enter data vectors q and b as column vectors
% solves: min J(x) = q'x + (1/2) x'Qx  subject to: x  <= b
% Q is assumed to be symmetric positive definite
% KKT:    Qx + s + t + q = 0; a <= x <= b, s >=0, t<= 0, s'(x-b) = 0, t'(a-x) = 0;
% input:    (Q,q,a,b) problem data,
%           cA=0 ... initial active set empty set
%           cA=1 ... initial acitve set full set
%           cA=2 ... initial active set defined by A
%           A ... guess on initial active set, e.g. A=(1:n)???; 
% output:   (x,s,t) optimal solution,
%           (A1,A2) = active sets at opt. sol. 
%           I2 = inactive set at the opt. sol.
%           time1 ... time in seconds
%           k ... # iterations
%           objnew ... optimal objective value
time1 = 0;
flag = 1;
n = length(d);                  % problem size
done = 0;                        % not yet done
k = 0;
nosub = 0;
innerittotal = 0;
if(nargin == 4)
    cA = 0; A1 = []; A2 = []; I1 = []; I2 = []; floor= 0; tol = 10^-10; sil = 1;
end

if (cA == 1 && floor == 0) % set all variables on the upper bound
    tic;
    s = -d - Q*b;
    A1 = find(s>0);
    I1 = ones(n,1);
    I1(A1) = 0;
    I1 = find(I1>0);
    I2 = [];
    k = 0;
    solves = 1;
else if (cA == 0 & floor == 0)
        tic;
        A1 = [];
        A2 = [];
        I1 = 1:n;
        I2 = [];
        solves = 0;
        else if(cA == 2 & floor == 0) % initial active set is given by the input data
            tic;
            I1 = ones(n,1);
            I1(A1) = 0;
            I1(A2) = 0;
            I1 = find(I1>0);
            I2 = [];
            solves = 0;
            end
    end
end


I = [I1;I2];
N = [A1;A2;I];
N = sort(N);
objnew = 10^20;

% main loop 
while done < 1;                  % while not optimal
     k = k + 1;                  % start a new iteration
     innerdone = 0;
     Aold1 = A1;
     Aold2 = A2;
     A11sum = [];
     A12sum = [];
     A21sum = [];
     A22sum = [];
     innerit = 0;
     objold = objnew;
% solve system KKT(A):
    while innerdone < 1
     innerit = innerit + 1;
     if innerit > 1
         help = zeros(n,1);
         help(I1) = 1;
         help(I1check) = 0;
         I1 = find(help>0);
         help = zeros(n,1);
         help(I2) = 1;
         help(I2check) = 0;
         I2 = find(help>0);
     end
     I = [I1;I2];
     I = sort(I);
     x = zeros(n,1);                      % x( A) = b( A); 
     x(A1) = b(A1);
     x(A2) = a(A2);
     avg_inact = avg_inact + length(I);
     solves = solves + 1;
     if isempty(I) == 0;          
         rhs = -d(I);
         if isempty(A1) == 0;       % update right hand side rhs
           rhs = rhs - Q(I,A1)*b(A1);   
         end; 	   
         if isempty(A2) == 0;       % update right hand side rhs
           rhs = rhs - Q(I,A2)*a(A2);   
         end; 	
         x(I) = Q(I,I) \ rhs;
     end;
     innerdone = max(x-b) <= 0 & max(a-x) <= 0;  % if there are variables that are primal infeasible then add them to the according active set and reoptimize
     if innerit >= 100; flag = 0; s = []; t = []; floor = floor - 1;% emergency exit to avoid cycling
         if(sil==0)
            fprintf('max number of inner iterations reached. \n');
         end
         return
     end;
     if (innerdone < 1)
         help1 = x>=b;
         help2 = zeros(n,1);
         help2(I1) = 1;
         I1check1 = help1.*help2;
         I1check1 = find(I1check1>0);
         A11sum = [A11sum;I1check1];
         help2 = zeros(n,1);
         help2(I2) = 1;
         I2check1 = help1.*help2;
         I2check1 = find(I2check1>0);
         A12sum = [A12sum;I2check1];
         Anew = [I1check1;I2check1];
         A1 = [A1;Anew];
         help1 = x<=a;
         help2 = zeros(n,1);
         help2(I1) = 1;
         I1check2 = help1.*help2;
         I1check2 = find(I1check2>0);
         A21sum = [A21sum;I1check2];
         help2 = zeros(n,1);
         help2(I2) = 1;
         I2check2 = help1.*help2;
         I2check2 = find(I2check2>0);
         A22sum = [A22sum;I2check2];
         Anew = [I1check2;I2check2];
         A2 = [A2;Anew];
         I1check = [I1check1;I1check2];
         I2check = [I2check1;I2check2];
     end
    end
     I2 = [I2;I1];
     s = zeros( n,1);  
     if isempty(A1) == 0;            % backsubsitute for s(A), if |A|>0 
         s(A1) = -d(A1) - Q(A1,N)* x(N);
     end;
     
     t = zeros(n,1);
     if isempty(A2) == 0;            % backsubsitute for s(A), if |A|>0 
         t(A2) = -d(A2) - Q(A2,N)* x(N);
     end;
     


    objnew = 0.5*x(N)'*Q(N,N)*x(N)+d(N)'*x(N);   % compute objective value of the new point
    
     if(nosub == 0) % if number of constraints is reduced then "nosub = 1" as we start with a "new" problem
     if (objnew >= objold)  % solve a smaller subproblem if the objective value did not improve
         A1sum = [A11sum;A12sum];
         A2sum = [A21sum;A22sum];
         help11 = s < 0;
         help21 = zeros(n,1);
         help21(A1sum) = 1;
         help12 = t > 0;
         help22 = zeros(n,1);
         help22(A2sum) = 1;
          if(isempty(Aold1) & isempty(Aold2)) % if Aold1 and Aold2 (in paper B_s and D_s) are empty take one element from A1sum (in paper B_1 \cup B_2) 
             % or if A1sum is also empty then A2sum (D_1 \cup D_2 in the
             % paper) and if A2sum is also empty then from I_1 (in the
             % paper J_1 \cup J_2) which is the last part of I2.
             if(isempty(A1sum))
                 if(isempty(A2sum))
                    Add = I2(end);
                    I2 = I2(1:end-1);
                    Aleft1 = Add;
                    Aleft2 = [];
                 else
                    Add = A2sum(1);
                    help22(A2sum(1)) = 0;
                    Aleft2 = Add;
                    Aleft1 = [];
                 end
             else
                Add = A1sum(1);
                help21(A1sum(1)) = 0;
                Aleft1 = Add;
                Aleft2 = [];
             end
          else
              Aleft1 = Aold1;
              Aleft2 = Aold2;
          end
         I1inner1 = help11.*help21;
         I1inner1 = find(I1inner1>0);
         I1inner2 = help12.*help22;
         I1inner2 = find(I1inner2>0);
         I1inner = [I1inner1;I1inner2];
         help11 =  s>= 0;
         Ainner1 = help11.*help21;
         Ainner1 = find(Ainner1>0);
         help12 =  t<= 0;
         Ainner2 = help12.*help22;
         Ainner2 = find(Ainner2>0);
         floor = floor + 1;
         if(isempty(Aleft1))
             d_new = d+Q(:,Aleft2)*a(Aleft2);
         else if(isempty(Aleft2))
                 d_new = d+Q(:,Aleft1)*b(Aleft1);
             else
                 d_new = d+Q(:,Aleft1)*b(Aleft1)+Q(:,Aleft2)*a(Aleft2);
             end
         end                     
             [x, s, t, A12sum, A22sum, objnew, dummy, maxfloor, solves, avg_inact, I1, I2, floor, dummy, tol, sil, flag] = mkr_box( Q, d_new, a, b,cA, Ainner1, Ainner2, I1inner,I2, floor, tol, avg_inact, solves, maxfloor, sil); % solve a subproblem of smaller dimension
             if(flag == 0)
                 return;
             end
             floor = floor - 1;
             I = I2;  % I1 is empty by construction
             A1 = [Aleft1;A12sum];
             A2 = [Aleft2;A22sum];
             x = zeros(n,1);   % recompute x und s because d is different in the subproblem
             x(A1) = b(A1);
             x(A2) = a(A2);
             if isempty(I) == 0;          
                 ltri = chol((Q( I,I)));  % cholesky factor
                 rhs = -d(I);
                 if isempty(A1) == 0;       % update right hand side rhs
                   rhs = rhs - Q(I,A1)*b(A1);   
                 end; 	   
                 if isempty(A2) == 0;       % update right hand side rhs
                   rhs = rhs - Q(I,A2)*a(A2);   
                 end; 	
                 xtmp = (ltri') \ rhs; 
                 x( I) = ltri \  xtmp;   % solve for inactive variables
             end;
              s = zeros(n,1);
             if isempty(A1) == 0            % backsubsitute for s(A), if |A|>0 
                s(A1) = -d(A1) - Q(A1,N)*x(N);  
             end;
             t = zeros(n,1);
             if isempty(A2) == 0            % backsubsitute for s(A), if |A|>0 
                t(A2) = -d(A2) - Q(A2,N)*x(N);  
             end;
             objnew = 0.5*x(N)'*Q(N,N)*x(N)+d(N)'*x(N); % update objective value for new primal point
     end
     end
     
     Apre1 = A1;    % compute the new active sets
     Apre2 = A2;
     I11 = find(s<0);
     I12 = find(t>0);
     I1 = [I11;I12];
     I2h = find(x>=b | x<=a);
     I2 = zeros(n,1);
     I2(N) = 1;
     I2(I2h) = 0;
     I2 = find(I2>0);
     I = [I1;I2];
     I = sort(I);
     A1 = zeros(n,1);
     A1(N) = 1;
     A1(I) = 0;
     A1(Apre2) = 0;
     A1 = find(A1>0);
     A2 = zeros(n,1);
     A2(N) = 1;
     A2(I) = 0;
     A2(Apre1) = 0;
     A2 = find(A2>0);
     innerittotal = innerittotal + innerit;
     
     nosub = 0;
     
     if(isempty(A1) == 1 & isempty(A2) == 1 & length(Apre1) + length(Apre2) == 1) % in this special case, one bound can be removed
        if(isempty(Apre1)==0)
            b(Apre1) = 10^20;
        else
            a(Apre2) = -10^20;
        end
        nosub = 1; % "new problem" - in first iteration no subproblem allowed
     end
   if (sil == 0)
     err = log10(max( norm(Q(N,N)*x(N)+d(N)+s(N)+t(N),1), 1e-99));
     err_p1 = log10(max( max(x-b),1e-99));
     err_p2 = log10(max( max(a-x),1e-99));
     err_d1 = log10(max( -min(s), 1e-99));
     err_d2 = log10(max( max(t), 1e-99));
     fprintf(' %3.0d %3.0d %6.0f  %6.0f %6.0f  %6.0f  %10.2f %10.2f %10.2f %10.2f %10.2f  %10.10f \n',[k innerit length(I)   length(A1) length(A2) floor   err_p1  err_p2   err_d1  err_d2  err  objnew]);
   end
   I11tol = find(s<-tol);
   I12tol = find(t>tol);
   I1tol = [I11tol;I12tol];  
   done = isempty(I1tol); % if I1 is empty the new point is also dual feasible and hence optimal -
                             % we allow for a small tolerance tol to be also
                             % able to deal with not strictly complementary
                             % data - our default value for tol is 10^-10
   if (solves >= 1000); done = 1;       % emergency exit to avoid cycling
         if(sil==0)
            fprintf('max number of outer iterations reached. \n');
         end
         flag = 0;
   end;
          
end;                 % end while
maxfloor = max(maxfloor,floor);
    if(floor == 0)
        time1 = toc;
        avg_inact = avg_inact/solves;
    else
        time1 = 0;
    end
%k












