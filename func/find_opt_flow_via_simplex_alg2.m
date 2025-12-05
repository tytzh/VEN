
function [optimal_flow, cost,Aeq,beq,lb,ub,f,t_linprog] = find_opt_flow_via_simplex_alg2(supply, forward_star_rep,idx_matched)
% uses matlab's linear programming solver to find the optimal flow
% given the supply matrix and the forward_star_rep fo the network without
% the arcs (i,i) for the supply nodes

%%
n = size(supply,2); % number of nodes in G
original_m = size(forward_star_rep,1); %number of arcs in G before adding (i,i) for supply nodes

%% add arc (i,i) to forward_star_rep for suply nodes i
value_arc = 0.5;
for i = 1:n
    if supply(i) > 0
        forward_star_rep = [forward_star_rep; i i 0 inf value_arc]; 
    end
end
m = size(forward_star_rep,1); %number of arcs in G

%% create incidence matrix
incidence = sparse(n,m); %use sparse matrix to save memory, this is a n by m zero matrix
for col=1:m %modify each column in turn
    arc = forward_star_rep(col,[1,2]); %the arc represented by the column
    i = arc(1,1); %head
    j = arc(1,2); %tail
    mu = forward_star_rep(col,5); % the multiplier of the arc (i,j)
    incidence(i,col) = 1; % if i is the head, incidence = 1
    
    if  col <= original_m && i~=j
        incidence(j,col) = -mu; % 
    else
        incidence(j,col) = -mu+1; %arc (i,i) requires we add 1
    end
end


%% solve by linear programming tool
A = [];
b = [];
Aeq = incidence;
beq = supply.';
lb = zeros(m,1);
ub = forward_star_rep(:,4); %capacity
ub(idx_matched) = 0;
f = forward_star_rep(:,3); %cost
options = optimoptions('linprog', 'Display', 'none');


%% slack
idx_slack_rows = find(beq~=0);
[n1,n2] = size(Aeq);
k = length(idx_slack_rows); 
Aeq_slack = zeros(n1, k);

for i = 1:k
    Aeq_slack(idx_slack_rows(i), i) = -1; 
end
Aeq_new = [Aeq, Aeq_slack];
lb_new = [lb; zeros(k,1)];
ub_new = [ub; inf(k,1)];
penalty = 10;  
f_new = [f; penalty * ones(k,1)];

%% 求解
tic;
% [x, cost] = linprog(f,A,b,Aeq,beq,lb,ub,options);
[xs_opt, cost] = linprog(f_new, A, b, Aeq_new, beq, lb_new, ub_new,options);
t_linprog = toc;
%% put x to optimal_flow flow
x_opt     = xs_opt(1:n2); 
slack_opt = xs_opt(n2+1:end);
optimal_flow = [forward_star_rep(:,[1,2]) x_opt];

end