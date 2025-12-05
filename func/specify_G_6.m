function [supplyG,forward_star_rep,look_up] = specify_G_6(supply, veh_routes, zc, zd, W, T, yc, yd, inoutcap, storage_cap)
% this function specify_G creates a network G for time-varying VEN of
% with supply and forward_star_rep that specifies the supplies of each node in network G
% and arcs in network G respectively

% Note that the self-loop arcs for nodes with positive supplies are created
% in find_opt_flow_via_simplex_alg, not here

% It takes as input 
% 1)the supply of road junction in the road network
% 2) veh_routes
    % veh_routes is a cell structure that contains the vehicular routes and
    % flows on each route. It is a 3 by |R| structure where |R| is the number
    % of vehicular routes. Each route is represented as a cell in the top row
    % of veh_routes, with the cell containing a vector of the list of nodes
    % that represents the path of the route. The flow on that route is recorded
    % as a scalar or vector (for time-varying VEN) in the cell in the same column as the route in veh_routes
    % The third row contains the absolute arrival times at each node in the corresponding route.
% 3,4) zc, zd, which are the chargin and discharging effciencies respectively
% 5)W = amount of energy each EV can carry
% 6) T %number of seconds in each period = 5min
% 7-10) the charge and discharge efficiency of storage, the max power in and out of storage, the max storage  

num_period = size(cell2mat(veh_routes(2,1)),1);%number of period
R = size(veh_routes,2); % R = number of routes
n_D = size (supply,2); %number of road junctions

% look_up is a matrix that stores the numerical node assignment of the
% non-junction nodes that are created for G.
% Each non-junction node is a row in look_up.

% Each non-junction node has nomaclature J_m^n-n-m-t where
% n = the n^th route
% J_m^n = the junction node i that the non-junction node is connected
% to
% m = the number of node along the route (which is a path) r_n that
% J_m^n is
% t = time period associated with the non-junction node

% look_up has 5 columns: 1=J_m^n, 2=n , 3=m,4=t, 5=
% number assigned to the non-junction node

% look_up also stores the junction nodes at different time, i.e J_t, in which case
% the 5 columns are 1=J, 2=0, 3=0, 4=t, 5=number assigned

% look_up also stores the supersink and supersource nodes, in which case
% the 5 columns are 1=J, 2=0, 3=0, 4=0, 5=number assigned

% look_up also stores the nodes that represent storage at junctions, in which case
% the 5 columns are 1=J, 2=-1, 3=-1, 4=t, 5=number assigned

% first, fill look_up with junction nodes of different times
% assign numbers for the junction nodes at different times at look_up 
look_up = zeros(n_D*num_period,5);
look_up(:,5) = 1:n_D*num_period;
look_up(:,1) = repmat([1:n_D].',num_period,1); %i.e. 1234512345
for t = 1:num_period
    look_up(1+n_D*(t-1):n_D*t,4) = t;
end

%second, preallocate entries in look_up for non-junction nodes
num_assigned = n_D*num_period; %the number of nodes in G that has been assigned so far. 
% the numbering for those vehciualr routes nodes start at num_assigned+1

% calculate the totoal number of non-junction nodes at each period
summ = 0;
for n = 1:R %for the n^th route
    route = cell2mat(veh_routes(1,n));% the route considered in this iteration
    r_len = size(route,2);
    summ = summ+r_len;
end

look_up = [look_up; zeros(summ*num_period,5)];


%third, preallocate for nodes that represent storage at junctions.
%num_period of them for each junction node
look_up = [look_up; zeros(n_D*num_period,5)];

%size(look_up)
%% supplyG
% supply of nonjunction nodes are zero, supply of junctio nodes are the same
% as specified in supply
% create node for each source/ sink node
supplyG = zeros(1,size(look_up,1));
% the first n_D*num_period entries are given in supply
% the same junction node in different times are assigned numbers in look_up
% separated by n_D
for i = 1:n_D
    supplyG(1,[i:n_D:num_period*n_D]) = supply(:,i).';
end

supplyG = supplyG*T; %suppy has units of energy, since we now talk about energy demanded per period

%% arcs
% create the arcs in G and record them in forward_star_rep
% forward_star_rep is the matrix that stores all the arcs in G
% five columns are tail, head, cost, cap, gain

%preaooloate entries in forward_star_rep
%1.number of arcs linking between junction node and non-junction node in each period
summ = 0; %total number of nodes in all vehicular routes
for n = 1:R %for the n^th route
    route = cell2mat(veh_routes(1,n));% the route considered in this iteration
    r_len = size(route,2);
    summ = summ+ 2 + 2*(r_len-2);
end
forward_star_rep = zeros(summ*num_period,5);
%size(forward_star_rep)

B = inf; %capacity in many arcs

f_star_used = 0; %number of rows used

%2. preallocate arcs that link between non-junciton nodes across time
% total arcs to add
summ = 0; %total number of arcs in all vehicular routes
for n = 1:R %for the n^th route
    route = cell2mat(veh_routes(1,n));% the route considered in this iteration
    r_len = size(route,2)-1; %number of arcs
    summ = summ+ r_len;
end

forward_star_rep = [forward_star_rep; zeros(summ*(num_period-1),5)];


%3. also preallocate for arcs for storage. Note we add arcs that links the
%storage node of the last period to itself too. Other selflinking nodes for
%nodes with positive supplies are added in 'find_opt_flow_via_simplex_alg.
arcs_for_storage = n_D*(2*(num_period-1) + num_period-1)+n_D;
forward_star_rep = [forward_star_rep; zeros(arcs_for_storage,5)];

%size(forward_star_rep)

%% add arcs that link juncito nodes to nodes of vehicualr routes at each time period
for t = 1:num_period  
    
for n = 1:R %for the n^th route
    %disp(['------considering ' num2str(n) 'th route at time ' num2str(t) ' -----']);
    route = cell2mat(veh_routes(1,n));% the route considered in this iteration
    r_len = size(route,2); %length of the route, = no. of nodes visited
    start_node = look_up(all(look_up(:,1:4)==[route(1,1),0,0,t],2),5);%number assigned to start node
    end_node = look_up(all(look_up(:,1:4)==[route(1,end),0,0,t],2),5);%number assigned to end node
   
    
    %create two nodes corr to start and end of route
    look_up(num_assigned+1, :) = [route(1,1) n 1 t num_assigned+1];
    look_up(num_assigned+2, :) = [route(1,end) n r_len t num_assigned+2];

    %add directed arc to forward_star_rep
    forward_star_rep(f_star_used+1,:) = [start_node num_assigned+1 1-zc B zc];
    forward_star_rep(f_star_used+2,:) = [num_assigned+2 end_node 1-zd B zd];
    
    %we have added two arcs so adjust
    f_star_used = f_star_used+2;
    
    
    % adjust num_assigned since we have added two nodes
    num_assigned = num_assigned +2; 
    
    for m = 2:r_len-1
        %disp(['---' num2str(m) 'th node along route----']);
        junc_node = look_up(all(look_up(:,1:4)==[route(1,m),0,0,t],2),5); %the number assigned to J_t
        %create one nodes corr to each intermediate node in route
        look_up(num_assigned+1, :) = [route(1,m) n m t num_assigned+1];
        
        %create two arcs
        forward_star_rep(f_star_used+1,:) = [num_assigned+1 junc_node 1-zd B zd];
        forward_star_rep(f_star_used+2,:) = [junc_node num_assigned+1 1-zc B zc];
    
        %we have added two arcs so adjust
        f_star_used = f_star_used+2;
        
        % adjust num_assigned for use in next iteration since we have added two nodes
        num_assigned = num_assigned +1; 
        
    end
   
end

end

%size(look_up);
%size(forward_star_rep);

%% for connections bewteen non-junction nodes across time 
% Version 1
% for n = 1:R %for the n^th route
%     route = cell2mat(veh_routes(1,n));% the route considered in this iteration
%     r_len = size(route,2); %length of the route, = no. of nodes visited
%     flow_n = cell2mat(veh_routes(2,n));%flow on route n
%     T_n = cell2mat(veh_routes(3,n));
%     for m = 1:r_len-1
% 
%         tau = T_n(m+1)-T_n(m);
% 
%         for t = 1:num_period-tau
%             t1 = t+T_n(m)-1;
%             t2 = t+T_n(m+1)-1;
%             if t2 > num_period
%                 break
%             end
%             origin_node = look_up(all(look_up(:,1:4)==[route(1,m),n,m,t],2),5);
%             dest_node = look_up(all(look_up(:,1:4)==[route(1,m+1),n,m+1,t+tau],2),5);
%             cap = W*T*sum(flow_n( t1:t2));
%             % add arc        
%             forward_star_rep(f_star_used+1,:) = [origin_node dest_node 0 cap 1];
%     
%             %we have added 1 arc so adjust
%             f_star_used = f_star_used+1;
%         end
%         
%     end
% end

% Version 2
for r = 1:size(veh_routes,2) %for the n^th route
    routei = veh_routes{1,r};
    Ti = veh_routes{3,r};
    flowi = veh_routes{2,r};
    for t = 0:num_period - 1  
        for i = 1:length(routei)-1
            i1 = i;
            i2 = i+1;
            nodes = routei(i1);
            nodet = routei(i2);
            times = t + Ti(i1);
            timet = t + Ti(i2);

            if timet > num_period
                break
            end
            idx_nodes_xuni = look_up(all(look_up(:,1:4)==[nodes,r,i1,times],2),5);
            idx_nodet_xuni = look_up(all(look_up(:,1:4)==[nodet,r,i2,timet],2),5);
            capi = W*T*sum(flowi(times:timet));
            forward_star_rep(f_star_used+1,:) = [idx_nodes_xuni,idx_nodet_xuni,0,capi,1];
            f_star_used = f_star_used+1;
        end
    end

end


% Version 1


%% for arcs for storage
for i = 1:n_D
for t = 1:num_period
    %create storage node
    look_up(num_assigned+1,:) = [i, -1, -1, t, num_assigned+1];
    
    node_at_t = n_D*(t-1)+i; % the number assinged to junction node i at t
    
    if t==1
        forward_star_rep(f_star_used+1,:) = [node_at_t num_assigned+1 1-yc inoutcap*T yc];
        
        f_star_used = f_star_used+1;
    
    elseif t == num_period
        forward_star_rep(f_star_used+1,:) = [num_assigned+1 node_at_t 1-yd inoutcap*T yd];
        forward_star_rep(f_star_used+2,:) = [num_of_storage_at_prev_iteration num_assigned+1 0 storage_cap 1];
        forward_star_rep(f_star_used+3,:) = [num_assigned+1 num_assigned+1 0 inf 0.5]; %self-loop arc to allow energy to be left in storage
    
        f_star_used = f_star_used+3;
    else
        
        forward_star_rep(f_star_used+1,:) = [node_at_t num_assigned+1 1-yc inoutcap*T yc];
        forward_star_rep(f_star_used+2,:) = [num_assigned+1 node_at_t 1-yd inoutcap*T yd];
        forward_star_rep(f_star_used+3,:) = [num_of_storage_at_prev_iteration num_assigned+1 0 storage_cap 1];
        
        f_star_used = f_star_used+3;

    end
    
    num_of_storage_at_prev_iteration = num_assigned+1;
    num_assigned = num_assigned+1; % we have assigned one node in this iteration
    
end
end


forward_star_rep = forward_star_rep(~all(forward_star_rep == 0, 2), :);
end


%%
%{
figure(1)
G = digraph(forward_star_rep(:,1),forward_star_rep(:,2));
plot(G)

figure(2)

G2 = digraph(forward_star_rep([1:50,71:80],1),forward_star_rep([1:50,71:80],2));
plot(G2)

[optimal_flow, cost] = find_opt_flow_via_simplex_alg(supplyG, forward_star_rep);
%}