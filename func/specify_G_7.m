
function [supplyG ,forward_star_rep, look_up] = specify_G_7(supply, veh_routes, zc, zd, W, T, yc, yd, inoutcap, storage_cap)
 
%% 1. 初始化
num_period = size(cell2mat(veh_routes(2,1)),1);%number of period
R = size(veh_routes,2); % R = number of routes
n_D = size (supply,2); %number of road junctions


% 1. look_up
look_up = zeros(n_D*num_period,5);
look_up(:,5) = 1:n_D*num_period;
look_up(:,1) = repmat([1:n_D].',num_period,1); %i.e. 1234512345
for t = 1:num_period
    look_up(1+n_D*(t-1):n_D*t,4) = t;
end

summ = 0;
for n = 1:R %for the n^th route
    route = cell2mat(veh_routes(1,n));% the route considered in this iteration
    r_len = size(route,2);
    summ = summ+r_len;
end
look_up = [look_up; zeros(summ*num_period,5)];
look_up = [look_up; zeros(n_D*num_period,5)];


% 2. forward_star_rep
summ = 0; %total number of nodes in all vehicular routes
for n = 1:R %for the n^th route
    route = cell2mat(veh_routes(1,n));% the route considered in this iteration
    r_len = size(route,2);
    summ = summ+ 2 + 2*(r_len-2);
end
forward_star_rep = zeros(summ*num_period,5);
summ = 0; %total number of arcs in all vehicular routes
for n = 1:R %for the n^th route
    route = cell2mat(veh_routes(1,n));% the route considered in this iteration
    r_len = size(route,2)-1; %number of arcs
    summ = summ+ r_len;
end
forward_star_rep = [forward_star_rep; zeros(summ*(num_period-1),5)];
arcs_for_storage = n_D*(2*(num_period-1) + num_period-1)+n_D;
forward_star_rep = [forward_star_rep; zeros(arcs_for_storage,5)];


f_star_used = 0; 
l_up_used = n_D*num_period;  %
%% 2. 
B = inf;
for r = 1:R
    % 1. routei: 
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
            
            if i == 1
                look_up(l_up_used+1,:)= [nodes,r,i1,times,l_up_used+1];
                look_up(l_up_used+2,:)= [nodet,r,i2,timet,l_up_used+2];
                l_up_used = l_up_used  + 2;
            else
                look_up(l_up_used+1,:) = [nodet,r,i2,timet,l_up_used+1];
                l_up_used = l_up_used  + 1;
            end

            idx_nodes_xuni = look_up(all(look_up(:,1:4)==[nodes,r,i1,times],2),5);
            idx_nodet_xuni = look_up(all(look_up(:,1:4)==[nodet,r,i2,timet],2),5);           
            idx_nodes_zhenshi = look_up(all(look_up(:,1:4)==[nodes,0,0,times],2),5);
            idx_nodet_zhenshi = look_up(all(look_up(:,1:4)==[nodet,0,0,timet],2),5);
            

            % 1) nodes
            if i1 == 1
                forward_star_rep(f_star_used+1,:) = [idx_nodes_zhenshi,idx_nodes_xuni,1-zc,B,zc];
                f_star_used  = f_star_used+1;
            else
                forward_star_rep(f_star_used+1,:) =  [idx_nodes_zhenshi,idx_nodes_xuni,1-zc,B,zc];
                forward_star_rep(f_star_used+2,:) =  [idx_nodes_xuni,idx_nodes_zhenshi,1-zc,B,zc];
                f_star_used  = f_star_used+2;
            end

            % 2) nodet
            if i2 == length(routei)
                forward_star_rep(f_star_used+1,:) =  [idx_nodet_xuni,idx_nodet_zhenshi,1-zc,B,zc];
                f_star_used  = f_star_used+1;
            else
                forward_star_rep(f_star_used+1,:) =  [idx_nodet_zhenshi,idx_nodet_xuni,1-zc,B,zc];
                forward_star_rep(f_star_used+2,:) = [idx_nodet_xuni,idx_nodet_zhenshi,1-zc,B,zc];
                f_star_used  = f_star_used+2;
            end

            % 3) nodes --> nodet
            capi = W*T*sum(flowi(times:timet));
            forward_star_rep(f_star_used+1,:) = [idx_nodes_xuni,idx_nodet_xuni,0,capi,1];
            f_star_used  = f_star_used+1;

        end
    end
end


%% for arcs for storage
for i = 1:n_D
for t = 1:num_period
    %create storage node
    look_up(l_up_used+1,:) = [i, -1, -1, t,l_up_used+1] ;
    node_at_t = n_D*(t-1)+i; % the number assinged to junction node i at t

    if t==1
        forward_star_rep(f_star_used+1,:) = [node_at_t l_up_used+1 1-yc inoutcap*T yc];
        
    elseif t == num_period
        forward_star_rep(f_star_used+1,:) = [l_up_used+1 node_at_t 1-yd inoutcap*T yd];
        forward_star_rep(f_star_used+2,:) = [num_of_storage_at_prev_iteration l_up_used+1 0 storage_cap 1];
        forward_star_rep(f_star_used+3,:) = [l_up_used+1 l_up_used+1 0 inf 0.5]; %self-loop arc to allow energy to be left in storage
    
        f_star_used = f_star_used+3;
    else
        forward_star_rep(f_star_used+1,:) = [node_at_t l_up_used+1 1-yc inoutcap*T yc];
        forward_star_rep(f_star_used+2,:) = [l_up_used+1 node_at_t 1-yd inoutcap*T yd];
        forward_star_rep(f_star_used+3,:) = [num_of_storage_at_prev_iteration l_up_used+1 0 storage_cap 1];
        
        f_star_used = f_star_used+3;
    end
    num_of_storage_at_prev_iteration = l_up_used+1;
    l_up_used = l_up_used+1; % we have assigned one node in this iteration
    
end
end

forward_star_rep = forward_star_rep(~all(forward_star_rep == 0, 2), :);
look_up = look_up(~all(look_up == 0, 2), :);

%% 3. supplyG
n_nodes = size(look_up,1);
supplyG = zeros(1,n_nodes);

for i = 1:n_D
    supplyG(1,[i:n_D:num_period*n_D]) = supply(:,i).';
end
supplyG = supplyG*T;

end


