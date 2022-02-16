function tree = volume_filling_tree(xmin,xmax,ymin,ymax,zmin,zmax,n_seed)
%INPUTS

%FOR NOW

% xmin,xmax,ymin,ymax,zmin,zmax : define the bounds to grow the tree within
% n_seed : number of seedpoints 

%TO BE ADDED AS OPTIONS - JUST DEFINING FOR NOW

stem_pos=[500 500 0];


%options for tree growing
%number of generations to (try and) grow
n_gen=12;
%the ratio of the distance to the centre of mass that the diameter grows
branch_distance_ratio=0.5;
%maximum branching angle
max_branch_angle=100;

d_ratio = 0.8;

d_stem = 10;

%stem position
%stem diameter 
%max number of generations
%branch distance ratio - distance to the CoM that the branch grows
%max branch angle 
%diameter ratio (could do this later with a trees toolbox function)
%capillary adding options

%OUTPUT
%tree in trees toolbox format
% X, Y, Z : coordinates of the nodes
% dA : adjacency matrix
%
%
%
%
%


% DO FOR X,Y,Z,dA FIRST - SIMPLEST PARAMETERS!

%initialise arrays using the number of seed points 
%actual size will be smaller but can trim at the end
%(I think!) the maximum number of nodes is 2*n_seed
max_nodes = 2*n_seed;
%node positions in X,Y,Z form
X = zeros(max_nodes,1);
Y = zeros(max_nodes,1);
Z = zeros(max_nodes,1);
%node positions in vector form
P = zeros(max_nodes,3);
%initialise adjancency matrix
dA = sparse(max_nodes, max_nodes);

X(1) = stem_pos(1);
Y(1) = stem_pos(2);
Z(1) = stem_pos(3);
P = [X Y Z];


%counter encoding the current number of nodes
n_nodes = 1;

%set up the seed points
%X,Y,Z form
Xseed = xmin + rand(n_seed,1).*(xmax-xmin);
Yseed = ymin + rand(n_seed,1).*(ymax-ymin);
Zseed = zmin + rand(n_seed,1).*(zmax-zmin);
%vector form
Pseed = [Xseed Yseed Zseed];

%matrix storing which node each seed point is associated with 
seed_node_index = zeros(n_gen,n_seed);
seed_node_index(1,:) = ones(1,n_seed);

%matrix storing which side of the plane (defined by the branch and COM of seed points)
%each seed point lies on
seed_side_index = zeros(n_gen,n_seed);
%seed_side_index(1,:) 

%matrix storing which nodes are terminal
terminal_node_index = zeros(max_nodes,1);

%calculate the centre of mass of the seed points
seed_COM = [mean(Xseed) mean(Yseed) mean(Zseed)];

%first branch grows towards this centre of mass
%first node
P1 = [X(1) Y(1) Z(1)];
%vector between first node and centre of mass
V = seed_COM - P1;
%normalised vector
Vnorm = V / norm(V);

%second node is the end of the first branch
P2 = P1 + branch_distance_ratio*V;
%update node counter
n_nodes = n_nodes + 1;
%update X,Y,Z form
X(2) = P2(1);
Y(2) = P2(2);
Z(2) = P2(3);
%update vector form
P(2,:) = [X(2) Y(2) Z(2)];
%the two nodes are connected
dA(2,1) = 1;


% X(2) = X(1) + branch_distance_ratio * Xcom;
% Y(2) = Y(1) + branch_distance_ratio * Ycom;
% Z(2) = Z(1) + branch_distance_ratio * Zcom;

%store the current tree in trees toolbox format
tree.X = X(1:2);
tree.Y = Y(1:2);
tree.Z = Z(1:2);
tree.dA = dA(1:2,1:2);






%set up for growing the next generation of branches

%define a plane containing P1, P2, and COM of the seed points
plane_Vnorm = vector_normal_to_3_points(P1, P2, seed_COM);
%calculate on which side of the plane the seed points lie
for i=1:size(Xseed,1)
    seed_side_index(2,i) = check_which_side_of_plane(Pseed(i,:),plane_Vnorm,seed_COM); %THIS CAN BE VECTORISED!
end
%these seed points are associated with node 2
seed_node_index(2,:) = 2 * ones(1,n_seed);


%grow the branches

for gen = 2:n_gen
    %find NODES at the previous generation that: i) have no child nodes,
    %ii) are not terminal nodes    
    child_node_index = child_tree(tree);
    end_nodes = find(child_node_index(:) == 0 & terminal_node_index(1:n_nodes) == 0)';    
     
    for end_node_index = end_nodes %loop through end_nodes            
        %get the parent node index
        parent_node_index = find(tree.dA(end_node_index,:));
        %calculate the CoM of the seed points associated with this end node 
        seed_COM = mean(Pseed(seed_node_index(gen,:) == end_node_index, :));
        %define a plane which contains: 
        %the end node, its parent node, and the CoM of the associated seed points       
        %plane_Vnorm = vector_normal_to_3_points(P(end_node_index,:), P(parent_node_index,:), seed_COM);
        %define a plane using two vectors:
        %1)vector between current node's parent and COM of seed points
        %2)vector between current node and COM of seed points
        v1 = seed_COM - P(parent_node_index,:) ;
        v2 = seed_COM - P(end_node_index) ;
        plane_Vnorm = vector_normal_to_2_vectors(v1,v2);         
        %calculate on which side of the plane the seed points lie
        this_node_seed_points = find(seed_node_index(gen,:) == end_node_index);
        for i = this_node_seed_points
            seed_side_index(gen,i) = 1 + check_which_side_of_plane(Pseed(i,:),plane_Vnorm,seed_COM); %THIS CAN BE VECTORISED!
        end
        %now grow two nodes from the current node
        %each node grows towards the CoM of seed points on one side of the
        %plane
        COM_side = zeros(2,3);
        for side = 1:2
            %find all seed points that are: (i) in this region, and 
            %(ii) on this side of the plane               
            seed_region_and_side_flag = seed_node_index(gen,:)==end_node_index & seed_side_index(gen,:)==side;                  
            %if there are no seed points left to grow towards, make this
            %node a terminal node
            if sum(seed_region_and_side_flag) == 0
                terminal_node_index(n_nodes) = 1;
            else %otherwise can grow a new node
                %calculate the CoM of these seed points
                COM_side(side,:) = mean(Pseed(seed_region_and_side_flag, :));            
                %vector between this end node and this centre of mass
                V = COM_side(side,:) - P(end_node_index,:);
                %grow the new node
                new_node = P(end_node_index,:) + branch_distance_ratio*V;                         
                %update node counter
                n_nodes = n_nodes + 1;               
                %put into X,Y,Z
                X(n_nodes) = new_node(1);
                Y(n_nodes) = new_node(2);
                Z(n_nodes) = new_node(3);
                %put into vector format
                P(n_nodes,:) = [X(n_nodes) Y(n_nodes) Z(n_nodes)];
                %update the adjacency matrix
                dA(n_nodes,end_node_index) = 1;
                %store which seed points are associated with this node
                seed_node_index(gen+1,seed_region_and_side_flag) = n_nodes;

                %if there is only one seed point left in this region
                %then make the new node a terminal node
                if sum(seed_region_and_side_flag) == 1
                    terminal_node_index(n_nodes) = 1;
                end

            end
         
        end
        

    end

    %trim empty entries and store the current tree in trees toolbox format
    tree.X = X(1:n_nodes);
    tree.Y = Y(1:n_nodes);
    tree.Z = Z(1:n_nodes);   
    tree.dA = dA(1:n_nodes,1:n_nodes);
end


%assign diameters to the tree
D=zeros(size(X(1:n_nodes)));
%calculate the branching order
BO = BO_tree(tree);
for gen = unique(BO(:))'
    %assign the diameters for this generation - diameter reduces by a
    %factor of d_ratio for each generation    
    D(BO == gen) = (d_ratio^gen) * d_stem ;
end


tree.D = D;


%trim empty entries and store the final tree in trees toolbox format
tree.X = X(1:n_nodes);
tree.Y = Y(1:n_nodes);
tree.Z = Z(1:n_nodes);
tree.P = P(1:n_nodes,:);
tree.dA = dA(1:n_nodes,1:n_nodes);
%store the seed point index arrays
tree.seed_node_index = seed_node_index;
tree.seed_side_index = seed_side_index;
tree.Xseed = Xseed;
tree.Yseed = Yseed;
tree.Zseed = Zseed;
tree.Pseed = Pseed;
tree.terminal_node_index = terminal_node_index(1:n_nodes);







