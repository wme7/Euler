function [node,elem,edge,bound] = BuildUnstructuredMesh2d(vx,vy,EtoV,nE,nN,BC)
% Build Unstructured Mesh 2-D
%
%   The following data, needed for NCFV method, will be constructed:
%       1. elem(:)
%       2. node(:)
%       3. edge(:)
%       4. bound(:)
%
%   NOTE: bc_name is the name of the boundary condition.
%         Only four BCs are available in this version:
%
%	1. "freestream"
%       Roe flux with freestream condition on the right state.
%
%	2. "slip_wall"
%       Solid wall condition. Mass flux through the boundary is set zero.
%
%   3. "outflow_supersonic"
%       Just compute the boundary flux by the physical Euler flux
%       (equivalent to the interior-extrapolation condition.)
%
%	4. "outflow_back_pressure"
%       Fix the back pressure. This should work for subsonic flows in a
%       large enough domain.

% 1. Read element-connectivity information

% Element quantities  : elem(:).x, elem(:).y, elem(:).vol
%
%  o-----------o            
%   \          |            o
%    \    (x,y)|           / \
%     \   .    |          /   \
%      \       |         /  .  \    (x,y): centroid coordinates
%       \      |        / (x,y) \     vol: volume of element
%        o-----o       o---------o
%
% Triangles and Quads: vertices are ordered counterclockwise
%
%         v3                    v4________v3
%         /\                     /        |
%        /  \                   /         |
%       /    \                 /          |
%      /      \               /           |
%     /        \             /            |
%    /__________\           /_____________|
%   v1           v2        v1             v2

% initialize elements
elem(nE).nV = 0;    % number of vertices
elem(nE).v = 0;     % vertices ids
elem(nE).x = 0;     % x-coordinate of the centroid
elem(nE).y = 0;     % y-coordinate of the centroid
elem(nE).A = 0;     % Area of the element (2d)
elem(nE).vol= 0;    % Volume of the element (3d)
elem(nE).nghbr=0;   % List of element neighbors of each element
elem(nE).nnghbrs=0; % Number of element neighbors of each element

nTRI=0; nQUAD=0;
for e = 1:nE
    elem(e).nV = numel(EtoV(e,:));
    switch elem(e).nV
        case 3; % Triangular element
            nTRI =nTRI +1;
            elem(e).v(1) = EtoV(e,1);
            elem(e).v(2) = EtoV(e,2);
            elem(e).v(3) = EtoV(e,3);
            ax = vx(EtoV(e,1)); ay = vy(EtoV(e,1));
            bx = vx(EtoV(e,2)); by = vy(EtoV(e,2));
            cx = vx(EtoV(e,3)); cy = vy(EtoV(e,3));
            Area = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
            elem(e).A = Area;
            elem(e).x = (ax+bx+cx)/3;
            elem(e).y = (ay+by+cy)/3;
            elem(e).vol=Area;
        case 4; % Quad element
            nQUAD=nQUAD+1;
            elem(e).v(1) = EtoV(e,1);
            elem(e).v(2) = EtoV(e,2);
            elem(e).v(3) = EtoV(e,3);
            elem(e).v(4) = EtoV(e,4);
            ax = vx(EtoV(e,1)); ay = vy(EtoV(e,1));
            bx = vx(EtoV(e,2)); by = vy(EtoV(e,2));
            cx = vx(EtoV(e,3)); cy = vy(EtoV(e,3));
            dx = vx(EtoV(e,4)); dy = vy(EtoV(e,4));
            Area = 0.5*(ax-dx)*(cy-dy)-(bx-dx)*(by-dy)-(cx-dx)*(ay-dy);
            elem(e).A = Area;
            elem(e).x = (ax+bx+cx+dx)/4;
            elem(e).y = (ay+by+cy+dy)/4;
            elem(e).vol=Area;
        otherwise
            error('element not set in options');
    end
end
%disp([nE,nE == (nTRI+nQUAD)]);

% Loop over elements and construct the following data.
%
% 1. Surrounding elements: node(:)%nelms, node(:)%elm(:)
%
%    Example: Node i is surrounded by the eleemnts, 23, 101, 13, 41.
%             node(i)%nelms = 4
%             node(i)%elm(1) = 23
%             node(i)%elm(2) = 13
%             node(i)%elm(3) = 41
%             node(i)%elm(4) = 101
%
%        o-------o-------------o
%       /        |   .         |
%      /    23   |      41     |
%     o----------o-------------o
%      \        i \            |
%       \   101    \     13    |
%        \          \          | 
%         o----------o---------o

% 2. Read Node data:

% initialize arrays
node(nN).nnghbrs=0; % Number of node neighbors of each node
node(nN).nghbr= 0;  % List of node neighbors of each node
node(nN).nelms= 0;  % Number of adjacent elements of each node
node(nN).elms = 0;  % List of adjacent elements of each node
node(nN).vol = 0;   % Dual volume around each node
node(nN).x = 0;     % x-coordinate of vertex
node(nN).y = 0;     % y-coordinate of vertex

% Initialize (counter) lists and node coordinates
for i = 1:nN
    node(i).nelms = 0;
    node(i).vol = 0;
    node(i).x = vx(i);
    node(i).y = vy(i);
end

% Distribute element indexes to nodes
for e = 1:nE
    elem(e).nV = numel(EtoV(e,:));
    switch elem(e).nV
        case 3
            v1=elem(e).v(1); node(v1).nelms=node(v1).nelms+1; node(v1).elms(node(v1).nelms)=e;
            v2=elem(e).v(2); node(v2).nelms=node(v2).nelms+1; node(v2).elms(node(v2).nelms)=e;
            v3=elem(e).v(3); node(v3).nelms=node(v3).nelms+1; node(v3).elms(node(v3).nelms)=e;
            % Dual volume around every node
            x1=node(v1).x; y1=node(v1).y; x2=node(v2).x; y2=node(v2).y;
            x3=node(v3).x; y3=node(v3).y; xc=elem( e).x; yc=elem( e).y;
            % v1
            node(v1).vol = node(v1).vol + elem(e).vol/3;
            % v2
            node(v2).vol = node(v2).vol + elem(e).vol/3;
            % v3
            node(v3).vol = node(v3).vol + elem(e).vol/3;
        case 4
            v1=elem(e).v(1); node(v1).nelms=node(v1).nelms+1; node(v1).elms(node(v1).nelms)=e;
            v2=elem(e).v(2); node(v2).nelms=node(v2).nelms+1; node(v2).elms(node(v2).nelms)=e;
            v3=elem(e).v(3); node(v3).nelms=node(v3).nelms+1; node(v3).elms(node(v3).nelms)=e;
            v4=elem(e).v(4); node(v4).nelms=node(v4).nelms+1; node(v4).elms(node(v4).nelms)=e;
            % Dual volume around every node
            x1=node(v1).x; y1=node(v1).y;       x2=node(v2).x; y2=node(v2).y; 
            x3=node(v3).x; y3=node(v3).y;       x4=node(v4).x; y4=node(v4).y; 
            x12=0.5*(x1+x2); y12=0.5*(y1+y2);   x23=0.5*(x2+x3); y23=0.5*(y2+y3);
            x34=0.5*(x3+x4); y34=0.5*(y3+y4);   x41=0.5*(x4+x1); y41=0.5*(y4+y1);
            xc=elem( e).x; yc=elem( e).y;
            % v1
            ax = x1; bx = x12; cx = xc; dx = x41;
            ay = y1; by = y12; cy = yc; dy = y41;
            node(v1).vol = node(v1).vol + 0.5*(ax-dx)*(cy-dy)-(bx-dx)*(by-dy)-(cx-dx)*(ay-dy);
            % v2
            ax = x12; bx = x2; cx = x23; dx = xc;
            ay = y12; by = y2; cy = y23; dy = yc;
            node(v2).vol = node(v2).vol + 0.5*(ax-dx)*(cy-dy)-(bx-dx)*(by-dy)-(cx-dx)*(ay-dy);
            % v3
            ax = xc; bx = x23; cx = x3; dx = x34;
            ay = yc; by = y23; cy = y3; dy = y34;
            node(v3).vol = node(v3).vol + 0.5*(ax-dx)*(cy-dy)-(bx-dx)*(by-dy)-(cx-dx)*(ay-dy);
            % v4
            ax = x41; bx = xc; cx = x34; dx = x4;
            ay = y41; by = yc; cy = y34; dy = y4;
            node(v4).vol = node(v4).vol + 0.5*(ax-dx)*(cy-dy)-(bx-dx)*(by-dy)-(cx-dx)*(ay-dy);
    end
end

% Loop over elements 2: to build element-neighbor data
%
%  Allocate elm(:)%nghbr(:) : elm(:)%nnghrs, elm(:)%nghr(:)
%  Construct element nghbr data: elm(:)%nghbr(:)
%  Order of neighbor elements [e1,e2,e3,..] are closely related to
%  the order of vertices [v1,v2,v3,..] (see below).
%
%          o------o
%          |      |                
%        v4|  e1  |v3                     v3
%    o-----o------o------o      o---------o------------o
%    |     |      |      |       .      .   .        .
%    | e2  |      |  e4  |        . e2 .     . e1  .
%    o-----o------o------o         .  .       .  .
%       v1 |     .v2              v1 o---------o v2   
%          | e3 .                     .   e3  .
%          |   .                        .    .
%          |  .                           . .
%          | .                             o
%          o
%
% we will loop around the edges and match the vertex of the neigbour
% elements. Start by selectin a single edge/face:
%
%   Get the face of the element i:
%              o------o
%             vr      rl
%
%             vL      vR
%              o------o
%             /       |
%            /    i   |
%           o---------o
%
% and compare with the neigbour elements faces. (note the ordening of
% vertex when comparing edges% ) 
for e = 1:nE
    for k = 1:elem(e).nV    % also number of total edges
        % evaluate edges in counter clockwise direction
        vR = elem(e).v(k); if k < elem(e).nV; vL=elem(e).v(k+1); else vL=elem(e).v(1); end
        %disp([k,vR,vL]) % Perform edge matching of: vR o---o vL
        found=false;     % Start with a negative condition
        for eAround = 1:node(vR).nelms % for all elements around vR
            j = node(vR).elms(eAround);
            for i = 1:elem(j).nV    % for every edge of j-element
                % evaluate edges in clockwise direction
                vr = elem(j).v(i); if i > 1; vl=elem(j).v(i-1); else vl=elem(j).v(elem(j).nV); end
                %disp([j,vr,vl]) % Perform edge matching of: vR o---o vL
                if (vR==vr && vL==vl); found=true; im=i+1;
                    if im > elem(j).nV; im=im-elem(j).nV; end
                    break
                end
            end
            if found; break; end
        end % elements around vR
        in = k+2;
        if in > elem(e).nV; in=in-elem(e).nV; end 
        if found; 
            elem(e).nghbr(in) = j;
            elem(j).nghbr(im) = e;
        else
           elem(e).nghbr(in) = 0; % no neighbour
        end
    end
end

% 3. Edge-data for node-centered (edge-based) scheme.
%
% Loop over elements 3: Construct edge data: edge(:).n1, n2, e1, e2.
% Edge points from node n1 to node n2.
%
%      n2
%       o------------o
%     .  \         .
%    .    \   e2  .
%   .  e1  \    .
%  .        \ .         Directed area is positive: n1 -> n2
% o----------o         e1: left element
%             n1       e2: right element (e2 > e1 or e2 = 0)
%
% First count the number of edges.
%
% NOTE: Count edges only if the neighbor element number is
%       greater than the current element (i) to avoid double
%       count. Zero element number indicates that it is outside
%       the domain (boundary face).
%
% Count the total amount of edges
nedges = 0;
for e = 1:nE
    nEdge = numel(EtoV(e,:));
    for i = 1:nEdge
    if elem(e).nghbr(i)>e || elem(e).nghbr(i)==0; nedges=nedges+1; end
    end
end
edge(nedges).nEdges= nedges; %disp(nedges)

% Buggy snipet CAREFUL%%
% nedges = 0;
% for e = 1:nE
%     nEdge = numel(EtoV(e,:));
%     switch nEdge
%         case 3
%             face = [elem(e).v(3),elem(e).v(1);
%                     elem(e).v(2),elem(e).v(3);
%                     elem(e).v(1),elem(e).v(2)];
%         case 4
%             face = [elem(e).v(4),elem(e).v(1);
%                     elem(e).v(3),elem(e).v(4);
%                     elem(e).v(2),elem(e).v(3);
%                     elem(e).v(1),elem(e).v(2)];
%     end
%     disp(face);
%     for i = 1:nEdge
%         if elem(e).nghbr(i) > e || elem(e).nghbr(i) == 0
%             nedges = nedges + 1;
%             edge(nedges).nodes = face(i,:);
%             disp([e,elem(e).nghbr(i)]);
%             edge(nedges).elems = [e,elem(e).nghbr(i)];
%         end
%     end
% end

% Build Edges nodes and elements list 
% edge(nedges).nodes = [0,0]; % [n1,n2] End nodes of each edge (edge points n1 -> n2)
% edge(nedges).elems = [0,0]; % [e1,e2] Left and right elements of each edge
nedges = 0;
for i = 1:nE
    % vertices
    v1 = elem(i).v(1);
    v2 = elem(i).v(2);
    v3 = elem(i).v(3);
    nEdge = numel(EtoV(i,:));
    switch nEdge
        case 3 % Triangle element
            if elem(i).nghbr(3) > i || elem(i).nghbr(3)==0
                nedges = nedges + 1;
                edge(nedges).nodes = [v1,v2];
                edge(nedges).elems = [i,elem(i).nghbr(3)];
            end
            if  elem(i).nghbr(1) > i || elem(i).nghbr(1)==0
                nedges = nedges + 1;
                edge(nedges).nodes = [v2,v3];
                edge(nedges).elems = [i,elem(i).nghbr(1)];
            end
            if  elem(i).nghbr(2) > i || elem(i).nghbr(2)==0
                nedges = nedges + 1;
                edge(nedges).nodes = [v3,v1];
                edge(nedges).elems = [i,elem(i).nghbr(2)];
            end
        case 4 % quadrilateral element
            v4 = elem(i).v(4);
            if  elem(i).nghbr(3) > i || elem(i).nghbr(3) ==0
                nedges = nedges + 1;
                edge(nedges).nodes = [v1,v2];
                edge(nedges).elems = [i,elem(i).nghbr(3)];
            end
            if  elem(i).nghbr(4) > i || elem(i).nghbr(4) ==0
                nedges = nedges + 1;
                edge(nedges).nodes = [v2,v3];
                edge(nedges).elems = [i,elem(i).nghbr(4)];
            end
            if  elem(i).nghbr(1) > i || elem(i).nghbr(1) ==0
                nedges = nedges + 1;
                edge(nedges).nodes = [v3,v4];
                edge(nedges).elems = [i,elem(i).nghbr(1)];
            end
            if  elem(i).nghbr(2) > i || elem(i).nghbr(2) ==0
                nedges = nedges + 1;
                edge(nedges).nodes = [v4,v1];
                edge(nedges).elems = [i,elem(i).nghbr(2)];
            end
    end
end

% Loop over edges: Construct edge vector and directed area vector.
%
% Edge vector is a simple vector pointing froom n1 to n2.
% For each edge, add the directed area vector (dav) from
% the left and right elements.
%
%              n2
%   o-----------o-----------o
%   |     dav   |  dav      |
%   |       ^   |   ^       |
%   |       |   |   |       |
%   |   c - - - m - - -c    |
%   |           |           |
%   |           |           |    m: edge midpoint
%   |           |           |    c: element centroid
%   o-----------o-----------o
%                n1
%
for d = 1:nedges
    % base information
    n1 = edge(d).nodes(1);
    n2 = edge(d).nodes(2);
    e1 = edge(d).elems(1);
    e2 = edge(d).elems(2);
    xm = 0.5*( node(n1).x + node(n2).x );
    ym = 0.5*( node(n1).y + node(n2).y );
    
    edge(d).dav = 0;
    
    % Contribution from the left element
    if (e1 > 0)
        xc = elem(e1).x;
        yc = elem(e1).y;
        edge(d).dav(1) = -(ym-yc);
        edge(d).dav(2) =  (xm-xc);
    end
    
    % Contribution from the right element
    if (e2 > 0)
        xc = elem(e2).x;
        yc = elem(e2).y;
        edge(d).dav(1) = edge(d).dav(1) -(yc-ym);
        edge(d).dav(2) = edge(d).dav(2) +(xc-xm);
    end
    
    % edge(:).dav = Unit directed area vector of each edge
    % edge(:).da  = Magnitude of the directed area vector for each edge
    edge(d).da  = sqrt( edge(d).dav(1)^2 + edge(d).dav(2)^2 );
    edge(d).dav = edge(d).dav ./ edge(d).da;
    
    % edge(:).ev  = Unit edge vector of each edge (vector n1 -> n2)
    % edge(:).e   = Magnitude of the edge vector for each edge
    edge(d).ev(1) = node(n2).x - node(n1).x; 
    edge(d).ev(2) = node(n2).y - node(n1).y;
    edge(d).e     = sqrt( edge(d).ev(1)^2 + edge(d).ev(2)^2 );
    edge(d).ev    = edge(d).ev ./ edge(d).e;
end

% Construct node neighbor data: pointers to the neighbor nodes(o)
%
%        o     o
%         \   / 
%          \ /
%     o-----*-----o
%          /|
%         / |
%        /  o        *: node in interest
%       o            o: neighbors (edge-connected nghbrs)
%
% Initialize count neighbour list
for i = 1:nN; 
    node(i).nnghbrs=0; 
end

% Loop over edges and distribute the node numbers:
for d = 1:nedges
    % Basic data
    n1 = edge(d).nodes(1);
    n2 = edge(d).nodes(2);
    
    % (1) Add node1 to the neighbor list of n2
    node(n1).nnghbrs = node(n1).nnghbrs + 1;
    node(n1).nghbr(node(n1).nnghbrs) = n2;
    
    % (2) Add node2 to the neighbor list of n1
    node(n2).nnghbrs = node(n2).nnghbrs + 1;
    node(n2).nghbr(node(n2).nnghbrs) = n1;
end

% 4. Boudnary data
%    bound(:).bnx    = Outward normal at boundary nodes (x-component of unit vector)
%    bound(:).bny    = Outward normal at boundary nodes (y-component of unit vector)
%    bound(:).bn     = Magnitude of (bnx,bny)
%    bound(:).bfnx   = Outward normal at boundary nodes (x-component of unit vector)
%    bound(:).bfny   = Outward normal at boundary nodes (y-component of unit vector)
%    bound(:).bfn    = Magnitude of (bfnx,bfny)
%    bound(:).belm   = Element to which the boundary face belongs

nbound = 5;	% Number of boundary segments
bound(1).nbnodes = numel(BC.inflow);
bound(2).nbnodes = numel(BC.leftwall);
bound(3).nbnodes = numel(BC.bottomOutflow);
bound(4).nbnodes = numel(BC.rightOutflow);
bound(5).nbnodes = numel(BC.topwall);

bound(1).bnode = flipud(BC.inflow);
bound(2).bnode = flipud(BC.leftwall);
bound(3).bnode = BC.bottomOutflow;
bound(4).bnode = BC.rightOutflow;
bound(5).bnode = flipud(BC.topwall);

bound(1).bname = 'freestream';
bound(2).bname = 'slip_wall';
bound(3).bname = 'outflow_supersonic';
bound(4).bname = 'outflow_supersonic';
bound(5).bname = 'slip_wall';

% Boundary normal at nodes constructed by accumulating the contribution
% from each boundary face normal. This vector will be used to enforce
% the tangency condition, for example.
%
%        Interior domain      /
%                            o
%                  .        /
%                  .       /
% --o-------o-------------o
%           j   |  .  |   j+1
%               v  .  v
%
%        Left half added to the node j, and
%       right half added to the node j+1.
%
% Allocate and initialize the normal vector arrays
bound(nbound).bnx = 0;
bound(nbound).bny = 0;
bound(nbound).bn  = 0;
for i = 1:nbound
    for j = 1:bound(i).nbnodes
        bound(i).bnx(j) = 0;
        bound(i).bny(j) = 0;
        bound(i).bn( j) = 0;
    end
end

% Compute the outward normals
for i = 1:nbound
    for j = 1:bound(i).nbnodes-1;
        x1 = node(bound(i).bnode( j )).x;
        y1 = node(bound(i).bnode( j )).y;
        x2 = node(bound(i).bnode(j+1)).x;
        y2 = node(bound(i).bnode(j+1)).y;
        bound(i).bnx( j ) = bound(i).bnx( j ) + 0.5*( -(y1-y2) );
        bound(i).bny( j ) = bound(i).bny( j ) + 0.5*(   x1-x2  );
        bound(i).bnx(j+1) = bound(i).bnx(j+1) + 0.5*( -(y1-y2) );
        bound(i).bny(j+1) = bound(i).bny(j+1) + 0.5*(   x1-x2  );
    end
end

% Compute the magnitude and turn (bnx,bny) into a unit vector
for i = 1:nbound
    for j = 1: bound(i).nbnodes
        bound(i).bn(j)  = sqrt( bound(i).bnx(j)^2 + bound(i).bny(j)^2 );
        bound(i).bnx(j) = bound(i).bnx(j) / bound(i).bn(j);
        bound(i).bny(j) = bound(i).bny(j) / bound(i).bn(j);
    end
end

% Boundary face data
%
%      |     Domain      |
%      |                 |
%      o--o--o--o--o--o--o  <- Boundary segment
%   j= 1  2  3  4  5  6  7
%
%   In the above case, nbnodes = 7, nbfaces = 6
%
for i = 1:nbound
    bound(i).nbfaces = bound(i).nbnodes-1;
    bound(i).bfnx( 1:bound(i).nbfaces ) = 0;
    bound(i).bfny( 1:bound(i).nbfaces ) = 0;
    bound(i).bfn ( 1:bound(i).nbfaces ) = 0;
    bound(i).belm( 1:bound(i).nbfaces ) = 0;
end

% Boundary face vector: outward normal
for i = 1:nbound
    for j = 1:bound(i).nbfaces
        
        x1 = node(bound(i).bnode( j )).x;
        y1 = node(bound(i).bnode( j )).y;
        x2 = node(bound(i).bnode(j+1)).x;
        y2 = node(bound(i).bnode(j+1)).y;
        
        bound(i).bfn(j)  =  sqrt( (x1-x2)^2 + (y1-y2)^2 );
        bound(i).bfnx(j) = -(y1-y2) / bound(i).bfn(j);
        bound(i).bfny(j) =  (x1-x2) / bound(i).bfn(j);
    end
end

% Boundary normal vector at nodes: outward normal
for i = 1:nbound
    for j = 1:bound(i).nbnodes-1
        
        x1 = node(bound(i).bnode(j  )).x;
        y1 = node(bound(i).bnode(j  )).y;
        x2 = node(bound(i).bnode(j+1)).x;
        y2 = node(bound(i).bnode(j+1)).y;
        
        bound(i).bfn(j)  =  sqrt( (x1-x2)^2 + (y1-y2)^2 );
        bound(i).bfnx(j) = -(y1-y2) / bound(i).bfn(j);
        bound(i).bfny(j) =  (x1-x2) / bound(i).bfn(j);
    end
end

% Find element adjacent to the face: belm
%
%  NOTE: This is useful to figure out what element
%        each boundary face belongs to. Boundary flux needs
%        special weighting depending on the element.
%
%      |_________|_________|________|
%      |         |         |        | 
%      |         |         |        | 
%      |_________|_________|________|
%      |         |         |        |     <- Grid (e.g., quads)
%      |         | elmb(j) |        |
%   ---o---------o---------o--------o---  <- Boundary segment
%                 j-th face
%
% elmb(j) is the element number of the element having the j-th boundary face.
%
for i = 1:nbound
    for j = 1:bound(i).nbfaces
        % bface is defined by the nodes v1 and v2.
        v1 = bound(i).bnode(j); v2 = bound(i).bnode(j+1); %disp([v1,v2]);
        
        found = false;
        
        % Find the element having the bface from the elements around the node v1.
        for k = 1:node(v1).nelms
            ielm = node(v1).elms(k);
            for ii = 1:elem(ielm).nV
                in = ii;
                im = ii+1;
                if im > elem(ielm).nV; im = im-elem(ielm).nV; end %return to 1
                vt1 = elem(ielm).v(in); vt2 = elem(ielm).v(im); %disp([vt1,vt2]);
                if (vt1 == v1 && vt2 == v2); found = true; break; end
            end
            
            if (found); break; end
        end
        if (found)
            bound(i).belm(j) = ielm;
        else
            error(' Boundary-adjacent element not found. Error...');
        end
    end
end

% Construct least-squares matrix for node-centered schemes.
%
%        o     o
%         \   / 
%          \ /
%     o-----*-----o
%          /|
%         / |
%        /  o        *: node in interest
%       o            o: neighbors (edge-connected nghbrs)
%

% Check the number of neighbor nodes (must have at least 2 neighbors)

  ave_nghbr = node(1).nnghbrs;
  min_nghbr = node(1).nnghbrs;
  max_nghbr = node(1).nnghbrs;
       imin = 1;
       imax = 1;

 for i = 2:nN
   ave_nghbr = ave_nghbr + node(i).nnghbrs;
   if (node(i).nnghbrs < min_nghbr); imin = i; end
   if (node(i).nnghbrs > max_nghbr); imax = i; end
   min_nghbr = min(min_nghbr, node(i).nnghbrs);
   max_nghbr = max(max_nghbr, node(i).nnghbrs);
end

  fprintf('  ave_nghbr = %g\n',ave_nghbr/nN);
  fprintf('  min_nghbr = %g at node %g\n' ,min_nghbr,imin);
  fprintf('  max_nghbr = %g ay node %g\n' ,max_nghbr,imax);
  fprintf('\n');

% Now, compute the inverse of the LSQ matrix at each node.

for i = 1:nN
    node(i).invAtA = LSQinvMat2d(node(i).x,node(i).y,i,...
        [node(node(i).nghbr).x],[node(node(i).nghbr).y],node(i).nghbr);
end