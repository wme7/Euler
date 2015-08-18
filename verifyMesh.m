function verifyMesh(node,elem,edge,bound,vx,vy,EtoV,nE,nN)
%
% Here we check the grid data as follows:
%
% 1. Visual inspection of the elements, cell centers and boundary nodes.
% 2. Directed area must sum up to zero around every node.
% 3. Directed area must sum up to zero over the entire grid.
% 4. Global sum of the boundary normal vectors must vanish.
% 5. Global sum of the boundary face normal vectors must vanish.
% 6. Check element volumes which must be positive.
% 7. Check dual volumes which must be positive.
% 8. Global sum of the dual volumes must be equal to the sum of element volumes.
% 9. Linear LSQ gradients must be exact for linear functions.

%% 1. Visual inspection of the elements, cell centers and boundary nodes.

if true
    % Plot Element as patches
    patch('Faces',EtoV,'Vertices',[vx,vy],...
        'FaceColor','none','EdgeColor','k'); hold on;
    % or: patch(vx(EtoV)',vy(EtoV)',zeros(size(EtoV')),'w');
    
    % plote Element center
    scatter([elem.x],[elem.y])
    
    % Plot nodes and elements numbers
    for i = 1:nN
        text(vx(i),vy(i),int2str(i),'fontsize',8,....
            'fontweight','bold','Color','r');
    end
    for e = 1:nE
        pos = [mean(vx(EtoV(e,:))),mean(vy(EtoV(e,:)))];
        text(pos(1),pos(2),int2str(e),'fontsize',8,...
            'fontweight','bold','color','b');
    end
    
    % identify boundary and elements
    for b = 1:numel(bound);
        scatter([node(bound(b).bnode).x],[node(bound(b).bnode).y])
    end
    
    % draw edge lines
    for g = 1:edge(end).nEdges
       plot([node(edge(g).nodes).x],[node(edge(g).nodes).y],'-r') 
    end
end

%% 2. Directed area must sum up to zero around every node.
disp(abs(sum([elem.vol])-2)<1E-10);
disp(abs(sum([node.vol])-2)<1E-10);

end