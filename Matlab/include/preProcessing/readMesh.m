
function msh = readMesh(filename,meshType)
fprintf(1, 'READ MESH:      '); tic
msh.Types = { ...
    { 2, 1, 'LINES', 'nbLines'}, ... % 1
    { 3,  2, 'TRIANGLES', 'nbTriangles'}, ...
    { 4,  2, 'QUADS', 'nbQuads'}, ...  
    { 4,  3, 'TETS', 'nbTets'}, ...
    { 8,  3, 'HEXAS', 'nbHexas'}, ... %5
    { 6,  3, 'PRISMS', 'nbPrisms'}, ...
    { 5,  3, 'PYRAMIDS', 'nbPyramids'}, ...
    { 3,  1, 'LINES3', 'nbLines3'}, ...
    { 6,  2, 'TRIANGLES6', 'nbTriangles6'}, ...
    { 9,  2, 'QUADS9', 'nbQuads9'}, ... % 10
    { 10,  3, 'TETS10', 'nbTets10'}, ...
    { 27,  3, 'HEXAS27', 'nbHexas27'}, ...
    { 18,  3, 'PRISMS18', 'nbPrisms18'}, ...
    { 14,  3, 'PYRAMIDS14', 'nbPyramids14'}, ...
    { 1,  0, 'POINTS', 'nbPoints'}, ... % 15
    { 8,  3, 'QUADS8', 'nbQuads8'}, ...
    { 20,  3, 'HEXAS20', 'nbHexas20'}, ...
    { 15,  3, 'PRISMS15', 'nbPrisms15'}, ...
    { 13,  3, 'PYRAMIDS13', 'nbPyramids13'}, ...
};

switch meshType
    
    case 'QUADS'
        
        msh.ndim = 2;
        mesh_num = 3;
        
        
    case 'HEXAS'
        
        msh.ndim = 3;
        mesh_num = 5;
   
        
end

fid = fopen(filename, 'r');



% Find Nodes
iline = 0; findNodes = 0;

while findNodes == 0
   iline = iline + 1; 
   tline = fgetl(fid);
   if strcmp(tline, '$Nodes')
       findNodes = 1;
   end
end

nnodes = fscanf(fid, '%d', 1);
msh.nnode = nnodes;
tmp = fscanf(fid, '%g', [4 nnodes]);
msh.coords = tmp(2:4,:)';
fgetl(fid); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tline = fgetl(fid); tline = fgetl(fid);

% %%%%%%%%%%%%%%%%%% Find Elements %%%%%%%%%%%%%%%%%%%%%%%%
numElem = fscanf(fid, '%d' ,1);

tmp = fscanf(fid, '%g', inf); % Read it all in

cnt = 1;

elem_info.connectivity = zeros(numElem,msh.Types{mesh_num}{1});

for ie = 1:numElem
    elem_info.type(ie) = tmp(cnt+1);
    elem_info.tags(ie) = tmp(cnt+2);
    elem_info.region(ie) = tmp(cnt+4);
    elem_info.numNodes(ie) = msh.Types{elem_info.type(ie)}{1};
    first_node = cnt+2+elem_info.tags(ie)+1;
    last_node = cnt+2+elem_info.tags(ie) + elem_info.numNodes(ie);
    elem_info.connectivity(ie,1:elem_info.numNodes(ie)) = tmp(first_node:last_node);
    cnt = last_node + 1;
end

id = find(elem_info.type == mesh_num);

msh.nelem = length(id);


for je = 1:msh.nnode
    msh.nodes{je}.elements = [];
end
    

for ie =1:msh.nelem
    
    msh.elements{ie}.type = elem_info.type(id(ie));
    msh.elements{ie}.region = elem_info.region(id(ie));
    msh.elements{ie}.connectivity = elem_info.connectivity(id(ie),:);
    
    for j = 1:elem_info.numNodes(id(ie))
        je = msh.elements{ie}.connectivity(j);       
        msh.nodes{je}.elements = union(msh.nodes{je}.elements,ie);
    end
    
end

msh.nedof = msh.Types{mesh_num}{1}*msh.ndim;
msh.enode = msh.Types{mesh_num}{1};

msh.e2g = zeros(msh.nelem,msh.nedof);

for ie = 1:msh.nelem

    for j = 1:msh.ndim
        msh.e2g(ie,(j-1)*msh.enode + 1:j*msh.enode) = (j-1)*msh.nnode + msh.elements{ie}.connectivity;
    end

end

msh.tdof = msh.ndim*msh.nnode;
fprintf(1, [num2str(toc,'%8.6f'),'\n']);
end