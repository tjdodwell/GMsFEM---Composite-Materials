% need to hand code rigid body modes

%  == Code Visualisation of Partition of Unity
clear all;
addpath('include/metis-5.1.0/metismex/');
addpath('include/postProcessing/');
addpath('inhull/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Model Setup 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numEigs = 100; % number of eigenvalues in excess of 3 rigid transalations
sizeOverlap = 1;

L = 5; % ensure this agrees with the geo file

useMetis = false;
include_defect = true;
plotxi = false;

test_code = true;

BMatrix = 'GenEO'; %'GenEO'; % other option is BMatrix = 'GenEO'
partition_of_unity_method = 'SF'; % SF for shape functions
                                  % (partition of unity done by
                                  % linear shape functions of macro grid

verb = 1;

filename_coarse = 'layercake_coarse.msh'; % coarse mesh
filename = 'layercake_fine.msh'; % fine scale mesh
numLayers = 1;
ss = [0.0,-1,0.0,-1,0.0,-1,0.0];

layer_colours = [2.5] * ones(7,1);

% E_R = 10; %    GPa
% nu_R = 0.35;

%AS4/8552
% E1 = 162; %    GPa
% E2 = 10; %    GPa
% E3 = 10; %    GPa

% nu_21 = 0.35; 
% nu_31 = 0.35;
% nu_32 = 0.5;

% G_12 = 5.2;   %   GPa
% G_13 = 5.2;   %   GPa
% G_23 = 3.5; %   GPa

E_R = 5.0; %    GPa
nu_R = 0.4;

E1 = 135; %    GPa
E2 = 5.0; %    GPa
E3 = 5.0; %    GPa

nu_21 = 0.022; 
nu_31 = 0.022;
nu_32 = 0.4;

G_12 = 5.01;   %   GPa
G_13 = 5.01;   %   GPa
G_23 = 5.01; %   GPa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[isotropic,composite] = makeMaterials(E_R,nu_R,E1,E2,E3,nu_21,nu_31,nu_32,G_12,G_13,G_23);

mshC= readMesh(filename_coarse,'HEXAS');
msh = readMesh(filename,'HEXAS');

numParts = mshC.nnode; % must be odd to include defect in middle part

% msh_coarse = defineIPs(msh_coarse);
msh = defineIPs(msh);

% [msh_coarse.N, msh_coarse.dN] = ShapeFunctions(msh_coarse);
[msh.N,msh.dN] = ShapeFunctions(msh);

msh.lhs = find(msh.coords(:,1) < 1e-6);
msh.rhs = find(msh.coords(:,1) > L - 1e-6);

mshC.lhs = find(mshC.coords(:,1) < 1e-6);
mshC.rhs = find(mshC.coords(:,1) > L - 1e-6);

%%%% Create Overlapping Partition

msh.Omg = cell(numParts,1);
msh.Omg_p = cell(numParts,1);

for ie = 1 : mshC.nelem

  X = zeros(8);

  for i = 1 : 8
      x = mshC.coords(mshC.elements{ie}.connectivity(i),:);
      X(i,:) = [1, x(1), x(2), x(3), x(1)*x(2), x(1)*x(3), x(2)*x(3), x(1)*x(2)*x(3)];
  end

  mshC.elements{ie}.invX = inv(X);
  
  mshC.elements{ie}.tess = convhulln(mshC.coords(mshC.elements{ie}.connectivity,:));

end



% == Partition using linear shape functions

for ip = 1 : numParts
	msh.Omg{ip}.map = zeros(msh.nelem,1);
end

for ie = 1 : msh.nelem

	xm = mean(msh.coords(msh.elements{ie}.connectivity,:)); % Find the midpoint

	for i = 1 : mshC.nelem % For each element - find the first macro element it belongs
        if(inhull(xm,mshC.coords(mshC.elements{i}.connectivity,:),mshC.elements{i}.tess,1e-5))
                macroElement = i;
                break
        end                
    end


    for ip = 1 : numParts % Find out which partition this macro element belongs
    	id = find(macroElement == mshC.nodes{ip}.elements);
    	if (length(id) > 0) % Record it in that parition
    		msh.Omg{ip}.map(ie) = 1;
    	end
    end

end

tmp = 1 : msh.nelem;
for ip = 1 : numParts
	msh.Omg{ip}.elems = tmp(logical(msh.Omg{ip}.map)); msh.Omg{ip}.map = [];
	msh.Omg{ip}.nelem = length(msh.Omg{ip}.elems);
end


% ===  Construct Overlapping Partition of m-layers

M = sparse(msh.nnode,msh.nelem);

for ie = 1 : msh.nelem
    for j = msh.elements{ie}.connectivity
        M(j,ie) = 1;
    end  
end

R = cell(numParts,1);

for ip = 1 : numParts % For each partion
    
    %msh.Omg{ip}.elems = msh.Omg_p{ip}.elements; % Record non-overlapping elements

    for i = 1 : sizeOverlap    % For each layer of overlap
        e = zeros(msh.nelem,1); % Initial element marker to zero
        e(msh.Omg{ip}.elems) = ones(length(msh.Omg{ip}.elems),1);
        tmp_elements = msh.Omg{ip}.elems;
        n = M * e;
        msh.Omg{ip}.nodes = find(n>0);
        n(msh.Omg{ip}.nodes) = ones(length(msh.Omg{ip}.nodes),1); 
        e = M'*n;
        msh.Omg{ip}.elems = find(e > 0);
    end

    msh.Omg{ip}.dof = [msh.Omg{ip}.nodes; msh.nnode + msh.Omg{ip}.nodes; 2 * msh.nnode + msh.Omg{ip}.nodes]; % nodes not included on processor boundary
    
    n = M * e;
    msh.Omg{ip}.nodes = find(n>0); % include nodes on processor
                                   % boundary - used for Neumann Matrix

    % === Overlapping process complete for patition Omega_ip

    msh.Omg{ip}.nelem = length(msh.Omg{ip}.elems);
    msh.Omg{ip}.nnode = length(msh.Omg{ip}.nodes);
    msh.Omg{ip}.coords = msh.coords(msh.Omg{ip}.nodes,:);
    msh.Omg{ip}.dofbar = [msh.Omg{ip}.nodes; msh.nnode + msh.Omg{ip}.nodes; 2 * msh.nnode + msh.Omg{ip}.nodes];

    msh.Omg{ip}.processor_boundary = setdiff(msh.Omg{ip}.dofbar,msh.Omg{ip}.dof); % Degrees of freedom which live on processor boundary

    msh.Omg{ip}.tdof = length(msh.Omg{ip}.dofbar);

        % Need to define local numbering on processor to do assembly

    % Initialise connectivity matrix on this partition (with nodes in local number)

    msh.Omg{ip}.connectivity = zeros(msh.Omg{ip}.nelem,msh.enode);

    % This segment of code is constructuing connectivity matrix
    % for a partition in local partion node and element numbering
    
    tmp = zeros(msh.Omg{ip}.nnode,1);

    for i = 1 : msh.Omg{ip}.nnode % for each node in overlapping partition
        
        nodeNum = msh.Omg{ip}.nodes(i); % global node number
        
        elems = msh.nodes{nodeNum}.elements; % Elements which have nodeNum belongs to (some of the elements might not be in this partition)
        
        [elems,~,elem_ids_local] = intersect(elems,msh.Omg{ip}.elems); % Ensure elements belong to partition
        
        % For each of these, nodeNum should appear in their connectivity
        % list as i (local numbering)
        
        for j = 1 : length(elems)       
            node_position = find(msh.e2g(elems(j),1:msh.enode) == nodeNum); % Position in global connectivity matrix      
            msh.Omg{ip}.connectivity(elem_ids_local(j),node_position) = i;    
        end
        
        
    end % end for each node
        
        %[~,ID,~] = intersect(msh.Omg{ip}.dofbar,msh.Omg{ip}.dof); 
        
        R{ip} = sparse(1:length(msh.Omg{ip}.dofbar),msh.Omg{ip}.dofbar,ones(length(msh.Omg{ip}.dofbar),1),length(msh.Omg{ip}.dofbar),msh.tdof);
        
        for ie = 1 : msh.Omg{ip}.nelem
            
            msh.Omg{ip}.elements{ie}.connectivity = msh.Omg{ip}.connectivity(ie,:);   
            
        end
        
        
        
        % === Compute element to global for each partition
        
        msh.Omg{ip}.e2g = zeros(msh.Omg{ip}.nelem,msh.nedof);  
        msh.Omg{ip}.e2g(:,1:msh.enode) = msh.Omg{ip}.connectivity;
        msh.Omg{ip}.e2g(:,msh.enode + 1:2 * msh.enode) = msh.Omg{ip}.nnode + msh.Omg{ip}.connectivity;
        msh.Omg{ip}.e2g(:,2 * msh.enode + 1:msh.nedof) = 2 * msh.Omg{ip}.nnode + msh.Omg{ip}.connectivity;
        
        
        % === Global Dirichlet Constraints
        
        [bnd_nodes,~,id] = intersect(union(msh.lhs,msh.rhs),msh.Omg{ip}.nodes);
        
        msh.Omg{ip}.BND = [id; id + msh.Omg{ip}.nnode; id + 2 * msh.Omg{ip}.nnode]; % Global Dirichlet Degrees of Freedom in local numbering
        
        msh.Omg{ip}.FREE = (1 : length(msh.Omg{ip}.dofbar))'; % FREE degrees of freedom in local numbering
        
        msh.Omg{ip}.FREE(msh.Omg{ip}.BND) = []; % Remove BND
                                                % degrees of
                                                % freedom

end

    % ================================================================================
    %                         Construct Xi
    % ================================================================================
for ip = 1 : numParts
    
    msh.Omg{ip}.xi = zeros(msh.Omg{ip}.nnode,1);

    for j = 1 : msh.Omg{ip}.nnode % for each node in the partition
        
        
        x = msh.Omg{ip}.coords(j,:);
        
        
        flag = 0;
        
        for ie = mshC.nodes{ip}.elements
            if(inhull(x,mshC.coords(mshC.elements{ie}.connectivity, ...
                                    :),mshC.elements{ie}.tess,1e-5))
                flag = 1;
                break
            end                
        end   
            
        if (flag > 0)
            % Local number of ip in macro element
            macroLocNum = find(ip == mshC.elements{ie}.connectivity);
            
            msh.Omg{ip}.xi(j) = computeXi(x, ...
                                          mshC.elements{ie}.invX,macroLocNum);
        else
            msh.Omg{ip}.xi(j) = 0.0;
        end
        
    end
    
    % ================================================================================
        %                         Plot Xi
        % ================================================================================
        if (plotxi)
                scalar.name = 'Shape_functions';
                scalar.data = msh.Omg{ip}.xi; 

                output_file = strcat('shapefunctions_',num2str(ip),'.vtk');
                % need to get global node numbers for plotting to work
                
                matlab2vtk (output_file, 'solution', msh.Omg{ip}, 'hex' , scalar, [], [])
        end
    
        
    msh.Omg{ip}.xi = repmat(msh.Omg{ip}.xi,3,1); % repeat because 3
                                                 % dofs per node
                                                 
    msh.Omg{ip}.xi = spdiags(msh.Omg{ip}.xi, 0,length(msh.Omg{ip}.xi),length(msh.Omg{ip}.xi));

end


%%%%%% Construct Neumann matrix on each partition

K = cell(numParts,1);

M = cell(numParts,1);

for ip = 1 : numParts
    
    indx_j = repmat(1:msh.nedof,msh.nedof,1); indx_i = indx_j';
    Kindx.i = msh.Omg{ip}.e2g(:,indx_i(:)); Kindx.j = msh.Omg{ip}.e2g(:,indx_j(:));
    
    Ke_all = zeros(msh.nedof^2,msh.Omg{ip}.nelem);
    Me_all = zeros(msh.nedof^2,msh.Omg{ip}.nelem);
    
    for j = 1 : msh.Omg{ip}.nelem % for each element in partition
        
        ie = msh.Omg{ip}.elems(j);       
        
        % Constuct element stiffness matrix
        
        region = msh.elements{ie}.region;
        
        if ss(region) >= 0
            ElasticTensor = composite;
            ang = ss(region);
        else
            ElasticTensor = isotropic;
            ang = 0;
        end
        
        Ke_all(:,j) = elementStiffness(msh.coords(msh.elements{ie}.connectivity,:),msh.elements{ie},msh.dN,msh.nip,ElasticTensor,ang,msh.ip);
        
        Me_all(:,j) = elementMassStiffness(msh.coords(msh.elements{ie}.connectivity,:),msh.elements{ie},msh.N,msh.dN,msh.nip,ElasticTensor,ang,msh.ip);
        
    end
    
    K{ip} = sparse(Kindx.i',Kindx.j',Ke_all);
    M{ip} = sparse(Kindx.i',Kindx.j',Me_all);
  
    
end

% === Compute Generalised EigenValue Problem for each overlapping subdomain
tic
V = cell(numParts,1); % Cell container for set of eigenvalues on each subdomain
lambda = cell(numParts,1);

which_partition = [];
counter = 1;

mylhs = zeros(numParts,1);
myrhs = zeros(numParts,1);
for ip = 1 : numParts
    
    % Record if function macro dof lives on left and right boundary
    
    if (abs(mshC.coords(ip,1) - L) < 1.0e-6)
        myrhs(ip) = 1;
    end
    
    if (abs(mshC.coords(ip,1)) < 1.0e-6)
        mylhs(ip) = 1;
    end
    
    
    V{ip}.ALL = zeros(length(K{ip}),numEigs + 3);
    
    % First of all add zero energy modes to each partition (u,v,w) 
    
    cnst_function = ones(msh.Omg{ip}.nnode,1);
    my_zeros = zeros(msh.Omg{ip}.nnode,1);
    
    V{ip}.ALL(:,1) = vertcat(cnst_function,my_zeros,my_zeros);
    V{ip}.ALL(:,2) = vertcat(my_zeros,cnst_function,my_zeros);
    V{ip}.ALL(:,3) = vertcat(my_zeros,my_zeros,cnst_function);
    
   
%      Rx = [ 1 , 0, 0; 0,  0, -1; 0, 1, 0];
%      Ry = [ 0 , 0, 1; 0,  1, 0; -1, 0, 0];
%      Rz = [ 0, -1, 0; 1,  0, 0; 0, 0, 1];
%     
%     for i = 1 : msh.Omg{ip}.nnode
%         
%         x = msh.Omg{ip}.coords(i,:);
%         
%         V{ip}.ALL([i,i+msh.Omg{ip}.nnode,i + 2*msh.Omg{ip}.nnode],4) = Rx * x' - x';
%         V{ip}.ALL([i,i+msh.Omg{ip}.nnode,i + 2*msh.Omg{ip}.nnode],5) = Ry * x' - x';
%         V{ip}.ALL([i,i+msh.Omg{ip}.nnode,i + 2*msh.Omg{ip}.nnode],6) = Rz * x' - x'; 
%         
%     end
    
    % Compute non-zero energy modes
    
    if (numEigs > 0)
        
        % Check if any of the boundary is on global boundary
        
        elem_in_ip = mshC.nodes{ip}.elements;
        
        nodes_in_ip = unique(mshC.e2g(elem_in_ip,1:8));
        flag = false;
        computedEigs = numEigs;
        if(isempty(intersect(nodes_in_ip,union(mshC.lhs,mshC.rhs))))
            computedEigs = numEigs + 6;
            flag = true;
        end

        mshlocal = msh.Omg{ip};
        A = K{ip}(msh.Omg{ip}.FREE,msh.Omg{ip}.FREE);

        if (strcmp(BMatrix,'GenEO'))
            B = msh.Omg{ip}.xi * K{ip} * msh.Omg{ip}.xi;
        else
            B = M{ip};
        end

        B = B(msh.Omg{ip}.FREE,msh.Omg{ip}.FREE);

        [tmpEV,D] = eigs(A,B,computedEigs,1e-6); % Solve Generalised
                                                 % EigenValue Problem  
                                          
        if(flag)                                       
            [~,id] = sort(diag(D));
            id(1:6) = [];
        else
            id = 1 : computedEigs;
        end
                                                 
        %pos_eigvals = find(diag(D) > 1e-8); % First 6 will be zero energy modes if not constrainted by global dirichlet constraint
        
        %num_pos_eigs(ip) = length(pos_eigvals);

        V{ip}.ALL(msh.Omg{ip}.FREE,4: numEigs+3) = tmpEV(:,id);
    
    end
    
%     if(plot_eigenmodes)
%         
%         vector.name = 'displacements';
%         vector.data = zeros(msh.nnode,3);
% 
%         vector.data(:,1) = u(1:msh.nnode);
%         vector.data(:,2) = u(msh.nnode + 1: 2 * msh.nnode);
%         vector.data(:,3) = u(2 * msh.nnode + 1: 3 * msh.nnode);
% 
%         output_file = strcat('solution_composite_3D.vtk');
% 
%         matlab2vtk (output_file, 'solution', msh, 'hex' , [], vector, [])
%         
%         
%     end
   
    ip
end

% Now we construct the coarse space and do an inexpensive solve
% up to this point the code only needs to run once!

% construct R_H operator

total_eigs = 0;
for ip = 1 : numParts
    eigs_on_ip(ip) = length(V{ip}.ALL(1,:));
end
total_eigs = sum(eigs_on_ip);

KH = sparse(msh.tdof,msh.tdof);
RH = sparse(msh.tdof,total_eigs);
running_sum = 0;
for ip = 1 : numParts
    RH(:,(1:eigs_on_ip(ip)) + running_sum) = R{ip}'* msh.Omg{ip}.xi * ...
        V{ip}.ALL; % global dirichlet bnd conds included in V
    running_sum = running_sum + eigs_on_ip(ip);
    
    % let us assemble the full K for now
    indk = msh.Omg{ip}.dofbar;
    KH(indk,indk) = KH(indk,indk) + (K{ip} * msh.Omg{ip}.xi);
   %KH = KH + R{ip}' * msh.Omg{ip}.xi * K{ip} * ;
end


% if (test_code)
%     % KH_Test 
% 
%     indx_j = repmat(1:msh.nedof,msh.nedof,1); indx_i = indx_j';
%     Kindx.i = msh.Omg{ip}.e2g(:,indx_i(:)); Kindx.j = msh.Omg{ip}.e2g(:,indx_j(:));
% 
%     Ke_all = zeros(msh.nedof^2,msh.Omg{ip}.nelem);
% 
%         for j = 1 : msh.nelem % for each element in partition      
% 
%             % Constuct element stiffness matrix
% 
%             region = msh.elements{ie}.region;
% 
%             if ss(region) >= 0
%                 ElasticTensor = composite;
%                 ang = ss(region);
%             else
%                 ElasticTensor = isotropic;
%                 ang = 0;
%             end
% 
%             Ke_all(:,j) = elementStiffness(msh.coords(msh.elements{ie}.connectivity,:),msh.elements{ie},msh.dN,msh.nip,ElasticTensor,ang,msh.ip);
% 
%         end
% 
%     KHT = sparse(Kindx.i',Kindx.j',Ke_all);
% 
%     assert(sum(sum(KHT - KH)) < 1e-10,'Matrix Assembly of Full K is Wrong - likely bug in formulation of xi');
% end

% Bug is not xi

% Bug must be in Boundary conditions or in eigenvectors

opts.isreal = 1;
opts.issym = 1;

KH = RH' * KH * RH; % KH is now nM x nM

KH = 0.5 * (KH + KH');



 rbid = [1 + (which_partition - 1)*numEigs];

% bug is in the boundary conitions

%assert(abs(eigs(KH,1,'SM'))>1e-12,'KH is singular - before boundary conditions are applied');

cs = cumsum(eigs_on_ip);

FREE = 1 : total_eigs;
BND = zeros(total_eigs,1);
U = zeros(total_eigs,1);

myrhs = find(myrhs);
mylhs = find(mylhs);

for i = 1 : length(myrhs)
    
    j = myrhs(i);
    
    if j > 1
        start_point = cs(j-1);
    else
        start_point = 0;
    end  
    BND(start_point + (1:3)) = 1;
    U(start_point + 3) = 1;
end

for i = 1 : length(mylhs) 
    j = mylhs(i);
    if j > 1
        start_point = cs(j-1);
    else
        start_point = 0;
    end
    BND(start_point + (1:3)) = 1;  
end

BND = logical(BND);

FREE(BND) = [];

U(FREE) = KH(FREE,FREE) \ (-KH(FREE,BND) * U(BND));

%eigs(KH(FREE,FREE),10,'SM')



% for wp = 1: length(which_partition)
%     rbid(wp) = (which_partition(wp) - 1)*numEigs + wp; % because
%                                                        % we're
%                                                        % adding 1
%                                                        % rigid body
%                                                        % mode per
%                                                        % partition
%                                                        % in which_partition
% end
% KH(rbid,:) = 0;
% for eg = 1 : length(rbid)
%     KH(rbid(eg),rbid(eg)) = 1.0;
% end
% บบ
% %assert(abs(eigs(KH,1,'SM'))>0.0,'KH is singular - after boundary conditions are applied');
% 
% 
% F = zeros(length(KH),1);
% F(rbid) = 1.0;
% 
% % coarse solve
% U = KH \ F;

u = RH * U; % transform macro scale solution back to fine scale
solve_time = toc

%[ker,D] = eigs(KH,1,'SM');

%u = RH * ker;

vector.name = 'displacements';
vector.data = zeros(msh.nnode,3);

vector.data(:,1) = u(1:msh.nnode);
vector.data(:,2) = u(msh.nnode + 1: 2 * msh.nnode);
vector.data(:,3) = u(2 * msh.nnode + 1: 3 * msh.nnode);

output_file = strcat('solution_composite_3D.vtk');

matlab2vtk (output_file, 'solution', msh, 'hex' , [], vector, [])
