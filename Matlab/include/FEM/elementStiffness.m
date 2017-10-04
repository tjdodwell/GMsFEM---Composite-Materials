function Ke = elementStiffness(coords,elem,dN,nip,C,ang,intpoints)

    Ke = zeros(24);

	for ip = 1 : nip

            J	=	coords'*dN{ip};
	    dNdX	=	dN{ip}*inv(J); %#ok<MINV>

	    % For Composite Elements

		C = Rotation_Matrices(C,[0,0,ang]);

		% Local Coordinates to Global Coordinates

		th = zeros(3,1);

		C = Rotation_Matrices(C,th);
        
        % Calculate Strain Matrix
        
        B = zeros(6,24);

		B(1,1:8) = dNdX(:,1); % e_11 = u_1,1
		B(2,9:16) = dNdX(:,2); % e_22 = u_2,2 
		B(3,17:24) = dNdX(:,3); % e_33 = u_3,3

		B(4,9:16) = dNdX(:,3);	B(4,17:24) = dNdX(:,2);	% e_23 = u_2,3 + u_3,2
		B(5,1:8) = dNdX(:,3);	B(5,17:24) = dNdX(:,1);	% e_13 = u_1,3 + u_3,1
		B(6,1:8) = dNdX(:,2);	B(6,9:16) = dNdX(:,1);	% e_12 = u_1,2 + u_2,1
        
	    Ke = Ke + B'*C*B*intpoints.wgts(ip)*det(J);

	end	% for each ip

    Ke = 0.5*(Ke + Ke');
    
    if sum(sum(abs(Ke-Ke'))) > 1e-12
        error('Element Stiffness Matrix Not Symmetric');
    end
    
    Ke = Ke(:); % Convert to a column vector
