function Me = elementMassStiffness(coords,elem,N,dN,nip,C,ang,intpoints)

    Msub = zeros(8);
    
    Me = zeros(24);

	for ip = 1 : nip

	  	J	=	coords'*dN{ip};
	    dNdX	=	dN{ip}*inv(J); %#ok<MINV>
        
        Msub = Msub + N{ip} * N{ip}' * det(J) * intpoints.wgts(ip);

	end	% for each ip

    Msub = 0.5*(Msub + Msub');
    
    Me = blkdiag(Msub,Msub,Msub);
  
    
    Me = Me(:); % Convert to a column vector
    
end