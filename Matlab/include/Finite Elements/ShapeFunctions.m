function [N,dN] = ShapeFunctions(msh)


	N = cell(msh.nip,1); dN = cell(msh.nip,1);

	for ip = 1:msh.nip

		N{ip} = zeros(8,1); dN{ip} = zeros(8,3);

		xi = msh.ip.coords(ip,1);	eta = msh.ip.coords(ip,2);	mu = msh.ip.coords(ip,3);

		% Displacement Shape Functions
		
		N{ip}(1) = 0.125 * (1 - xi) * (1 - eta) * (1 - mu);
		N{ip}(2) = 0.125 * (1 + xi) * (1 - eta) * (1 - mu);
		N{ip}(3) = 0.125 * (1 + xi) * (1 + eta) * (1 - mu);
		N{ip}(4) = 0.125 * (1 - xi) * (1 + eta) * (1 - mu);
		N{ip}(5) = 0.125 * (1 - xi) * (1 - eta) * (1 + mu);
		N{ip}(6) = 0.125 * (1 + xi) * (1 - eta) * (1 + mu);
		N{ip}(7) = 0.125 * (1 + xi) * (1 + eta) * (1 + mu);
		N{ip}(8) = 0.125 * (1 - xi) * (1 + eta) * (1 + mu);

		% Derivative Shape Functions

		dN{ip}(1,1) = -0.125 * (1 - eta) * (1 - mu);
		dN{ip}(2,1) =  0.125 * (1 - eta) * (1 - mu);
		dN{ip}(3,1) =  0.125 * (1 + eta) * (1 - mu);
		dN{ip}(4,1) = -0.125 * (1 + eta) * (1 - mu);
		dN{ip}(5,1) = -0.125 * (1 - eta) * (1 + mu);
		dN{ip}(6,1) =  0.125 * (1 - eta) * (1 + mu);
		dN{ip}(7,1) =  0.125 * (1 + eta) * (1 + mu);
		dN{ip}(8,1) = -0.125 * (1 + eta) * (1 + mu);


		dN{ip}(1,2) = -0.125 * (1 - xi) * (1 - mu);
		dN{ip}(2,2) = -0.125 * (1 + xi) * (1 - mu);
		dN{ip}(3,2) =  0.125 * (1 + xi) * (1 - mu);
		dN{ip}(4,2) =  0.125 * (1 - xi) * (1 - mu);
		dN{ip}(5,2) = -0.125 * (1 - xi) * (1 + mu);
		dN{ip}(6,2) = -0.125 * (1 + xi) * (1 + mu);
		dN{ip}(7,2) =  0.125 * (1 + xi) * (1 + mu);
		dN{ip}(8,2) =  0.125 * (1 - xi) * (1 + mu);


		dN{ip}(1,3) = -0.125 * (1 - xi) * (1 - eta);
		dN{ip}(2,3) = -0.125 * (1 + xi) * (1 - eta);
		dN{ip}(3,3) = -0.125 * (1 + xi) * (1 + eta);
		dN{ip}(4,3) = -0.125 * (1 - xi) * (1 + eta);
		dN{ip}(5,3) =  0.125 * (1 - xi) * (1 - eta);
		dN{ip}(6,3) =  0.125 * (1 + xi) * (1 - eta);
		dN{ip}(7,3) =  0.125 * (1 + xi) * (1 + eta);
		dN{ip}(8,3) =  0.125 * (1 - xi) * (1 + eta);

	end	%	end	N and dN

end	% 3DShapeFunctions