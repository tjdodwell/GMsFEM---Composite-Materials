function msh = defineIPs(msh)

    msh.nip = 8;
    
    P   =   0.577350269189626;
    
    msh.ip.coords(1,:) = [-P,-P,-P];
    msh.ip.coords(2,:) = [P,-P,-P];
    msh.ip.coords(3,:) = [P,P,-P];
    msh.ip.coords(4,:) = [-P,P,-P];
    
    msh.ip.coords(5,:) = [-P,-P,P];
    msh.ip.coords(6,:) = [P,-P,P];
    msh.ip.coords(7,:) = [P,P,P];
    msh.ip.coords(8,:) = [-P,P,P];

    msh.ip.wgts(1:8) = ones(8,1);
    
    % cohesive element ips
    msh.coh.nip = 4;
    
    msh.coh.ip.coords(1,:) = [-P,-P,0];
    msh.coh.ip.coords(2,:) = [P,-P,0];
    msh.coh.ip.coords(3,:) = [P,P,0];
    msh.coh.ip.coords(4,:) = [-P,P,0];
    
    msh.coh.wgts(1:4) = ones(4,1);

end
