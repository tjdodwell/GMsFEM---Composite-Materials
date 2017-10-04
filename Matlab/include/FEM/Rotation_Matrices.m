function C = Rotation_Matrices(C,th)

R = cell(3,1);

% Rotation about x_1

c = cos(th(1)); s = sin(th(1));

R{1} = zeros(6);

R{1}(1,1) = 1;
R{1}(2,2) = c * c; 
R{1}(2,3) = s * s;	
R{1}(2,4) = 2 * c * s;
R{1}(3,2) = s * s;
R{1}(3,3) = c * c;
R{1}(3,4) = -2 * c * s;
R{1}(4,2) = -c * s;
R{1}(4,3) = c * s;
R{1}(4,4) = (c * c) - (s * s);
R{1}(5,5) = c;
R{1}(5,6) = -s;
R{1}(6,5) = s;
R{1}(6,6) = c;

% Rotation about x_2

c = cos(th(2)); s = sin(th(2));

R{2} = zeros(6);

R{2}(1,1) = c * c;
R{2}(1,3) = s * s;
R{2}(1,5) = 2 * c * s;
R{2}(2,2) = 1;
R{2}(3,1) = s * s;
R{2}(3,3) = c * c;
R{2}(3,5) = -2 * c * s;
R{2}(4,4) = c;
R{2}(4,6) = -s;
R{2}(5,1) = -c * s;
R{2}(5,3) = c * s;
R{2}(5,5) = (c * c) - (s * s);
R{2}(6,4) = s;
R{2}(6,6) = c;

% Rotation about x_3

c = cos(th(3)); s = sin(th(3));

R{3} = zeros(6);

R{3}(1,1) = c * c;
R{3}(1,2) = s * s;
R{3}(1,6) = 2 * c * s;
R{3}(2,1) = s * s;
R{3}(2,2) = c * c;
R{3}(2,6) = -2 * c * s;
R{3}(3,3) = 1;
R{3}(4,4) = c;
R{3}(4,5) = s;
R{3}(5,4) = -s;
R{3}(5,5) = c;
R{3}(6,1) = -c * s;
R{3}(6,2) = c * s;
R{3}(6,6) = (c * c) - (s * s);

T = R{3}*R{2}*R{1};

C = T*C*T';


end