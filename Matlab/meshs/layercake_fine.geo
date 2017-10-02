elem_size = 1;

numLayers = 7;

Nodes_X = 10 * 5;
Nodes_Y = 10;
Nodes_Z = 4;

Lx = 5.0;
Ly = 1.0;

th = {0.2,0.02,0.2,0.02,0.2,0.02,0.2};

p = 0;
i = 0;
vv = 1;

For ( 1 : numLayers )

	If (i == 0)
		// For First Layer Only

		Point(1) = {0,0,0,elem_size};
		Point(2) = {Lx,0,0,elem_size};
		Point(3) = {Lx,Ly,0,elem_size};
		Point(4) = {0,Ly,0,elem_size};

		Line(1) = {1,2};
		Line(2) = {2,3};
		Line(3) = {3,4};
		Line(4) = {4,1};

		Line Loop(1) = {1,2,3,4};
		Plane Surface(1) = {1};

		Transfinite Line{1,3} = Nodes_X;
		Transfinite Line{2,4} = Nodes_Y;

		Transfinite Surface {1} = {1:4};

		topLastLayer = 0;
		p = 4;
		l = 4;
		ll = 1;

	EndIf

	Point(p+1) = {0,0,topLastLayer+th[i],elem_size};
	Point(p+2) = {Lx,0,topLastLayer+th[i],elem_size};
	Point(p+3) = {Lx,Ly,topLastLayer+th[i],elem_size};
	Point(p+4) = {0,Ly,topLastLayer+th[i],elem_size};

	Line(l+1) = {p-3,p+1};
	Line(l+2) = {p-2,p+2};
	Line(l+3) = {p-1,p+3};
	Line(l+4) = {p,p+4};

	Line(l+5) = {p+1,p+2};
	Line(l+6) = {p+2,p+3};
	Line(l+7) = {p+3,p+4};
	Line(l+8) = {p+4,p+1};

	Line Loop (ll+1) = {l-3,l+2,-(l+5),-(l+1)};
	Line Loop (ll+2) = {l-2,l+3,-(l+6),-(l+2)};
	Line Loop (ll+3) = {l-1,l+4,-(l+7),-(l+3)};
	Line Loop (ll+4) = {l  ,l+1,-(l+8),-(l+4)};
	Line Loop(ll+5) = {l+5,l+6,l+7,l+8};

	Plane Surface(ll+1) = {ll+1};
	Plane Surface(ll+2) = {ll+2};
	Plane Surface(ll+3) = {ll+3};
	Plane Surface(ll+4) = {ll+4};
	Plane Surface(ll+5) = {ll+5};

	Surface Loop(vv) = {ll:ll+5};

	Volume(vv) = {vv};

	Physical Volume(i) = {vv};

	Transfinite Line{l+5,l+7} = Nodes_X;
	Transfinite Line{l+6,l+8} = Nodes_Y;

	Transfinite Line{l+1:l+4} = Nodes_Z;

	Transfinite Surface(ll+1) = {p-3,p-2,p+2,p+1};
	Transfinite Surface(ll+2) = {p-2,p-1,p+3,p+2};
	Transfinite Surface(ll+3) = {p-1,p,p+4,p+3};
	Transfinite Surface(ll+4) = {p,p-3,p+4,p+1};
	Transfinite Surface(ll+5) = {p+1,p+2,p+3,p+4};

	Transfinite Volume{vv} = {p-3:p+4};

	Recombine Surface{ll:ll+5};

	p = p + 4;
	l = l + 8;
	ll = ll + 5;
	vv = vv + 1;


	topLastLayer = topLastLayer + th[i];
	i = i + 1;

EndFor