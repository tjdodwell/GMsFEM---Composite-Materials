


function xi = computeXi(x,invX,j)

    phi= zeros(8,1); phi(j) = 1;

    C = invX * phi;

    xi = [1, x(1), x(2), x(3), x(1)*x(2), x(1)*x(3), x(2)*x(3), ...
          x(1)*x(2)*x(3)] * C;

end

