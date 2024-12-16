function L2 = ComputeLinearOperator_1D_Fisher(nx, h)
    % Neumann 2nd derivative matrix, 2nd order
    e = ones(nx, 1);
    L2 = spdiags([e -2*e e], -1:1, nx, nx) / h^2;
    L2(1, 2) = 2 / h^2;       % Left Neumann BC
    L2(end, end-1) = 2 / h^2; % Right Neumann BC
end