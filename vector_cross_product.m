function [dir, mag] = vector_cross_product(v1, v2)
    % CROSS Vector Cross Product
    %
    % [dir, mag] = vector_cross_product(v1, v2)
    %
    % INPUT ARGUMENTS
    % ================
    % v1 vector 1 : [x y z]
    % v2 vector 2 : [x y z]
    %
    % OUTPUT ARGUMENTS
    % ================
    % dir Direction of the output vector
    % .x component in x-direction
    % .y component in y-direction
    % .z component in z-direction
    % mag Magnitude of the output vector
    dir.x = v1(2)*v2(3)-v1(3)*v2(2);
    dir.y = v1(3)*v2(1)-v1(1)*v2(3);
    dir.z = v1(1)*v2(2)-v1(2)*v2(1);
    mag = sqrt(dir.x^2+dir.y^2+dir.z^2);
end