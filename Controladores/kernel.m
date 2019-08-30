function [k] = kernel(sigma_quad, x1, x2, l, n)
    for i=1:n
            k(:,i) = sigma_quad*exp(-((x1 - x2(:,i)).^2)/(2*(l^2)));
    end
    
    k = k';
end