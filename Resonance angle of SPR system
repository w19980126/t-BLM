function R = MultiLayerReflect(lambda,theta0,n,d)
    
% The input parameters include wavelength(lambda), incident Angle(theta0), 
% refractive index(n) and thickness(d) of each layer
% The unit of wavelength and thickness is nm
% The output parameter R is the reflection of the multilayer sysmtem.

    d = d/10^9; 
    lambda = lambda/10^9;
    theta = zeros(length(n),1);
    theta(1) = theta0;
    for ii = 2:length(n)
        theta(ii) = asin(n(1)/n(ii)*sin(theta(1)));   
        % Calculating the angle of reflection of each layer by Fresnell's law.
        % Note: Use the incident Angle of the first layer, otherwise it is easy 
        % to oscillate due to data overflow during calculation.
    end
    M = [1 0;0 1];
    temp_M = M;
    
    for ii = 2:length(n)-1
        delta = 2*pi/lambda*n(ii)*d(ii)*cos(theta(ii));
        eta = n(ii)/cos(theta(ii));
        temp_M(1) = cos(delta);
        temp_M(2) = -1i*eta*sin(delta);
        temp_M(3) = -1i/eta*sin(delta);
        temp_M(4) = temp_M(1);
        M = M*temp_M;
    end

    eta_in = n(1)/cos(theta(1));
    eta_out = n(end)/cos(theta(end));
    M = M*[1;eta_out];
    Y = M(2)/M(1);
    r = (eta_in - Y)/(eta_in + Y);
    R = (abs(r))^2;
    
end
