function y = LineCut(I,x1,x2)
    
    x1 = round(x1);
    x2 = round(x2);
    temp = [x1,x2]';
    x1 = temp(:,1);
    x2 = temp(:,2);
    [n,dim] = max(abs(x2 - x1));
    k = (x2(2)-x1(2))/(x2(1)-x1(1));
    y = zeros(n,1);
    if dim == 1
        for ii = 1:n
            loc = round([x1(1)+ii;x1(2)+k*ii]);
            loc = sub2ind(size(I),loc(2),loc(1));
            y(ii) = mean(I([loc,loc+1,loc-1]));
            I(loc) = 0;
        end
    else
        for ii = 1:n
            loc = round([x1(1)+ii/k;x1(1)+ii]);
            loc = sub2ind(size(I),loc(2),loc(1));
            y(ii) = mean(I([loc,loc+1,loc-1]));
            I(loc) = 0;
        end
    end
    imagesc(I)
    axis image
    axis off
    figure
    plot(y,'linewidth',1.5);
    xlabel('Pixel');
    ylabel('Intensity');
    title('剖面强度');
    set(gca,'fontsize',15,'fontweight','bold');
    
end