function [C] = sufdpc(dist,ent,pages)

    ND = size(dist, 1);
    % Compute dc value
%     percent=2.0/exp(ND/pages);
    percent=2.0;

    position=round(ND*(ND-1)/2*percent/100);
    if position == 0
        position = 1;
    end


    temp = dist(find(tril(dist)~=0));
    sda = sort(temp);
    %%%%%
    dini=sda(position);
    k=10;
    dc = dini/exp(ND/pages);
    % Compute the rho factor
    rho = zeros(ND,1);
    for i = 1:ND
        for j = 1:ND
            if i~=j
                rho(i) =rho(i)+ exp(-(dist(i,j)/dc)^2);
            end
        end
    end
    % Compute the delta factor
    [rho_sorted,ordrho]=sort(rho,'descend');
    maxd=max(max(dist(ordrho(1),:)));
    delta = zeros(ND,1);
    delta(ordrho(1))=-1;
    maxD = max(max(dist));
    for i = 2:ND
        delta(ordrho(i))=maxD;
        for j=1:i-1
            if dist(ordrho(i),ordrho(j))<delta(ordrho(i))
                delta(ordrho(i)) = dist(ordrho(i),ordrho(j));
            end
        end
    end
    delta(ordrho(1)) = max(delta);
    % normalize the factors
    rho = (rho-min(rho(:)))/(max(rho(:))-min(rho(:)));
    delta = (delta-min(delta(:)))/(max(delta(:))-min(delta(:)));
    
    % The final importance
    gamma =rho.*delta.*ent';
    
    % Find the selected bands
    [~,order_band] = sort(gamma,'descend');
%     C = order_band(1:k);
    C = order_band(1:ceil(ND/pages*k));
    
end
