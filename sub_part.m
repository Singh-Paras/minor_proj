function [locs, pks] = sub_part(X,band)
for i = 1:band-1
    V(i) = X(i,i+1);
end

%     figure();
%     findpeaks(V);
[pks,locs] = findpeaks(V);
%     figure();

V_S = smoothdata(V,'gaussian',18);

%     plot(1:band-1, V_S);
%     figure();
%     findpeaks(V_S);
[pks,locs] = findpeaks(V_S);
end