function [ d ] = func_distance( M, K )
%func_distance ����M��AP��K���û����Եľ���
s = [size(M,1),size(K,1)];
d = zeros(s);
for i = s(1):-1:1
    d(i,:) = sqrt( (M(i,1)-K(:,1)).^2 + (M(i,2)-K(:,2)).^2 );
end

