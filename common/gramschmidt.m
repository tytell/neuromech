function V = gramschmidt(U)

k = size(U,2);
V = U;
for i = 1:k
    V(:,i) = V(:,i) / sqrt(sum(V(:,i).^2));
    for j = i+1:k
        proj = sum(V(:,i) .* V(:,j)) .* V(:,i);
        V(:,j) = V(:,j) - proj;
    end
end

    