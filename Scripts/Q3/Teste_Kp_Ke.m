S = load('L04_1v_versao2.mat'); sigs = S.input.signals.values;
u = sigs(:,1);  alphae = sigs(:,3);
i = find(diff(u) > 1);                     % subidas do comando (-1 V -> +1 V)
i = i(i > 20 & i < length(u) - 20);        % evita as pontas
d = zeros(size(i));
for k = 1:numel(i)
    d(k) = mean(alphae(i(k)+1 : i(k)+10)) - mean(alphae(i(k)-10 : i(k)));
end
fprintf('Variacao media de alphae apos subida do comando: %+.3f V\n', mean(d));