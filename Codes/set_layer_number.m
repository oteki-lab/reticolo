function y=set_layer_number(n,zou,Nb_couches)

Number = zeros(length(n(:,1)));
for index = 1:length(n(:,1))
    if length(n(index,:))==1; Number(index)=n(index); else; Number(index)=n(index, zou); end
end

y=Number(1:Nb_couches);
