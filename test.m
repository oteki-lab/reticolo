function Z=test(X, Y)
%[X,Y] = meshgrid(0:.5:32);
%Z1 = exp(-((X-9).^2+(Y-9).^2)/(2*2.^2));
%Z2 = 2*exp(-((X-2).^2+(Y-20).^2)/(2*10.^2));
%Z3 = 2*exp(-((X-25).^2+(Y-5).^2)/(2*10.^2));
%Z = Z1 + Z2 + Z3;
%figure
%mesh(X,Y,Z)

if X > Y
    T = readtable('DB_map_max.csv');
    rows = (T.Eg==X & T.Eic==Y);
    Z = T(rows,:).('eff');
else
    Z= 0.0;
end