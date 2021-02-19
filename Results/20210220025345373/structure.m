function [u,a] = structure(n_layer,period,N,Nm,k0,diameter_x,diameter_y,Ntre,init,op_retcouche)

u=[];
a=[];
for az=1:n_layer
    if Nm(az)==0
        u{az}=retu(period,{N(az),k0});
    else
        dx=11; dy=11;
        if N(az)==1
            structure_array = {};
            for px=-(dx-1)/2:(dx-1)/2
                for py=-(dy-1)/2:(dy-1)/2
                    structure_array = [structure_array, [period(1)*1000*px/dx,period(2)*1000*py/dy,diameter_x(az),diameter_y(az),Nm(az),Ntre]];
%                    structure_array = [structure_array, [period(1)*1000*px/dx,period(2)*1000*py/dy,period(1)/dx*0.9,0.9*period(2)/dy,Nm(az),Ntre]];
                end
            end
            texture=[{N(az)},structure_array(:)',{k0}];
        else
            structure_layer=[0,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
            texture={N(az),structure_layer,k0};
        end

        u{az}=retu(period,texture);
    end
    a{az}=retcouche(init,u{az},op_retcouche);
end