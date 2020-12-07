function [u,a] = structure(Nb_couches,period,N,Nm,k0,diameter_x,diameter_y,Ntre,init,op_retcouche)

u=[];
a=[];
for az=1:Nb_couches
    if Nm(az)==0
        u{az}=retu(period,{N(az),k0});
    else
        % u{az}=retu(period,{N(az),[0,0,diameter_x,diameter_y,Nm(az),Ntre],[-diameter_x/2+w_rectangle/2,diameter_y/2+h_rectangle/2,w_rectangle,h_rectangle,Nm(az),Ntre],[diameter_x/2+h_rectangle/2,-diameter_y/2+w_rectangle/2,h_rectangle,w_rectangle,Nm(az),Ntre ],k0});
        if N(az)==1    %Nm(az)~=0
            structure_array = {};
            for px=-5:5
                for py=-5:5
                    structure_array = [structure_array, [1200*2*px/11,1200*2*py/11,diameter_x(az),diameter_y(az),Nm(az),Ntre]];
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