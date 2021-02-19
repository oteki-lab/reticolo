function filename = calc_cross_section(in, x0, y0, z0, init,struct0,sh,sb,inc, tab0, xx, yy, zou)
    if isempty(x0)==1&&isempty(z0)==1
        cs=y0;
        [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{xx,y0},tab0,[],(1:6)+7.25i,1,1,1:6);
    elseif isempty(y0)==1&&isempty(z0)==1
        cs=x0;
        [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{x0,yy},tab0,[],(1:6)+7.25i,1,1,1:6);
    elseif isempty(x0)==1&&isempty(y0)==1
        cs=z0;
        tab0(:,3)=0;
        HH=cumsum(tab0(:,1));
        numz=1;
        while (HH(end)-z0)>HH(numz);numz=numz+1; end
        tab1=[tab0(1:numz-1,:);[tab0(numz,1)-(z0-sum(tab0(numz+1:end,1))),numz,0];[0,numz,1];[z0-sum(tab0(numz+1:end,1)),numz,0];tab0(numz+1:end,:)];
        [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{xx,yy},tab1,[],(1:6)+7.25i,1,1,1:6);
    end

    for ii=1:3
        o0(:,:,:,ii+3)=o0(:,:,:,ii+3)./o0(:,:,:,ii);
        o0(:,:,:,ii)=1;
    end
    indice=squeeze(sqrt(o0(:,:,:,4)));
    Ex=squeeze(e0(:,:,:,1)); Ey=squeeze(e0(:,:,:,2)); Ez=squeeze(e0(:,:,:,3));
    Hx=squeeze(e0(:,:,:,4)); Hy=squeeze(e0(:,:,:,5)); Hz=squeeze(e0(:,:,:,6));
    E = tif(in.pol==0, Ey, Ex);

    % Plot a cross section with trace_champ=1, x0=[], y0=0, z0=[]
    if isempty(x0)==1&&isempty(z0)==1
        filename = save_cross_section(in, xx, zz, indice, E, 'x', 'z', 'y', [-in.period_x/2,in.period_x/2], [0 inf], zou);
    elseif isempty(y0)==1&&isempty(z0)==1
        filename = save_cross_section(in, yy, zz, indice, E, 'y', 'z', 'x', [-in.period_y/2,in.period_y/2], [0 inf], zou);
    elseif isempty(x0)==1&&isempty(y0)==1
        filename = save_cross_section(in, yy, xx, indice, E, 'y', 'x', 'z', [-in.period_y/2,in.period_y/2], [-in.period_x/2,in.period_x/2], zou);
    end