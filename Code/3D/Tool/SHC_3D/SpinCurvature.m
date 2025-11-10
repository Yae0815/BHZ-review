function [SpinCurvature,Ek] = SpinCurvature(kvec,ftn58sparse)
    % calculate the SpinCurvature of a ftn58sparse at a given k-point

    %% --- Initialization --- %%%
    BR   = ftn58sparse.BR;
    abc  = ftn58sparse.abc; 
    Norb = ftn58sparse.norb;
    ii   = ftn58sparse.ij(:,1);
    jj   = ftn58sparse.ij(:,2);
    tt   = ftn58sparse.tt;
    dd   = ftn58sparse.dd;
    gam  = 1e-20;                % broadening factor to avoid divergence at degeneracies
    sx   = (1/2)*[[0 1];[1 0]];
    sy   = (1/2)*[[0 -1i];[1i 0]];
    sz   = (1/2)*[[1 0];[0 -1]];

    %% transform to CARTESIAN coordinates !
    T    = [BR(1,:)*abc(1); BR(2,:)*abc(2); BR(3,:)*abc(3)]; 
    DD   = dd*T;
    DX   = DD(:,1);
    DY   = DD(:,2);
    DZ   = DD(:,3);
    
    %%% --- Actual Procedure --- %%%
    SpinCurvature = zeros(Norb,27);  
    H0            = full(sparse(ii,jj,exp(1i*dd*kvec').*tt,Norb,Norb));
    H0            = (H0 + H0')/2;
    [V,D]         = eig(H0);
    Ek            = diag(D);   
                
    %%% velocities (according to Hellmann-Feynman thm.) in CARTESIAN coordinates 
    vx  = 1i*full(sparse(ii,jj,DX.*exp(1i*dd*kvec').*tt,Norb,Norb)); % dH/dkx
    vy  = 1i*full(sparse(ii,jj,DY.*exp(1i*dd*kvec').*tt,Norb,Norb)); % dH/dky
    vz  = 1i*full(sparse(ii,jj,DZ.*exp(1i*dd*kvec').*tt,Norb,Norb)); % dH/dkz
    vx  = (vx + vx')/2;
    vy  = (vy + vy')/2;
    vz  = (vz + vz')/2;

    %% Matrix element for velocities %%%
    V_x = V'*vx*V; 
    V_x = (V_x + V_x')/2;  % make sure the velocity is Hermitian
    V_y = V'*vy*V; 
    V_y = (V_y + V_y')/2;  
    V_z = V'*vz*V; 
    V_z = (V_z + V_z')/2;  

    %%% spin currents %%%
    %% sx
    jxx = ( kron(sx,eye(Norb/2))*vx + vx*kron(sx,eye(Norb/2)) )/2;
    jyx = ( kron(sx,eye(Norb/2))*vy + vy*kron(sx,eye(Norb/2)) )/2;
    jzx = ( kron(sx,eye(Norb/2))*vz + vz*kron(sx,eye(Norb/2)) )/2;
    %% sy
    jxy = ( kron(sy,eye(Norb/2))*vx + vx*kron(sy,eye(Norb/2)) )/2;
    jyy = ( kron(sy,eye(Norb/2))*vy + vy*kron(sy,eye(Norb/2)) )/2;
    jzy = ( kron(sy,eye(Norb/2))*vz + vz*kron(sy,eye(Norb/2)) )/2;
    %% sz
    jxz = ( kron(sz,eye(Norb/2))*vx + vx*kron(sz,eye(Norb/2)) )/2;
    jyz = ( kron(sz,eye(Norb/2))*vy + vy*kron(sz,eye(Norb/2)) )/2;
    jzz = ( kron(sz,eye(Norb/2))*vz + vz*kron(sz,eye(Norb/2)) )/2;
    
    %%% Matrix element for spin currents %%%
    %% sx
    J_xx = V'*jxx*V;
    J_xx = (J_xx + J_xx')/2;
    J_yx = V'*jyx*V;
    J_yx = (J_yx + J_yx')/2;
    J_zx = V'*jzx*V;
    J_zx = (J_zx + J_zx')/2;
    %% sy
    J_xy = V'*jxy*V;
    J_xy = (J_xy + J_xy')/2;
    J_yy = V'*jyy*V;
    J_yy = (J_yy + J_yy')/2;
    J_zy = V'*jzy*V;
    J_zy = (J_zy + J_zy')/2;
    %% sz
    J_xz = V'*jxz*V;
    J_xz = (J_xz + J_xz')/2;
    J_yz = V'*jyz*V;
    J_yz = (J_yz + J_yz')/2;
    J_zz = V'*jzz*V;
    J_zz = (J_zz + J_zz')/2;
    
    curvature_xxx = zeros(Norb);
    curvature_xxy = zeros(Norb); 
    curvature_xxz = zeros(Norb); 
    curvature_yxx = zeros(Norb);
    curvature_yxy = zeros(Norb); 
    curvature_yxz = zeros(Norb); 
    curvature_zxx = zeros(Norb);
    curvature_zxy = zeros(Norb); 
    curvature_zxz = zeros(Norb); 
    curvature_xyx = zeros(Norb);
    curvature_xyy = zeros(Norb); 
    curvature_xyz = zeros(Norb); 
    curvature_yyx = zeros(Norb);
    curvature_yyy = zeros(Norb); 
    curvature_yyz = zeros(Norb); 
    curvature_zyx = zeros(Norb);
    curvature_zyy = zeros(Norb); 
    curvature_zyz = zeros(Norb);
    curvature_xzx = zeros(Norb);
    curvature_xzy = zeros(Norb); 
    curvature_xzz = zeros(Norb); 
    curvature_yzx = zeros(Norb);
    curvature_yzy = zeros(Norb); 
    curvature_yzz = zeros(Norb); 
    curvature_zzx = zeros(Norb);
    curvature_zzy = zeros(Norb); 
    curvature_zzz = zeros(Norb);
    %%  Calulate the Berry Curvature, # Grosso, VIII (62) %%%   
    for nb1 = 1:Norb                         % no contribution for nb2=nb1 
        for nb2 = nb1+1:Norb                 % only upper-triangular part of the matrix
            diff_energy = Ek(nb1) - Ek(nb2); 
    
            %% avoid degeneracies, where Abelian curvature should be vanished instead of being diverged
            %% # B. Andrei Bernevig et al., TI & TSC, ch.2.3
            denominator_p   = 1/(diff_energy^2 + 1i*gam);
            denominator_n   = 1/(diff_energy^2 - 1i*gam);
          
            %% every component of the Im[ <nb1|v|nb2> x <nb2|v|nb1> ]
            curvature_xxx(nb1,nb2) = imag( J_xx(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_xxy(nb1,nb2) = imag( J_xx(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_xxz(nb1,nb2) = imag( J_xx(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_yxx(nb1,nb2) = imag( J_yx(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_yxy(nb1,nb2) = imag( J_yx(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_yxz(nb1,nb2) = imag( J_yx(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_zxx(nb1,nb2) = imag( J_zx(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_zxy(nb1,nb2) = imag( J_zx(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_zxz(nb1,nb2) = imag( J_zx(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );
            
            curvature_xyx(nb1,nb2) = imag( J_xy(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_xyy(nb1,nb2) = imag( J_xy(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_xyz(nb1,nb2) = imag( J_xy(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_yyx(nb1,nb2) = imag( J_yy(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_yyy(nb1,nb2) = imag( J_yy(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_yyz(nb1,nb2) = imag( J_yy(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_zyx(nb1,nb2) = imag( J_zy(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_zyy(nb1,nb2) = imag( J_zy(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_zyz(nb1,nb2) = imag( J_zy(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );

            curvature_xzx(nb1,nb2) = imag( J_xz(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_xzy(nb1,nb2) = imag( J_xz(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_xzz(nb1,nb2) = imag( J_xz(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_yzx(nb1,nb2) = imag( J_yz(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_yzy(nb1,nb2) = imag( J_yz(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_yzz(nb1,nb2) = imag( J_yz(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_zzx(nb1,nb2) = imag( J_zz(nb1,nb2)*V_x(nb2,nb1)*( denominator_p + denominator_n ) );  
            curvature_zzy(nb1,nb2) = imag( J_zz(nb1,nb2)*V_y(nb2,nb1)*( denominator_p + denominator_n ) );
            curvature_zzz(nb1,nb2) = imag( J_zz(nb1,nb2)*V_z(nb2,nb1)*( denominator_p + denominator_n ) );
        end  
    end
    
    %% compensate the lower-triangular part ( "-" because 1i' == -1i ) 
    curvature_xxx = curvature_xxx - curvature_xxx.'; % .' => transport (no conjugate)
    curvature_xxy = curvature_xxy - curvature_xxy.';
    curvature_xxz = curvature_xxz - curvature_xxz.';
    curvature_yxx = curvature_yxx - curvature_yxx.';
    curvature_yxy = curvature_yxy - curvature_yxy.';
    curvature_yxz = curvature_yxz - curvature_yxz.';
    curvature_zxx = curvature_zxx - curvature_zxx.';
    curvature_zxy = curvature_zxy - curvature_zxy.';
    curvature_zxz = curvature_zxz - curvature_zxz.';

    curvature_xyx = curvature_xyx - curvature_xyx.';
    curvature_xyy = curvature_xyy - curvature_xyy.';
    curvature_xyz = curvature_xyz - curvature_xyz.';
    curvature_yyx = curvature_yyx - curvature_yyx.';
    curvature_yyy = curvature_yyy - curvature_yyy.';
    curvature_yyz = curvature_yyz - curvature_yyz.';
    curvature_zyx = curvature_zyx - curvature_zyx.';
    curvature_zyy = curvature_zyy - curvature_zyy.';
    curvature_zyz = curvature_zyz - curvature_zyz.';

    curvature_xzx = curvature_xzx - curvature_xzx.';
    curvature_xzy = curvature_xzy - curvature_xzy.';
    curvature_xzz = curvature_xzz - curvature_xzz.';
    curvature_yzx = curvature_yzx - curvature_yzx.';
    curvature_yzy = curvature_yzy - curvature_yzy.';
    curvature_yzz = curvature_yzz - curvature_yzz.';
    curvature_zzx = curvature_zzx - curvature_zzx.';
    curvature_zzy = curvature_zzy - curvature_zzy.';
    curvature_zzz = curvature_zzz - curvature_zzz.';
    
    %% sum over the orbital 2 to get the complete Berry Curvature (xxx means: cross(Jxx,Vx) )
    SpinCurvature(:,1)  = sum(curvature_xxx,2);
    SpinCurvature(:,2)  = sum(curvature_xxy,2);
    SpinCurvature(:,3)  = sum(curvature_xxz,2);
    SpinCurvature(:,4)  = sum(curvature_yxx,2);
    SpinCurvature(:,5)  = sum(curvature_yxy,2);
    SpinCurvature(:,6)  = sum(curvature_yxz,2);
    SpinCurvature(:,7)  = sum(curvature_zxx,2);
    SpinCurvature(:,8)  = sum(curvature_zxy,2);
    SpinCurvature(:,9)  = sum(curvature_zxz,2);

    SpinCurvature(:,10) = sum(curvature_xyx,2);
    SpinCurvature(:,11) = sum(curvature_xyy,2);
    SpinCurvature(:,12) = sum(curvature_xyz,2);
    SpinCurvature(:,13) = sum(curvature_yyx,2);
    SpinCurvature(:,14) = sum(curvature_yyy,2);
    SpinCurvature(:,15) = sum(curvature_yyz,2);
    SpinCurvature(:,16) = sum(curvature_zyx,2);
    SpinCurvature(:,17) = sum(curvature_zyy,2);
    SpinCurvature(:,18) = sum(curvature_zyz,2);

    SpinCurvature(:,19) = sum(curvature_xzx,2);
    SpinCurvature(:,20) = sum(curvature_xzy,2);
    SpinCurvature(:,21) = sum(curvature_xzz,2);
    SpinCurvature(:,22) = sum(curvature_yzx,2);
    SpinCurvature(:,23) = sum(curvature_yzy,2);
    SpinCurvature(:,24) = sum(curvature_yzz,2);
    SpinCurvature(:,25) = sum(curvature_zzx,2);
    SpinCurvature(:,26) = sum(curvature_zzy,2);
    SpinCurvature(:,27) = sum(curvature_zzz,2);
end
