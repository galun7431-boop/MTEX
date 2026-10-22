function run_symmetry_validation
% 不依赖 MTEX 的生产矩阵算法测试 A-D 及真实 5 度回归。
% MATLAB: addpath('tests'); run_symmetry_validation
    root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'scripts'));
    s=zeros(3,3,12);
    for k=1:6
        s(:,:,k)=rz((k-1)*pi/3); s(:,:,k+6)=s(:,:,k)*rx(pi);
    end
    ny=7; nx=11; n=ny*nx; idx=(1:n)'; g=ones(n,1);
    o=struct('maxEdgeAngle',3*pi/180,'localOutlierAngle',2*pi/180, ...
        'minLocalNeighbors',3,'excludeLocalOutliers',false, ...
        'blockedX',false(ny,nx-1),'blockedY',false(ny-1,nx));
    for total=[0,5]*pi/180
        r=zeros(3,3,n); truth=r;
        for col=1:nx
            for row=1:ny
                p=sub2ind([ny,nx],row,col);
                truth(:,:,p)=rz(total*(col-1)/(nx-1));
                r(:,:,p)=truth(:,:,p)*s(:,:,mod(p-1,12)+1);
            end
        end
        f=prepareSymmetryContinuousField(r,g,idx,[ny,nx],s,true(n,1),zeros(n,1),o);
        assert(f.maxInvarianceError<1e-7);
        assert(all(f.edgeX(:)) && all(f.edgeY(:)));
        assert(all(diff(f.energyHistory)<=1e-10));
        for p=1:n
            assert(norm(f.Rcontinuous(:,:,p)-f.Rraw(:,:,p)*s(:,:,f.symmetryIndex(p)),'fro')<1e-12);
        end
        delta=f.Rcontinuous(:,:,end)*f.Rcontinuous(:,:,1)';
        assert(norm(delta-rz(total),'fro')<1e-10,'真实总取向变化被改变');
        w=nan(ny,nx,3);
        for p=1:n
            [i,j]=ind2sub([ny,nx],p); w(i,j,:)=reshape(localRotmLog(f.Rcontinuous(:,:,p)),1,1,3);
        end
        [theta,kx,~]=densityFromOmega(w,f.grainGrid,f.validGrid,f.edgeX,f.edgeY,1,1);
        expected=total/(nx-1); z=kx(:,:,3);
        assert(max(abs(z(:)-expected))<1e-10,'梯度不符');
        assert(all(isfinite(theta(:))) && max(abs(theta(:)))<1e-10,'产生虚假向错密度');
        fprintf('A/B/5度场 total=%.1f deg：通过，theta最大值 %.3g\n',total*180/pi,max(abs(theta(:))));
    end
    % 局部真实 5 度跃迁：角度本身不构成删边理由。
    for p=1:n
        [~,j]=ind2sub([ny,nx],p); r(:,:,p)=rz((j>5)*5*pi/180)*s(:,:,mod(p-1,12)+1);
    end
    f=prepareSymmetryContinuousField(r,g,idx,[ny,nx],s,true(n,1),zeros(n,1),o);
    assert(all(f.edgeX(:)) && all(f.validPoint));
    assert(any(f.largePhysicalX(:)));
    % C：两晶粒独立，真实 20 度约化错配不消失。
    gg=reshape(g,ny,nx); gg(:,6:end)=2;
    for p=1:n
        r(:,:,p)=rz((gg(p)==2)*20*pi/180)*s(:,:,mod(p-1,12)+1);
    end
    f=prepareSymmetryContinuousField(r,gg(:),idx,[ny,nx],s,true(n,1),zeros(n,1),o);
    assert(~any(f.edgeX(:,5)) && max(f.componentGrid(:))==2);
    a=f.Rcontinuous(:,:,1); b=f.Rcontinuous(:,:,end);
    angles=zeros(12,1);
    for k=1:12, angles(k)=acos(max(-1,min(1,(trace(a'*(b*s(:,:,k)))-1)/2))); end
    assert(abs(min(angles)-20*pi/180)<1e-10);
    w=nan(ny,nx,3);
    for p=1:n
        [i,j]=ind2sub([ny,nx],p); w(i,j,:)=reshape(localRotmLog(f.Rcontinuous(:,:,p)),1,1,3);
    end
    theta=densityFromOmega(w,f.grainGrid,f.validGrid,f.edgeX,f.edgeY,1,1);
    assert(all(isfinite(theta(:))) && max(abs(theta(:)))<1e-10);
    fprintf('C 两晶粒独立且保留晶界错配、晶界两侧同侧差分：通过\n');
    % D：第二相只使用其二次旋转群，空洞/相界无导数边。
    phase=ones(ny,nx); phase(:,6:end)=2;
    groups={s,cat(3,eye(3),rx(pi))};
    for ph=1:2
        selected=find(phase(:)==ph); ss=groups{ph}; rr=zeros(3,3,numel(selected));
        for p=1:numel(selected), rr(:,:,p)=rz(0.12)*ss(:,:,mod(p-1,size(ss,3))+1); end
        op=o; op.phaseId=ph*ones(numel(selected),1);
        f=prepareSymmetryContinuousField(rr,ones(numel(selected),1),selected,[ny,nx], ...
            ss,true(numel(selected),1),zeros(numel(selected),1),op);
        assert(~any(f.edgeX(:,5)) && max(f.symmetryIndex)<=size(ss,3));
        assert(f.maxInvarianceError<1e-7);
    end
    fprintf('D 不同晶相分别使用自身对称群：通过\n');
end
function r=rz(a)
    r=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
end
function r=rx(a)
    r=[1,0,0;0,cos(a),-sin(a);0,sin(a),cos(a)];
end
