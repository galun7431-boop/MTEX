function field = prepareSymmetryContinuousField(R, grainId, linIdx, gridSize, ...
        symOps, qualityGood, seedScore, options)
% PREPARESYMMETRYCONTINUOUSFIELD 在矩阵对数前选择晶粒内连续的等价取向。
% R: 3x3xN，晶体到试样的主动旋转；候选严格为 R(:,:,p)*symOps(:,:,k)。
% 仅使用晶体正旋转群，不施加试样对称性、不平均、不平滑原始取向。
% grainId/linIdx/qualityGood/seedScore: N 个测点；gridSize=[ny nx]。
% linIdx 必须由调用方验证的轴对齐矩形扫描坐标构造，不能按导入顺序生成。
% options: maxEdgeAngle（弧度，仅诊断），localOutlierAngle（弧度），
% minLocalNeighbors，excludeLocalOutliers，blockedX(ny,nx-1)、blockedY(ny-1,nx)。
% 输出边掩膜同时供两阶差分使用；断开的连通域独立选种子。
% 确定性图传播及离散坐标下降不保证全局最优；回路不一致的边
% 会被标记并禁用，不能为得到光滑图像而跨真实晶界/孪晶界连接。
    % 本接口处理一个已明确选择的相；禁止混相后共用 symOps。
    defaults = struct('maxIterations',50,'energyTolerance',1e-12, ...
        'ambiguityTolerance',1e-10,'invarianceTolerance',1e-7);
    names=fieldnames(defaults);
    for k=1:numel(names)
        if ~isfield(options,names{k}), options.(names{k})=defaults.(names{k}); end
    end
    assert(all(isfinite([options.energyTolerance,options.ambiguityTolerance,options.invarianceTolerance])) ...
        && options.energyTolerance>=0 && options.ambiguityTolerance>=0 && options.invarianceTolerance>0);
    assert(options.maxIterations>=1 && options.maxIterations==floor(options.maxIterations));
    if isfield(options,'phaseId')
        assert(numel(options.phaseId)==size(R,3) && numel(unique(options.phaseId))==1, ...
            '本函数只接受一个晶相；请按相分别调用并传入该相的晶体对称群。');
    end
    assert(~options.excludeLocalOutliers, ...
        '代表元选择不执行异常点剔除；请在独立质量模块中验证并处理。');
    n = size(R,3); ny = gridSize(1); nx = gridSize(2);
    grainId = grainId(:); linIdx = linIdx(:);
    qualityGood = logical(qualityGood(:)); seedScore = seedScore(:);
    assert(size(R,1)==3 && size(R,2)==3 && numel(grainId)==n ...
        && numel(linIdx)==n && numel(qualityGood)==n && numel(seedScore)==n, ...
        '取向、晶粒编号、坐标与质量数组的长度必须一致。');
    assert(numel(unique(linIdx))==n && all(linIdx>=1 & linIdx<=ny*nx & linIdx==floor(linIdx)), ...
        '网格索引必须唯一且不越界。');
    assert(isscalar(options.maxEdgeAngle) && isfinite(options.maxEdgeAngle) ...
        && options.maxEdgeAngle>0 && options.maxEdgeAngle<pi, '增量诊断阈值须在 (0,pi) 内。');
    assert(isscalar(options.localOutlierAngle) && isfinite(options.localOutlierAngle) ...
        && options.localOutlierAngle>0 && isscalar(options.minLocalNeighbors) ...
        && isfinite(options.minLocalNeighbors) && options.minLocalNeighbors>=2 ...
        && options.minLocalNeighbors<=4 && options.minLocalNeighbors==floor(options.minLocalNeighbors), '局部异常点参数无效。');
    assert(isequal(size(options.blockedX),[ny,nx-1]) ...
        && isequal(size(options.blockedY),[ny-1,nx]), '保护边界掩膜尺寸不匹配。');
    assert(all(ismember(options.blockedX(:),[0,1])) ...
        && all(ismember(options.blockedY(:),[0,1])), '保护边界只允许逻辑值。');
    assert(all(isfinite(symOps(:))) && size(symOps,1)==3 && size(symOps,2)==3 ...
        && size(symOps,3)>=1, '晶体正旋转群为空或维度错误。');
    for k = 1:size(symOps,3)
        s = symOps(:,:,k);
        assert(norm(s'*s-eye(3),'fro')<1e-8 && abs(det(s)-1)<1e-8, ...
            '对称操作必须为 SO(3) 正旋转，不能包含反射。');
    end
    field.Rraw = nan(size(R)); field.Rcontinuous = nan(size(R));
    field.invalidRotation = false(n,1);
    for p = 1:n
        rp = R(:,:,p);
        % 只修正舍入误差；不把严重异常矩阵投影成看似有效的取向。
        if any(~isfinite(rp(:))) || norm(rp'*rp-eye(3),'fro')>1e-5 ...
                || abs(det(rp)-1)>1e-5
            field.invalidRotation(p) = true;
            continue;
        end
        % 原样保留输入，不通过 SO(3) 投影改变测量取向。
        field.Rraw(:,:,p) = rp;
    end
    initialGood = qualityGood & ~field.invalidRotation & isfinite(grainId) & grainId>0;
    pointMap = zeros(ny,nx); pointMap(linIdx) = (1:n)';
    field.grainGrid = zeros(ny,nx);
    validGrain = isfinite(grainId) & grainId>0;
    field.grainGrid(linIdx(validGrain)) = grainId(validGrain);
    % 先在原始取向上求晶体对称约化后的邻接角，不依赖任何对数。
    angleX = nan(ny,nx-1); angleY = nan(ny-1,nx);
    for row = 1:ny
        for col = 1:nx
            p = pointMap(row,col);
            if p==0 || ~initialGood(p), continue; end
            if col<nx && ~options.blockedX(row,col)
                q = pointMap(row,col+1);
                if q>0 && initialGood(q) && grainId(p)==grainId(q)
                    [~,~,angleX(row,col)] = closestEquivalent( ...
                        field.Rraw(:,:,q),field.Rraw(:,:,p),symOps);
                end
            end
            if row<ny && ~options.blockedY(row,col)
                q = pointMap(row+1,col);
                if q>0 && initialGood(q) && grainId(p)==grainId(q)
                    [~,~,angleY(row,col)] = closestEquivalent( ...
                        field.Rraw(:,:,q),field.Rraw(:,:,p),symOps);
                end
            end
        end
    end
    % 只比较同晶粒且非保护边界的邻点。异常标记是待检查线索，不是质量真值。
    field.suspectLocal = false(n,1);
    for p = find(initialGood)'
        [row,col] = ind2sub([ny,nx],linIdx(p));
        a = [];
        if col>1, a(end+1) = angleX(row,col-1); end %#ok<AGROW>
        if col<nx, a(end+1) = angleX(row,col); end %#ok<AGROW>
        if row>1, a(end+1) = angleY(row-1,col); end %#ok<AGROW>
        if row<ny, a(end+1) = angleY(row,col); end %#ok<AGROW>
        a = a(isfinite(a));
        field.suspectLocal(p) = numel(a)>=options.minLocalNeighbors ...
            && median(a)>options.localOutlierAngle;
    end
    field.validPoint = initialGood;
    field.validGrid = false(ny,nx); field.validGrid(linIdx) = field.validPoint;
    field.edgeX = isfinite(angleX) ...
        & field.validGrid(:,1:end-1) & field.validGrid(:,2:end);
    field.edgeY = isfinite(angleY) ...
        & field.validGrid(1:end-1,:) & field.validGrid(2:end,:);
    field.symmetryIndex = zeros(n,1); field.componentGrid = zeros(ny,nx);
    % 高质量点优先作种子；相同质量按测点顺序，便于复现。
    seedScore(~isfinite(seedScore)) = -Inf;
    [~,seedOrder] = sort(seedScore,'descend');
    visited = false(n,1); queue = zeros(n,1); component = 0;
    for seed = seedOrder(:)'
        if visited(seed) || ~field.validPoint(seed), continue; end
        component = component+1;
        % 种子取离单位旋转最近的等价元，尽量远离主值对数的 pi 分支。
        [field.Rcontinuous(:,:,seed),field.symmetryIndex(seed)] = ...
            closestEquivalent(field.Rraw(:,:,seed),eye(3),symOps);
        head = 1; tail = 1; queue(1) = seed; visited(seed) = true;
        field.componentGrid(linIdx(seed)) = component;
        while head<=tail
            p = queue(head); head = head+1;
            [row,col] = ind2sub([ny,nx],linIdx(p));
            neighbors = zeros(1,4); count = 0;
            if col>1 && field.edgeX(row,col-1)
                count=count+1; neighbors(count)=pointMap(row,col-1);
            end
            if col<nx && field.edgeX(row,col)
                count=count+1; neighbors(count)=pointMap(row,col+1);
            end
            if row>1 && field.edgeY(row-1,col)
                count=count+1; neighbors(count)=pointMap(row-1,col);
            end
            if row<ny && field.edgeY(row,col)
                count=count+1; neighbors(count)=pointMap(row+1,col);
            end
            for q = neighbors(1:count)
                if visited(q), continue; end
                [field.Rcontinuous(:,:,q),field.symmetryIndex(q)] = ...
                    closestEquivalent(field.Rraw(:,:,q),field.Rcontinuous(:,:,p),symOps);
                visited(q)=true; tail=tail+1; queue(tail)=q;
                field.componentGrid(linIdx(q))=component;
            end
        end
    end
    % 建立真实矩形网格四邻接图；边已经过坐标、同晶粒与质量检查。
    neighbors=cell(n,1);
    [rr,cc]=find(field.edgeX);
    for e=1:numel(rr)
        p=pointMap(rr(e),cc(e)); q=pointMap(rr(e),cc(e)+1);
        neighbors{p}(end+1)=q; neighbors{q}(end+1)=p;
    end
    [rr,cc]=find(field.edgeY);
    for e=1:numel(rr)
        p=pointMap(rr(e),cc(e)); q=pointMap(rr(e)+1,cc(e));
        neighbors{p}(end+1)=q; neighbors{q}(end+1)=p;
    end
    % 离散坐标下降：只枚举 Rraw*S，最小化整个邻接图的角平方和。
    % 固定每个连通域的种子以移除整体对称规范自由度；不平均取向。
    anchors=false(n,1);
    for c=1:component
        members=seedOrder(field.componentGrid(linIdx(seedOrder))==c);
        anchors(members(1))=true;
    end
    field.energyHistory=zeros(options.maxIterations+1,1);
    field.energyHistory(1)=graphEnergy(field.Rcontinuous,neighbors);
    field.nonconverged=false(n,1); changed=false(n,1);
    for iteration=1:options.maxIterations
        changed(:)=false;
        for p=find(field.validPoint & ~anchors)'
            nb=neighbors{p}; if isempty(nb), continue; end
            costs=candidateCosts(field.Rraw(:,:,p),symOps,field.Rcontinuous,nb);
            [best,k]=min(costs); current=costs(field.symmetryIndex(p));
            if current-best>options.energyTolerance
                field.symmetryIndex(p)=k;
                field.Rcontinuous(:,:,p)=field.Rraw(:,:,p)*symOps(:,:,k);
                changed(p)=true;
            end
        end
        field.energyHistory(iteration+1)=graphEnergy(field.Rcontinuous,neighbors);
        if ~any(changed), break; end
    end
    field.iterations=iteration;
    field.energyHistory=field.energyHistory(1:iteration+1);
    if any(changed)
        % 达迭代上限时屏蔽受影响连通域，而非伪称收敛。
        badComponents=unique(field.componentGrid(linIdx(changed)));
        field.nonconverged=ismember(field.componentGrid(linIdx),badComponents);
    end
    field.ambiguous=false(n,1);
    for p=find(field.validPoint & ~anchors)'
        if isempty(neighbors{p}), continue; end
        costs=sort(candidateCosts(field.Rraw(:,:,p),symOps,field.Rcontinuous,neighbors{p}));
        if numel(costs)>1
            field.ambiguous(p)=costs(2)-costs(1)<=options.ambiguityTolerance;
        end
    end
    % 核查所有边，包括传播树未经过的回路边；残余跳变保留并禁用差分。
    field.residualX = nan(ny,nx-1); field.residualY = nan(ny-1,nx);
    [rows,cols] = find(field.edgeX);
    for e = 1:numel(rows)
        i=rows(e); j=cols(e); p=pointMap(i,j); q=pointMap(i,j+1);
        field.residualX(i,j)=rotationDistance(field.Rcontinuous(:,:,p),field.Rcontinuous(:,:,q));
    end
    [rows,cols] = find(field.edgeY);
    for e = 1:numel(rows)
        i=rows(e); j=cols(e); p=pointMap(i,j); q=pointMap(i+1,j);
        field.residualY(i,j)=rotationDistance(field.Rcontinuous(:,:,p),field.Rcontinuous(:,:,q));
    end
    field.loopConflictX = field.edgeX & field.residualX>angleX+options.invarianceTolerance;
    field.loopConflictY = field.edgeY & field.residualY>angleY+options.invarianceTolerance;
    field.edgeX(field.loopConflictX)=false; field.edgeY(field.loopConflictY)=false;
    field.reducedAngleX=angleX; field.reducedAngleY=angleY;
    field.reducedAfterX=nan(size(angleX)); field.reducedAfterY=nan(size(angleY));
    field.rawIncrementX=nan(size(angleX)); field.rawIncrementY=nan(size(angleY));
    for row=1:ny
        for col=1:nx
            p=pointMap(row,col); if p==0 || ~field.validPoint(p), continue; end
            if col<nx && isfinite(angleX(row,col))
                q=pointMap(row,col+1);
                if field.validPoint(q)
                    [~,~,field.reducedAfterX(row,col)]=closestEquivalent( ...
                        field.Rcontinuous(:,:,q),field.Rcontinuous(:,:,p),symOps);
                    field.rawIncrementX(row,col)=rotationDistance(field.Rraw(:,:,p),field.Rraw(:,:,q));
                end
            end
            if row<ny && isfinite(angleY(row,col))
                q=pointMap(row+1,col);
                if field.validPoint(q)
                    [~,~,field.reducedAfterY(row,col)]=closestEquivalent( ...
                        field.Rcontinuous(:,:,q),field.Rcontinuous(:,:,p),symOps);
                    field.rawIncrementY(row,col)=rotationDistance(field.Rraw(:,:,p),field.Rraw(:,:,q));
                end
            end
        end
    end
    differences=[field.reducedAfterX(:)-angleX(:);field.reducedAfterY(:)-angleY(:)];
    differences=differences(isfinite(differences));
    field.maxInvarianceError=max([0;abs(differences)]);
    assert(field.maxInvarianceError<=options.invarianceTolerance, ...
        '对称约化错配角发生变化，停止计算，请检查对称群与复合约定。');
    % 大物理错配只诊断，不据此删点或删边。
    field.largePhysicalX=angleX>options.maxEdgeAngle;
    field.largePhysicalY=angleY>options.maxEdgeAngle;
    unreliable=false(ny,nx);
    unreliable(linIdx)=field.ambiguous | field.nonconverged;
    field.edgeX=field.edgeX & ~unreliable(:,1:end-1) & ~unreliable(:,2:end);
    field.edgeY=field.edgeY & ~unreliable(1:end-1,:) & ~unreliable(2:end,:);
    field.validGrid=field.validGrid & ~unreliable;
    field.validPoint=field.validPoint & ~unreliable(linIdx);
end

function [candidate,index,a] = closestEquivalent(r,reference,symOps)
% MTEX 取向将晶体坐标映射到试样坐标，故晶体对称操作右乘。
    candidates = reshape(r*reshape(symOps,3,[]),9,[]);
    scores = sum(candidates.*reshape(reference,9,1),1);
    [best,index] = max(scores);
    candidate = reshape(candidates(:,index),3,3);
    a = acos(max(-1,min(1,(best-1)/2)));
end

function a = rotationDistance(aMatrix,bMatrix)
% 未作对称约化的 SO(3) 距离，用于检测已选代表元的残余不连续。
    a = acos(max(-1,min(1,(sum(sum(aMatrix.*bMatrix))-1)/2)));
end

function costs = candidateCosts(r,symOps,selected,neighbors)
% 矩阵距离天然不受四元数 q/-q 的符号影响。
    costs=zeros(1,size(symOps,3));
    for k=1:numel(costs)
        candidate=r*symOps(:,:,k);
        for q=neighbors
            costs(k)=costs(k)+rotationDistance(candidate,selected(:,:,q))^2;
        end
    end
end

function energy = graphEnergy(selected,neighbors)
    energy=0;
    for p=1:numel(neighbors)
        for q=neighbors{p}
            if q>p, energy=energy+rotationDistance(selected(:,:,p),selected(:,:,q))^2; end
        end
    end
end
