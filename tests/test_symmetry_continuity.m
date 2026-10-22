function tests = test_symmetry_continuity
% 合成矩阵回归测试；不需要 EBSD 文件或 MTEX，但需要 MATLAB。
% 执行状态见 scripts/README.md；此文件需 MATLAB unittest。
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testCase.TestData.originalPath = path;
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'scripts'));
end

function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end

function testEquivalentRepresentativesBecomeContinuous(testCase)
% 同一恒定物理取向交替保存为不同 HCP 对称代表元，不能产生假梯度。
    [r,s,g,idx,good,score,o] = fixture;
    base = rx(17*pi/180)*rz(8*pi/180);
    for p = 1:9, r(:,:,p)=base*s(:,:,mod(p-1,12)+1); end
    f = prepareSymmetryContinuousField(r,g,idx,[3,3],s,good,score,o);
    for p = 1:9
        verifyLessThan(testCase,norm(f.Rcontinuous(:,:,p)-f.Rcontinuous(:,:,1),'fro'),1e-10);
        % 每次修正必须属于该测点原始取向的晶体对称轨道。
        k=f.symmetryIndex(p);
        verifyLessThan(testCase,norm(f.Rcontinuous(:,:,p)-f.Rraw(:,:,p)*s(:,:,k),'fro'),1e-10);
    end
    verifyTrue(testCase,all(f.edgeX(:)) && all(f.edgeY(:)));
    verifyEqual(testCase,max(f.componentGrid(:)),1);
end

function testPhysicalGradientIsNotSmoothed(testCase)
% 真实 0.5 度/像素梯度必须保留，不能被平均成常量。
    [r,s,g,idx,good,score,o] = fixture;
    for col = 1:3
        for row = 1:3
            p=sub2ind([3,3],row,col);
            r(:,:,p)=rz((col-1)*0.5*pi/180)*s(:,:,mod(p,12)+1);
        end
    end
    f=prepareSymmetryContinuousField(r,g,idx,[3,3],s,good,score,o);
    for row = 1:3
        p=sub2ind([3,3],row,1); q=sub2ind([3,3],row,2);
        delta=f.Rcontinuous(:,:,q)*f.Rcontinuous(:,:,p)';
        verifyLessThan(testCase,norm(delta-rz(0.5*pi/180),'fro'),1e-10);
    end
end

function testTrueAndTwinBoundariesRemainSeparate(testCase)
% 用约 86 度 a 轴旋转表示 Mg 孪晶量级的实际取向差；不进行孪晶分类。
    [r,s,g,idx,good,score,o] = fixture;
    g(4:9)=2;
    for p=4:9, r(:,:,p)=rx(86*pi/180); end
    f=prepareSymmetryContinuousField(r,g,idx,[3,3],s,good,score,o);
    verifyFalse(testCase,any(f.edgeX(:,1)));
    verifyEqual(testCase,max(f.componentGrid(:)),2);
    verifyEqual(testCase,f.grainGrid,reshape(g,3,3));
    a=f.Rcontinuous(:,:,1); b=f.Rcontinuous(:,:,4);
    % 独立换代表元不能把不等价的晶粒取向变成相同旋转。
    verifyGreaterThan(testCase,norm(a-b,'fro'),0.1);
end

function testLargePhysicalJumpOnlyFlagged(testCase)
    [r,s,g,idx,good,score,o] = fixture;
    for p=4:9, r(:,:,p)=rx(5*pi/180); end
    f=prepareSymmetryContinuousField(r,g,idx,[3,3],s,good,score,o);
    verifyTrue(testCase,all(f.edgeX(:)));
    verifyTrue(testCase,any(f.largePhysicalX(:)));
    verifyEqual(testCase,max(f.componentGrid(:)),1);
    verifyLessThan(testCase,f.maxInvarianceError,1e-7);
end

function testQualityHoleCannotBeUsedAsBridge(testCase)
    [r,s,g,idx,good,score,o] = fixture;
    good(4:6)=false;
    f=prepareSymmetryContinuousField(r,g,idx,[3,3],s,good,score,o);
    verifyFalse(testCase,any(f.validGrid(:,2)));
    verifyFalse(testCase,any(f.edgeX(:)));
    verifyEqual(testCase,max(f.componentGrid(:)),2);
    verifyTrue(testCase,all(isnan(reshape(f.Rcontinuous(:,:,4:6),[],1))));
end

function testExplicitBoundaryBlocksSmallJump(testCase)
% 保护边界不受取向差大小影响，即使差值很小也不得跨边选代表元。
    [r,s,g,idx,good,score,o] = fixture;
    o.blockedX(:,1)=true;
    f=prepareSymmetryContinuousField(r,g,idx,[3,3],s,good,score,o);
    verifyFalse(testCase,any(f.edgeX(:,1)));
    verifyEqual(testCase,max(f.componentGrid(:)),2);
end

function testInvalidMatricesAreExcludedRatherThanRepaired(testCase)
    [r,s,g,idx,good,score,o] = fixture;
    r(:,:,5)=diag([1,1,-1]); r(:,:,2)=NaN;
    f=prepareSymmetryContinuousField(r,g,idx,[3,3],s,good,score,o);
    verifyTrue(testCase,f.invalidRotation(5) && f.invalidRotation(2));
    verifyFalse(testCase,f.validPoint(5) || f.validPoint(2));
end

function testLocalSuspectIsFlaggedWithoutAutomaticSmoothing(testCase)
    [r,s,g,idx,good,score,o] = fixture;
    r(:,:,5)=rx(2.5*pi/180);
    f=prepareSymmetryContinuousField(r,g,idx,[3,3],s,good,score,o);
    verifyTrue(testCase,f.suspectLocal(5));
    verifyTrue(testCase,f.validPoint(5));
    % 代表元模块不允许根据局部角度剔除点。
    verifyLessThan(testCase,f.maxInvarianceError,1e-7);
end

function [r,s,g,idx,good,score,o] = fixture
% 六次 c 轴旋转及六次基面二次旋转组成 6/mmm 的正旋转子群 D6。
    s=zeros(3,3,12);
    for k=1:6
        s(:,:,k)=rz((k-1)*pi/3);
        s(:,:,k+6)=rz((k-1)*pi/3)*rx(pi);
    end
    r=repmat(eye(3),1,1,9); g=ones(9,1); idx=(1:9)';
    good=true(9,1); score=zeros(9,1);
    o.maxEdgeAngle=3*pi/180; o.localOutlierAngle=2*pi/180;
    o.minLocalNeighbors=3; o.excludeLocalOutliers=false;
    o.blockedX=false(3,2); o.blockedY=false(2,3);
end

function r = rz(a)
    r=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
end

function r = rx(a)
    r=[1,0,0;0,cos(a),-sin(a);0,sin(a),cos(a)];
end
