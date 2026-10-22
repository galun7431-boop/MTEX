function validate_mtex_symmetry_api
% 需要已启动的 MTEX。用真实对象核实右乘约定与约化角，而非仅手工群。
    root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'scripts'));
    cs=crystalSymmetry('6/mmm');
    ori=orientation.byEuler(17*degree,23*degree,11*degree,cs);
    proper=cs.properGroup; ops=double(matrix(proper.rot));
    r=double(matrix(ori));
    for k=1:size(ops,3)
        equivalent=ori*proper.rot(k);
        assert(norm(double(matrix(equivalent))-r*ops(:,:,k),'fro')<1e-12);
        assert(angle(ori,equivalent)<1e-7);
    end
    % 用 MTEX 官方 angle 对照生产函数保存的约化角。
    q=orientation.byEuler(19*degree,23*degree,11*degree,cs);
    rr=cat(3,r,double(matrix(q))*ops(:,:,end),r,double(matrix(q)));
    o=struct('maxEdgeAngle',3*degree,'localOutlierAngle',2*degree, ...
        'minLocalNeighbors',2,'excludeLocalOutliers',false, ...
        'blockedX',false(2,1),'blockedY',false(1,2));
    f=prepareSymmetryContinuousField(rr,ones(4,1),(1:4)',[2,2],ops,true(4,1),zeros(4,1),o);
    assert(abs(f.reducedAngleY(1,1)-angle(ori,q))<1e-7);
    fprintf('MTEX 对象右乘约定及官方约化角对照：通过\n');
end
