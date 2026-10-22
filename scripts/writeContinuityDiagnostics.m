function writeContinuityDiagnostics(f,validResult,g,x,y,outDir,baseName,degree)
% 诊断不是正确性证明：约化角应保持不变，而不是越小越好。
    fig=figure('Color','w','Name','代表元选择及计算有效性');
    tiledlayout(2,3);
    before=[f.reducedAngleX(:);f.reducedAngleY(:)]/degree;
    after=[f.reducedAfterX(:);f.reducedAfterY(:)]/degree;
    nexttile; histogram(before(isfinite(before))); hold on;
    histogram(after(isfinite(after))); hold off;
    legend('原始约化错配','处理后约化错配'); xlabel('度');
    title(sprintf('约化角不变性误差 %.3g rad',f.maxInvarianceError));
    raw=[f.rawIncrementX(:);f.rawIncrementY(:)]/degree;
    selected=[f.residualX(:);f.residualY(:)]/degree;
    nexttile; histogram(raw(isfinite(raw))); hold on;
    histogram(selected(isfinite(selected))); hold off;
    legend('原始旋转增量','选定代表元增量'); xlabel('度');
    nexttile; imagesc(x,y,double(validResult)); axis image; set(gca,'YDir','normal');
    title('最终向错密度有效像素'); colorbar;
    % 未选相/未索引/缺失测量以 g=0 表示；具体原因保留输入相报告与质量报告。
    mask=zeros(size(g)); mask(g==0)=1; mask(g>0 & ~f.validGrid)=2;
    mask(f.validGrid & ~validResult)=3;
    nexttile; imagesc(x,y,mask); axis image; set(gca,'YDir','normal'); colorbar;
    title('0可用 1缺失/未选相 2质量/代表元 3模板不足');
    nexttile; imagesc(x,y,f.phaseGrid); axis image; set(gca,'YDir','normal'); colorbar;
    title('原始相 ID（含未处理相；0为无测点）');
    % 显式显示晶界/相界边的位置；两端都不跨边求导。
    boundary=zeros(size(g));
    boundary(:,1:end-1)=boundary(:,1:end-1)+double(f.grainBoundaryX);
    boundary(1:end-1,:)=boundary(1:end-1,:)+double(f.grainBoundaryY);
    boundary(:,1:end-1)=boundary(:,1:end-1)+2*double(f.phaseBoundaryX);
    boundary(1:end-1,:)=boundary(1:end-1,:)+2*double(f.phaseBoundaryY);
    nexttile; imagesc(x,y,boundary); axis image; set(gca,'YDir','normal'); colorbar;
    title('边界位置：晶界权重1，相界权重2（可叠加）');
    saveas(fig,fullfile(outDir,[baseName '_Symmetry_diagnostics.png']));
    ids=unique(g(g>0)); usable=zeros(size(ids)); masked=usable; anomalousEdges=usable;
    for k=1:numel(ids)
        inside=g==ids(k); usable(k)=nnz(inside & validResult);
        masked(k)=nnz(inside & ~validResult);
        anomalousEdges(k)=nnz(f.loopConflictX & inside(:,1:end-1)) ...
            +nnz(f.loopConflictY & inside(1:end-1,:));
    end
    stats=table(ids,usable,masked,anomalousEdges, ...
        'VariableNames',{'grainId','validPixels','maskedPixels','conflictEdges'});
    writetable(stats,fullfile(outDir,[baseName '_Grain_diagnostics.csv']));
end
