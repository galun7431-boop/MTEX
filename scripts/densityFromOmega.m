function [theta,kappaX,kappaY] = densityFromOmega(w,g,valid,edgeX,edgeY,dx,dy)
% 原有李代数截断公式保持不变；两个差分阶段都遵守同一保护边。
    dX = maskedDerivative(w,g,valid,2,dx,edgeX,edgeY);
    dY = maskedDerivative(w,g,valid,1,dy,edgeX,edgeY);
    kappaX = dX+0.5*crossField(w,dX); kappaY = dY+0.5*crossField(w,dY);
    curl = maskedDerivative(kappaY,g,valid,2,dx,edgeX,edgeY) ...
        -maskedDerivative(kappaX,g,valid,1,dy,edgeX,edgeY);
    theta = curl+0.5*crossField(w,curl);
end
