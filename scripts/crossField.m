function C = crossField(A,B)
% 对两个 ny×nx×3 向量场在每个像素处逐点计算三维叉乘。
    % 预分配结果数组 C，大小与 A 相同。
    C = nan(size(A));
    % 拆出 A 在三个轴向的二维分量。
    ax=A(:,:,1); ay=A(:,:,2); az=A(:,:,3);
    % 拆出 B 在三个轴向的二维分量。
    bx=B(:,:,1); by=B(:,:,2); bz=B(:,:,3);
    % 按向量叉乘公式计算第一分量 ay*bz-az*by。
    C(:,:,1) = ay.*bz-az.*by;
    % 按向量叉乘公式计算第二分量 az*bx-ax*bz。
    C(:,:,2) = az.*bx-ax.*bz;
    % 按向量叉乘公式计算第三分量 ax*by-ay*bx。
    C(:,:,3) = ax.*by-ay.*bx;
% 结束逐像素叉乘函数。
end
