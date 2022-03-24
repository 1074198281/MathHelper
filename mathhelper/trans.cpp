    //根据三个角度求旋转变换矩阵
    //绕pmatrix的三个轴
    Matx33d AngleTransMatrixp(double dAngle[3], Matx33d& pmatrix)
    {
        Matx13d vec0(pmatrix(0, 0), pmatrix(0, 1), pmatrix(0, 2));
        Matx13d vec1(pmatrix(1, 0), pmatrix(1, 1), pmatrix(1, 2));
        Matx13d vec2(pmatrix(2, 0), pmatrix(2, 1), pmatrix(2, 2));

        //旋转变换矩阵
        Matx33d Trans0 = g_initMat33;
        Matx33d Trans1 = g_initMat33;
        Matx33d Trans2 = g_initMat33;
        Matx33d Trans = g_initMat33;

        //旋转角度
        double dRot[3];
        double dcos0, dcos1, dcos2, dsin0, dsin1, dsin2;
        for (int i = 0; i < 3; i++) {
            dRot[i] = D2R(dAngle[i]);
        }
        dcos0 = cos(dRot[0]);
        dcos1 = cos(dRot[1]);
        dcos2 = cos(dRot[2]);
        dsin0 = sin(dRot[0]);
        dsin1 = sin(dRot[1]);
        dsin2 = sin(dRot[2]);

        Trans0(0, 0) = dcos0 + pow(vec0(0, 0), 2) * (1 - dcos0);
        Trans0(0, 1) = vec0(0, 0) * vec0(0, 1) * (1 - dcos0) - vec0(0, 2) * dsin0;
        Trans0(0, 2) = vec0(0, 0) * vec0(0, 2) * (1 - dcos0) + vec0(0, 1) * dsin0;
        Trans0(1, 0) = vec0(0, 1) * vec0(0, 0) * (1 - dcos0) + vec0(0, 2) * dsin0;
        Trans0(1, 1) = dcos0 + pow(vec0(0, 1), 2) * (1 - dcos0);
        Trans0(1, 2) = vec0(0, 1) * vec0(0, 2) * (1 - dcos0) - vec0(0, 0) * dsin0;
        Trans0(2, 0) = vec0(0, 2) * vec0(0, 0) * (1 - dcos0) - vec0(0, 1) * dsin0;
        Trans0(2, 1) = vec0(0, 2) * vec0(0, 1) * (1 - dcos0) + vec0(0, 0) * dsin0;
        Trans0(2, 2) = dcos0 + pow(vec0(0, 2), 2) * (1 - dcos0);

        Trans1(0, 0) = dcos1 + pow(vec1(0, 0), 2) * (1 - dcos1);
        Trans1(0, 1) = vec1(0, 0) * vec1(0, 1) * (1 - dcos1) - vec1(0, 2) * dsin1;
        Trans1(0, 2) = vec1(0, 0) * vec1(0, 2) * (1 - dcos1) + vec1(0, 1) * dsin1;
        Trans1(1, 0) = vec1(0, 1) * vec1(0, 0) * (1 - dcos1) + vec1(0, 2) * dsin1;
        Trans1(1, 1) = dcos1 + pow(vec1(0, 1), 2) * (1 - dcos1);
        Trans1(1, 2) = vec1(0, 1) * vec1(0, 2) * (1 - dcos1) - vec1(0, 0) * dsin1;
        Trans1(2, 0) = vec1(0, 2) * vec1(0, 0) * (1 - dcos1) - vec1(0, 1) * dsin1;
        Trans1(2, 1) = vec1(0, 2) * vec1(0, 1) * (1 - dcos1) + vec1(0, 0) * dsin1;
        Trans1(2, 2) = dcos1 + pow(vec1(0, 2), 2) * (1 - dcos1);

        Trans2(0, 0) = dcos2 + pow(vec2(0, 0), 2) * (1 - dcos2);
        Trans2(0, 1) = vec2(0, 0) * vec2(0, 1) * (1 - dcos2) - vec2(0, 2) * dsin2;
        Trans2(0, 2) = vec2(0, 0) * vec2(0, 2) * (1 - dcos2) + vec2(0, 1) * dsin2;
        Trans2(1, 0) = vec2(0, 1) * vec2(0, 0) * (1 - dcos2) + vec2(0, 2) * dsin2;
        Trans2(1, 1) = dcos2 + pow(vec2(0, 1), 2) * (1 - dcos2);
        Trans2(1, 2) = vec2(0, 1) * vec2(0, 2) * (1 - dcos2) - vec2(0, 0) * dsin2;
        Trans2(2, 0) = vec2(0, 2) * vec2(0, 0) * (1 - dcos2) - vec2(0, 1) * dsin2;
        Trans2(2, 1) = vec2(0, 2) * vec2(0, 1) * (1 - dcos2) + vec2(0, 0) * dsin2;
        Trans2(2, 2) = dcos2 + pow(vec2(0, 2), 2) * (1 - dcos2);


        MatMul(Trans, Trans2.t(), Trans);
        MatMul(Trans, Trans1.t(), Trans);
        MatMul(Trans, Trans0.t(), Trans);

        return Trans;