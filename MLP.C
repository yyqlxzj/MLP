TMatrixD MLP(double posinX_Z,double posoutX_Z,double posinY_Z,double posoutY_Z, double xkin, double xkout, double ykin, double ykout, double u1) {
    //need to set parameters by hand:u0,u2,X0,kin,kout,the defined domain of the function
    //unit:cm

    double para1=0.0, para2=0.0;//积分前的参数
    double u0 = 9.50, u2 = 19.50;//入射、出射的z
    double X0 = 36.10, E = 13.6 * 13.6;//X0为辐射长度、E和下式中的0.038为经验常数
    para1 = E * (1.0 + 0.038 * log((u1 - u0) / X0)) * (1.0 + 0.038 * log((u1 - u0) / X0));
    para2 = E * (1.0 + 0.038 * log((u2 - u1) / X0)) * (1.0 + 0.038 * log((u2 - u1) / X0));

    double func1(double* t, double* p);
    double func2(double* t, double* p);
    double func3(double* t, double* p);

    TF1 f1("func1", func1, 9.0, 20.50, 3);
    TF1 f2("func2", func2, 9.0, 20.50, 2);
    TF1 f3("func3", func3, 9.0, 20.50, 3);
    TF1 f4("func4", func1, 9.0, 20.50, 3);
    TF1 f5("func5", func2, 9.0, 20.50, 2);
    TF1 f6("func6", func3, 9.0, 20.50, 3);

    f1.SetParameters(para1, u1, X0);
    f2.SetParameters(para1, X0);
    f3.SetParameters(para1, u1, X0);
    f4.SetParameters(para2, u2, X0);
    f5.SetParameters(para2, X0);
    f6.SetParameters(para2, u2, X0);


    auto gamat1 = f1.Integral(u0, u1);
    auto gamac1 = f2.Integral(u0, u1);
    auto gamatc1 = f3.Integral(u0, u1);
    auto gamat2 = f4.Integral(u1, u2);
    auto gamac2 = f5.Integral(u1, u2);
    auto gamatc2 = f6.Integral(u1, u2);//位移、角度和（角度+位移）的协方差

    TMatrixD sigema1(2, 2), sigema2(2, 2), R0(2, 2), R1(2, 2),  matrixI(2, 2), matrixi(2, 1), Xmlp(2, 1), Ymlp(2, 1);//sigema为散射矩阵
    TMatrixD Y0(2, 1), Y2(2, 1), X_0(2, 1), X_2(2, 1);

    Y0(0, 0) = posinY_Z, Y0(1, 0) = atan(ykin), Y2(0, 0) = posoutY_Z, Y2(1, 0) = atan(ykout);
    X_0(0, 0) = posinX_Z, X_0(1, 0) = atan(xkin), X_2(0, 0) = posoutX_Z, X_2(1, 0) = atan(xkout);

    sigema1(0, 0) = gamat1, sigema1(0, 1) = gamatc1, sigema1(1, 0) = gamatc1, sigema1(1, 1) = gamac1;
    sigema2(0, 0) = gamat2, sigema2(0, 1) = gamatc2, sigema1(1, 0) = gamatc2, sigema2(1, 1) = gamac2;
    R0(0, 0) = 1.0, R0(0, 1) = u1 - u0, R0(1, 0) = 0.0, R0(1, 1) = 1.0;
    R1(0, 0) = 1.0, R1(0, 1) = u2 - u1, R1(1, 0) = 0.0, R1(1, 1) = 1.0;

    matrixI(0, 0) = 1.0, matrixI(1, 1) = 1.0, matrixi(0, 0) = 1.0, matrixi(1, 0) = 1.0;
    Xmlp = matrixi, Ymlp = matrixi;//将矩阵初始化，可删

    TMatrixD sigema1_inv = sigema1, sigema2_inv = sigema2, R1_T = R1, part1 = matrixI, partY2 = matrixi, partX2 = matrixi;//inv：求逆矩阵，T：求转置

    sigema1_inv.Invert(), sigema2_inv.Invert(), R1_T.T();
    part1 = sigema1_inv + R1_T * sigema2_inv * R1;
    part1.Invert();

    partY2 = sigema1_inv * R0 * Y0 + R1_T * sigema2_inv * Y2;
    partX2 = sigema1_inv * R0 * X_0 + R1_T * sigema2_inv * X_2;

    Ymlp = part1 * partY2;
    Xmlp = part1 * partX2;

    TMatrixD mlp(2,2);
    mlp(0,0) = Xmlp(0,0);
    mlp(1,0) = Xmlp(1,0);
    mlp(0,1) = Ymlp(0,0);
    mlp(1,1) = Ymlp(1,0);

    return mlp;
}

double func1(double* t, double* p) {
    double x = t[0];
    double a = p[0];
    double b = p[1];
    double c = p[2];

    return a * (b - x) * (b - x) / ((7.457e-6 + 4.548e-7 * x - 5.777e-8 * x * x + 1.301e-8 * x * x * x - 9.228e-10 * x * x * x * x + 2.687e-11 * x * x * x * x * x) * c);
}

double func2(double* t, double* p) {
    double x = t[0];
    double a = p[0];
    double c = p[1];

    return a / ((7.457e-6 + 4.548e-7 * x - 5.777e-8 * x * x + 1.301e-8 * x * x * x - 9.228e-10 * x * x * x * x + 2.687e-11 * x * x * x * x * x) * c);
}

double func3(double* t, double* p) {
    double x = t[0];
    double a = p[0];
    double b = p[1];
    double c = p[2];

    return a * (b - x) / ((7.457e-6 + 4.548e-7 * x - 5.777e-8 * x * x + 1.301e-8 * x * x * x - 9.228e-10 * x * x * x * x + 2.687e-11 * x * x * x * x * x) * c);
}
