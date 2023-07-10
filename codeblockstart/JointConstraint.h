#pragma once

#include "HigherOrderBeamLib.h"
#include "Edge.h"
#include "DOF.h"
#include "PSI.h"
#include "LinearAlgebra.h"

tuple<vector<VectorXi>, vector<VectorXi>, vector<VectorXd>> jntmat(vector<MatrixXd> csc, vector<MatrixXi> cscc, vector<vector<MatrixXd>> coef_s, vector<vector<MatrixXd>> coef_n, vector<vector<MatrixXd>> coef_z, MatrixXd pnt, MatrixXi pnt_con, MatrixXd xdir, vector<VectorXd> zcoord, VectorXi Nel, ArrayXi csnum, VectorXi Nmode, ArrayXi NW1, MatrixXd besz)
{
    VectorXd zcoord_end = zcoord[zcoord.size() - 1], zcoord1, zcoord2, zco1, zco2, Sall1;
    VectorXi enddof = noddof3(Nel, Nmode, zcoord.size(), zcoord_end.size()), dof1i, dof2i, dofi, Si1, Sj1;
    int smsize = enddof(Eigen::last), cs1, c1, b1, e1, b2, e2, cs2, c2, nod1, nod2, count1, count2 = 0;
    MatrixXd pc1, csc1, S, pc2, csc2, D1, R1, D2, R2, S0, Si, SS;//, epcoord1(besz.rows(), 3), jpcoord1(besz.rows(), 3), epcoord2(besz.rows(), 3), jpcoord2(besz.rows(), 3);
    //vector<MatrixXd> epcoord, jpcoord;
    double s1, zc1, yc1, xc1, s2, zc2, yc2, xc2, r1, r2;
    Vector3d zdir1, xdir1, ydir1, zdir2, xdir2, ydir2;
    MatrixXi cscc1, cscc2;
    vector<MatrixXd> coef_s1, coef_z1, coef_n1, coef_s2, coef_z2, coef_n2;
    ArrayXd alpha1, L1, alpha2, L2;
    Vector3d cc1, CC1, cc2, CC2, JP, cp1, cp2;
    Matrix3d T1, T2;
    vector<int> sir, sit, sjt;
    vector<VectorXi> dof, Sii, Sjj;
    vector<VectorXd> Sall;
    vector<double> sir2, sall;
    vector<vector<double>> sir4;
    vector<vector<int>> sir3;
    SparseMatrix<double> SM;
    typedef Triplet<double> Trip;
    vector<Trip> trp;
    Edge edge;


    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    for (int i = 0; i < besz.rows(); i++)
    {
        // Data - Beam 1 - Basic Data
        count1 = 0;
        b1 = besz(i, 0);
        e1 = besz(i, 1);
        s1 = besz(i, 2);
        zc1 = besz(i, 3);
        pc1 = pnt(pnt_con(b1, all), all);

        // Beam Direction

        zdir1 = normalize(pc1(1, all) - pc1(0, all));
        xdir1 = normalize(xdir(b1, all));
        ydir1 = normalize(zdir1.cross(xdir1));

        // Cross Section Data

        cs1 = csnum(b1);
        csc1 = csc[cs1];
        cscc1 = cscc[cs1];
        coef_s1 = coef_s[cs1];
        coef_n1 = coef_n[cs1];
        coef_z1 = coef_z[cs1];

        // Corner Coordinate

        alpha1 = edge.Angle(csc1, cscc1);
        L1 = edge.Length(csc1, cscc1);
        c1 = cscc1(e1, 0);
        xc1 = csc1(c1, 0) + s1 * cos(alpha1(e1));
        yc1 = csc1(c1, 1) + s1 * sin(alpha1(e1));
        cc1 << xc1, yc1, zc1;
        CC1 = axistrans(cc1, pc1(0, all), zdir1, xdir1, 0);

        // Joint Dof

        zcoord1 = zcoord[b1];
        zco1 = (zcoord1 - zc1 * VectorXd::Ones(zcoord1.size())).cwiseAbs();
        for (int j = 0; j < zcoord1.size(); j++)
        {
            if (abs(zcoord1(j) - zc1) == zco1.minCoeff())
                nod1 = j;
        }
        dof1i = noddof3(Nel, Nmode, b1 + 1, nod1 + 1);
        dof1i -= VectorXi::Ones(dof1i.size());

        // Transformation Matrix

        T1 << xdir1.transpose(),
            ydir1.transpose(),
            zdir1.transpose();
        T1 = T1.inverse().eval();

        // Data - Beam 2 - Basic Data

        b2 = besz(i, 4);
        e2 = besz(i, 5);
        s2 = besz(i, 6);
        zc2 = besz(i, 7);
        pc2 = pnt(pnt_con(b2, all), all);

        // Beam Direction

        zdir2 = normalize(pc2(1, all) - pc2(0, all));
        xdir2 = normalize(xdir(b2, all));
        ydir2 = normalize(zdir2.cross(xdir2));

        // Cross Section Data

        cs2 = csnum(b2);
        csc2 = csc[cs2];
        cscc2 = cscc[cs2];
        coef_s2 = coef_s[cs2];
        coef_n2 = coef_n[cs2];
        coef_z2 = coef_z[cs2];

        // Corner Coordinate

        alpha2 = edge.Angle(csc2, cscc2);
        L2 = edge.Length(csc2, cscc2);
        c2 = cscc2(e2, 0);
        xc2 = csc2(c2, 0) + s2 * cos(alpha2(e2));
        yc2 = csc2(c2, 1) + s2 * sin(alpha2(e2));
        cc2 << xc2, yc2, zc2;
        CC2 = axistrans(cc2, pc2(0, all), zdir2, xdir2, 0);

        // Joint Dof

        zcoord2 = zcoord[b2];
        zco2 = (zcoord2 - zc2 * VectorXd::Ones(zcoord2.size())).cwiseAbs();
        for (int j = 0; j < zcoord2.size(); j++)
        {
            if (abs(zcoord2(j) - zc2) == zco2.minCoeff())
                nod2 = j;
        }
        dof2i = noddof3(Nel, Nmode, b2 + 1, nod2 + 1);
        dof2i -= VectorXi::Ones(dof2i.size());

        dofi.setZero(dof1i.size() + dof2i.size());
        dofi.block(0, 0, dof1i.size(), 1) = dof1i;
        dofi.block(dof1i.size(), 0, dof2i.size(), 1) = dof2i;
        dof.push_back(dofi);

        // Transformation Matrix

        T2 << xdir2.transpose(),
            ydir2.transpose(),
            zdir2.transpose();
        T2 = T2.inverse().eval();

        // Joint Condition

        tie(r1, r2) = intersection(CC1, zdir1, CC2, zdir2);
        if (r1 == INT_MAX)
            JP = (CC1 + CC2) / 2;
        else
            JP = CC1 + zdir1 * r1;

        cp1 = T1.inverse() * (JP - CC1);
        cp2 = T2.inverse() * (JP - CC2);
        tie(D1, R1) = jntdisp(csc1, cscc1, coef_s1, coef_n1, coef_z1, e1, s1, cp1, NW1(cs1));
        tie(D2, R2) = jntdisp(csc2, cscc2, coef_s2, coef_n2, coef_z2, e2, s2, cp2, NW1(cs2));

        S0.setZero(2 * T1.rows(), D1.cols() + D2.cols());
        S0.block(0, 0, T1.rows(), D1.cols()) = T1 * D1;
        S0.block(0, D1.cols(), T2.rows(), D2.cols()) = -T2 * D2;
        S0.block(T1.rows(), 0, T1.rows(), R1.cols()) = T1 * R1;
        S0.block(T1.rows(), D1.cols(), T2.rows(), R2.cols()) = -T2 * R2;

        Si.setZero(S0.rows(), S0.cols());
        S.setZero(S0.rows(), S0.cols());
        for (int j = 0;j < S0.rows(); j++)
        {
            Si.block(0, 0, count1, S.cols()) = S.block(0, 0, count1, S.cols());
            Si.block(count1, 0, 1, S0.cols()) = S0(j, all);
            FullPivLU<Matrix<double, -1, -1>> Lu1(Si);
            FullPivLU<Matrix<double, -1, -1>> Lu2(S);
            if (Lu1.rank() > Lu2.rank())
            {
                S.block(0, 0, count1 + 1, S0.cols()) = Si.block(0, 0, count1 + 1, S0.cols());
                count1++;
            }
        }

        SS.setZero(S.rows(), S.cols());
        for (int i = 0; i < S.rows(); i++)
        {
            SS(i, all) = S(i, all) / S(i, all).cwiseAbs().maxCoeff();
        }

        sir3.clear();
        sir4.clear();
        for (int j = 0; j < SS.cols(); j++)
        {
            sir.clear();
            sir2.clear();
            for (int k = 0; k < SS.rows(); k++)
            {
                if (SS(k, j) != 0)
                {
                    sir.push_back(count2 + k);
                    sir2.push_back(SS(k, j));
                }
            }
            sir3.push_back(sir);
            sir4.push_back(sir2);
        }

        sit.clear(), sjt.clear(), sall.clear();
        for (int j = 0; j < dofi.size(); j++)
        {
            for (int k = 0; k < sir3[j].size(); k++)
            {
                sit.push_back(sir3[j][k]);
                sjt.push_back(dofi(j));
                sall.push_back(sir4[j][k]);
                trp.push_back(Trip(sir3[j][k], dofi(j), sir4[j][k]));
            }
        }
        count2 += SS.rows();
        Si1.setZero(sit.size());
        Sj1.setZero(sjt.size());
        Sall1.setZero(sall.size());

        for (int k = 0; k < sit.size(); k++)
        {
            Si1(k) = sit[k];
            Sj1(k) = sjt[k];
            Sall1(k) = sall[k];
        }
        Sii.push_back(Si1);
        Sjj.push_back(Sj1);
        Sall.push_back(Sall1);

        /*// End Point Coordinate
        epcoord1(i, all) = CC1;
        epcoord2(i, all) = CC2;
        // Joint Point Coordinate
        jpcoord1(i, all) = JP;
        jpcoord2(i, all) = JP;*/
    }
    SM.resize(count2, smsize);
    SM.setFromTriplets(trp.begin(), trp.end());
    SM.makeCompressed();

    //epcoord.push_back(epcoord1);
    //epcoord.push_back(epcoord2);
    //jpcoord.push_back(jpcoord1);
    //jpcoord.push_back(jpcoord2);

    return make_tuple(Sii, Sjj, Sall);
}

tuple<MatrixXd, MatrixXd> jntdisp(MatrixXd csc, MatrixXi cscc, vector<MatrixXd> coef_s, vector<MatrixXd> coef_n, vector<MatrixXd> coef_z, int e11, double sc1, Vector3d r, int Nw1)
{
    double tol = 1e-8, a1, a2;
    int Ni = coef_s[0].cols(), No = coef_z[0].cols(), Nmode = Ni + No, c, e1, e2;
    ArrayXd alpha, Le;
    vector<Vector2i> es;
    ArrayXi e;
    ArrayXd s, s1(1), s2(1), sc;
    MatrixXd D, R, F, dF, A1F, A2F, ATF, dF1, dF2, dFz, R0, Rarm, D0;
    Matrix3d A1, A2, A3, A4, T1;
    Matrix2d T2;
    Vector2i ref, tt, m_r;
    Vector3d zz1, zz2;
    VectorXi m_w, tt2;
    Edge edge;

    alpha = edge.Angle(csc, cscc);
    Le = edge.Length(csc, cscc);

    // At a Corner?

    if (sc1 < Le(e11) * tol || abs(Le(e11) - sc1) < Le(e11) * tol)
    {
        if (sc1 < Le(e11) * tol)
            c = cscc(e11, 0);
        else
            c = cscc(e11, 1);

        for (int j = 0; j < cscc.cols(); j++)
        {
            for (int i = 0; i < cscc.rows(); i++)
            {
                if (cscc(i, j) == c)
                    es.push_back({ i,j });
            }
        }
        e.setZero(es.size());
        s.setZero(es.size());
        for (int i = 0; i < es.size(); i++)
        {
            e(i) = es[i](0);
            s(i) = es[i](1);
        }
        sc = Le(e).cwiseProduct(s);
    }
    else
    {
        sc.setZero(1);
        e.setZero(1);
        sc = sc1;
        e = e11;
    }

    A1 << 0, 0, 0,
        0, 0, 1,
        0, -1, 0;
    A2 << 0, 1, 0,
        0, 0, 0,
        0, 0, 0;
    A3 << 0, r(2), 0,
        -r(2), 0, 0,
        0, 0, 0;
    A4 << 0, 0, -r(1),
        0, 0, r(0),
        r(1), -r(0), 0;

    e1 = e(0);
    T1 << cos(alpha(e1)), sin(alpha(e1)), 0,
        sin(alpha(e1)), -cos(alpha(e1)), 0,
        0, 0, 1;

    // Rotation
    if (e.size() == 1)
    {
        F = PSI(coef_s, coef_n, coef_z, e(0), sc, 0);
        dF = PSI(coef_s, coef_n, coef_z, e(0), sc, 1);
        A1F = A1 * dF;
        A2F = A2 * F;
        ATF.setZero(A1F.rows(), A1F.cols() + A2F.cols());
        ATF.block(0, 0, A1F.rows(), A1F.cols()) = A1F;
        ATF.block(0, A1F.cols(), A2F.rows(), A2F.cols()) = A2F;
        R = T1 * ATF;
    }
    else if (e.size() >= 2)
    {
        if (abs(sin(alpha(e(0)) - alpha(e(1)))) < 1e-10)
            ref << 0, 2; // Reference Range
        else
            ref << 0, 1; // Reference Range

        // Reference Coordinate
        e1 = e(ref(0));
        s1 = sc(ref(0));
        a1 = alpha(e1);
        dF1 = PSI(coef_s, coef_n, coef_z, e1, s1, 1);
        e2 = e(ref(1));
        s2 = sc(ref(1));
        a2 = alpha(e2);
        dF2 = PSI(coef_s, coef_n, coef_z, e2, s2, 1);

        T2 << cos(a2) / sin(a1 - a2), -cos(a1) / sin(a1 - a2),
            sin(a2) / sin(a1 - a2), -sin(a1) / sin(a1 - a2);
        zz1 << 0, 0, 1;
        dFz.setZero(2, dF1.cols());
        dFz.block(0, 0, 1, dF1.cols()) = zz1.transpose() * dF1;
        dFz.block(1, 0, 1, dF2.cols()) = zz1.transpose() * dF2;

        zz2 << 0, -1, 0;
        R0.setZero(3, dFz.cols());
        R0.block(0, 0, 2, dFz.cols()) = T2 * dFz;
        R0.block(2, 0, 1, dFz.cols()) = zz2.transpose() * dF1;
        R.setZero(R0.rows(), R0.cols() + Nmode);
        R.block(0, 0, R0.rows(), R0.cols()) = R0;
    }

    // ------------ Arm Rotation ------------ \\
    // Angle By Bending Rotation

    tt << 2, 3;
    m_r = tt + Ni * Vector2i::Ones();
    Rarm = bendrot(csc, cscc, coef_s, coef_n, coef_z, m_r);

    // 1st set Warping Modes

    if (Nw1 + 1 > 3)
    {
        tt2.setLinSpaced(Nw1 - 2, 3, Nw1 + 1);
        m_w = tt2 + (Ni - 1) * VectorXi::Ones(tt2.size());
        Rarm(all, m_w) = R(all, m_w);
    }
    else if (Nw1 + 1 == 3) {
        Rarm(all, 3 + Ni) = R(all, 3 + Ni);
    }

    // ------------ Displacement ------------ \\
    // Displacement at Beam End

    e1 = e(0);
    s1 = sc(0);
    F = PSI(coef_s, coef_n, coef_z, e1, s1, 0);
    D0.setZero(T1.rows(), F.cols() + Nmode);
    D0.block(0, 0, T1.rows(), F.cols()) = T1 * F;

    // Displacement at Joint
    D = D0 + A3 * Rarm + A4 * R;

    return make_tuple(D, R);
}

MatrixXd bendrot(MatrixXd csc, MatrixXi cscc, vector<MatrixXd> coef_s, vector<MatrixXd> coef_n, vector<MatrixXd> coef_z, Vector2i m_r)
{
    int Ni = coef_s[0].cols(), No = coef_z[0].cols(), Nmode = Ni + No;
    ArrayXd alpha, Le;
    MatrixXd R, df1, dF1, df2, dF2, dFz, R0;
    int i = 1, nci, c, e1, e2;
    double a1, a2;
    vector<int> nclf;
    vector<Vector2i> es;
    VectorXi e;
    ArrayXd s, s1(1), s2(1);
    VectorXd sc;
    Vector2i ref, m;
    Matrix2d T;
    Edge edge;

    alpha = edge.Angle(csc, cscc);
    Le = edge.Length(csc, cscc);

    m = m_r - 1 * Vector2i::Ones();

    while (i)
    {
        for (int j = 0; j < cscc.size(); j++)
        {
            if (cscc(j) == i - 1)
                nclf.push_back(j);
        }
        nci = nclf.size();
        if (nci >= 2)
        {
            c = i - 1;
            i = 0;
        }
        else
            i += 1;
    }

    // Continuity Data

    for (int j = 0; j < cscc.cols(); j++)
    {
        for (int k = 0; k < cscc.rows(); k++)
        {
            if (cscc(k, j) == c)
                es.push_back({ k,j });
        }
    }
    e.setZero(es.size());
    s.setZero(es.size());
    sc.setZero(e.size());
    for (int j = 0; j < es.size(); j++)
    {
        e(j) = es[j](0);
        s(j) = es[j](1);
    }
    for (int j = 0; j < e.size(); j++)
    {
        if (s(j) == 0)
            sc(j) = 0;
        else
            sc(j) = Le(e(j));
    }

    if (abs(sin(alpha(e(0)) - alpha(e(1)))) < 1e-10)
        ref << 0, 2; // Reference Range
    else
        ref << 0, 1; // Reference Range

    // Reference Coordinate

    e1 = e(ref(0));
    s1 = sc(ref(0));
    a1 = alpha(e1);
    df1 = PSI(coef_s, coef_n, coef_z, e1, s1, 1);
    dF1.setZero(1, Nmode);
    dF1(0, m) = df1(2, m);

    e2 = e(ref(1));
    s2 = sc(ref(1));
    a2 = alpha(e2);
    df2 = PSI(coef_s, coef_n, coef_z, e2, s2, 1);
    dF2.setZero(1, Nmode);
    dF2(0, m) = df2(2, m);

    T << cos(a2) / sin(a1 - a2), -cos(a1) / sin(a1 - a2),
        sin(a2) / sin(a1 - a2), -sin(a1) / sin(a1 - a2);

    dFz.setZero(dF1.rows() + dF2.rows(), dF1.cols());
    dFz.block(0, 0, dF1.rows(), dF1.cols()) = dF1;
    dFz.block(dF1.rows(), 0, dF1.rows(), dF1.cols()) = dF2;

    R0.setZero(T.rows() + 1, dFz.cols());
    R0.block(0, 0, T.rows(), dFz.cols()) = T * dFz;

    R.setZero(R0.rows(), R0.cols() + Nmode);
    R.block(0, 0, R0.rows(), R0.cols()) = R0;

    return R;
}