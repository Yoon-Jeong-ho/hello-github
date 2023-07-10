// Generates Cross-sectional Modes

#include "HigherOrderBeamLib.h"
#include "Edge.h"
#include "Center.h"
#include "Calculus.h"
#include "Inertia.h"
#include "Continuity.h"
#include "Probtune.h"
#include "ScaleFactor.h"
#include "nDirectional.h"

tuple<vector<MatrixXd>, vector<MatrixXd>, vector<MatrixXd>, MatrixXi> Modegen(Matrix<cpp_dec_float_50, -1, -1> csc, MatrixXi cscc, int Mset, MatrixXi settab)
{
    Matrix<double50, -1, -1> cR, cT;
    Array<double50, -1, 1> alpha, Le;
    double50 betaR, betaT;
    ArrayXi stw, stx, ste;
    Edge edge;
    Center center;
    Index NE = cscc.rows();
    Tensor<double, 3> pp1, pp2, pp3;

    alpha = edge.Anglemp(csc, cscc);
    Le = edge.Lengthmp(csc, cscc);
    cR = center.Rot.Coordinates(csc, cscc, alpha, Le);
    betaR = center.Rot.Angle(csc, cscc, alpha, Le, cR);
    betaT = center.Tor.Angle(csc, cscc, alpha, Le, betaR);
    cT = center.Tor.Coordinates(csc, cscc, alpha, Le, betaT);

    vector<Matrix<double50, -1, -1>> coef_s, coef_n, coef_z, coef_u, coef_sWb, coef_nWb;
    coef_sWb.resize(NE), coef_nWb.resize(NE);

    tie(coef_s, coef_n, coef_z, coef_u) = Rigid(NE, csc, cscc, alpha, betaT, cT, cR, betaR);
    settab(0, all) << 1, 3, 0;

    for (int i = 0; i < Mset; i++)
    {
        if (i == 0)
            tie(coef_s, coef_n, coef_z, stw) = LW_InXD(csc, cscc, coef_s, coef_n, coef_z, coef_u, Le);
        else
            std::tie(coef_z, stw) = NLW(csc, cscc, coef_s, coef_z, coef_u, Le);

        tie(coef_s, stx) = ExtDist(csc, cscc, coef_s, coef_z, coef_u, Le);
        coef_n = nDirectional(csc, cscc, coef_s, coef_n);

        tie(coef_nWb, ste) = WallBend(csc, cscc, coef_nWb, coef_n, coef_u, Le, i);

        settab(i + 1, all) = stw.transpose() + stx.transpose() + ste.transpose();
    }

    Index NMs = coef_s[0].cols(), NMswb = coef_nWb[0].cols(), NC = coef_s[0].rows(), NMz = coef_z[0].cols();
    vector<MatrixXd> Cs, Cn, Cz; Cs.resize(NE), Cn.resize(NE), Cz.resize(NE);

    for (int i = 0; i < NE; i++)
    {
        coef_s[i].conservativeResize(NC, NMs + NMswb);
        coef_n[i].conservativeResize(NC, NMs + NMswb);

        coef_s[i].block(0, NMs, NC, NMswb).setZero();
        coef_n[i].block(0, NMs, NC, NMswb) = coef_nWb[i];
        Cs[i].setZero(NC, NMs + NMswb);
        Cn[i].setZero(NC, NMs + NMswb);
        Cz[i].setZero(NC, NMz);
    }

    for (int i = 0; i < NE; i++)
    {
        for (int j = 0; j < NC; j++)
        {
            for (int k = 0; k < NMs + NMswb; k++)
            {
                Cs[i](j, k) = static_cast<double>(coef_s[i](j, k));
                Cn[i](j, k) = static_cast<double>(coef_n[i](j, k));
            }
            for (int k = 0; k < NMz; k++)
                Cz[i](j, k) = static_cast<double>(coef_z[i](j, k));
        }
    }

    return make_tuple(Cs, Cn, Cz, settab);
}

tuple<vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>> Rigid(Index NE, Matrix<double50, -1, -1> csc, MatrixXi cscc, Array<double50, -1, 1> alpha, double50 betaT, Array<double50, -1, 1> cT, Array<double50, -1, 1> cR, double50 betaR)
{
    //std::cout << "Deriving Rigid Body Modes" << endl << endl;
    int order = 25;
    //Tensor<double50, 3> coef_s(order + 1, 3, NE), coef_n(order + 1, 3, NE), coef_z(order + 1, 3, NE), coef_u(order + 1, NE, NE);
    vector<Matrix<double50, -1, -1>> coef_s, coef_n, coef_z, coef_u;
    coef_s.resize(NE), coef_n.resize(NE), coef_z.resize(NE), coef_u.resize(NE);
    for (int i = 0; i < NE; i++)
    {
        coef_s[i].setZero(order + 1, 3);
        coef_n[i].setZero(order + 1, 3);
        coef_z[i].setZero(order + 1, 3);
        coef_u[i].setZero(order + 1, NE);
    }
    Vector<double50, -1> Xi, Yi;
    Xi = csc(cscc(all, 0), 0), Yi = csc(cscc(all, 0), 1);

    for (int i = 0; i < NE; i++)
    {
        coef_s[i](0, 0) = cos(alpha(i) - betaT);
        coef_s[i](0, 1) = sin(alpha(i) - betaT);
        coef_s[i](0, 2) = (Xi(i) - cT(0)) * sin(alpha(i)) - (Yi(i) - cT(1)) * cos(alpha(i));

        coef_n[i](0, 0) = sin(alpha(i) - betaT);
        coef_n[i](0, 1) = -cos(alpha(i) - betaT);
        coef_n[i](0, 2) = -(Xi(i) - cT(0)) * cos(alpha(i)) - (Yi(i) - cT(1)) * sin(alpha(i));
        coef_n[i](1, 2) = double50(-1);

        coef_z[i](0, 0) = 1;
        coef_z[i](0, 1) = -(Xi(i) - cR(0)) * sin(betaR) + (Yi(i) - cR(1)) * cos(betaR);
        coef_z[i](0, 2) = -(Xi(i) - cR(0)) * cos(betaR) - (Yi(i) - cR(1)) * sin(betaR);
        coef_z[i](1, 1) = sin(alpha(i) - betaR);
        coef_z[i](1, 2) = -cos(alpha(i) - betaR);

        coef_u[i](0, i) = double50(1);
    }

    //std::cout << "Rigid Body Modes Derived" << endl << endl;
    return make_tuple(coef_s, coef_n, coef_z, coef_u);
}

tuple<vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, ArrayXi> LW_InXD(Matrix<double50, -1, -1> csc, MatrixXi cscc, vector<Matrix<double50, -1, -1>> coef_s, vector<Matrix<double50, -1, -1>> coef_n, vector<Matrix<double50, -1, -1>> coef_z, vector<Matrix<double50, -1, -1>> coef_u, Array<double50, -1, 1> Le)
{
    //std::cout << "Deriving Linear Warping and In-Extensional Distortion Modes" << endl << endl;
    Index NE = cscc.rows(), NC = coef_s[0].rows(), NMs = coef_s[0].cols(), NMu = coef_u[0].cols(), NMz = coef_z[0].cols();
    int count1 = 0;

    ArrayXi stw(3); stw.setZero();
    Matrix<double50, -1, -1>  C, eye;
    vector<Matrix<double50, -1, -1>> Phi_W, w1, w2; Phi_W.resize(NE);
    Calculus operate;
    Continuity cont;

    C.setZero(NE, (2 * NMu + NMs)), eye.setIdentity(NE, NE);
    C.block(0, NMu, NE, NE) = eye;

    w1 = operate.Intcoef(coef_s), w2 = operate.Intcoef(coef_u);

    for (int i = 0;i < NE;i++)
    {
        Phi_W[i].setZero(NC, 2 * NMu + NMs);

        Phi_W[i].block(0, 0, NC, NMu) = coef_u[i];
        Phi_W[i].block(0, NMu, NC, NMu) = w2[i];
        Phi_W[i].block(0, 2 * NMu, NC, NMs) = w1[i];
    }

    Matrix<double50, -1, -1> Pw = Inertia(Phi_W, Phi_W, Le);
    Matrix<double50, -1, -1> Qx = Inertia(coef_s, coef_u, Le) * C;
    Matrix<double50, -1, -1> Qw = Inertia(coef_z, Phi_W, Le);

    vector<vector<int>> cont_edg;
    vector<Matrix<double50, -1, -1>> cont_Rxy;
    vector<Array<double50, -1, 1>> cont_sc;
    Matrix<double50, -1, -1> Rw, Rx;

    tie(cont_edg, cont_Rxy, cont_sc) = cont.Edges(csc, cscc, Phi_W, 0);
    Rw = cont.Constraint(csc, cscc, Phi_W, 0);
    Rx = cont.Constraint(csc, cscc, coef_u, 1);

    Matrix<double50, -1, -1> S0, S, Si, x_w;
    vector<int> c, d;

    if (Rx.size() == 0)
    {
        S0.setZero(Qw.rows() + Qx.rows() + Rw.rows(), Qw.cols());
        S0.block(0, 0, Qw.rows(), Qw.cols()) = Qw;
        S0.block(Qw.rows(), 0, Qx.rows(), Qw.cols()) = Qx;
        S0.block(Qw.rows() + Qx.rows(), 0, Rw.rows(), Qw.cols()) = Rw;
    }
    else
    {
        S0.setZero(Qw.rows() + Qx.rows() + Rw.rows() + Rx.rows(), Qw.cols());
        S0.block(0, 0, Qw.rows(), Qw.cols()) = Qw;
        S0.block(Qw.rows(), 0, Qx.rows(), Qw.cols()) = Qx;
        S0.block(Qw.rows() + Qx.rows(), 0, Rw.rows(), Qw.cols()) = Rw;
        S0.block(Qw.rows() + Qx.rows() + Rw.rows(), 0, Rx.rows(), Qw.cols()) = Rx * C;
    }
    Si.setZero(S0.rows(), S0.cols());
    S.setZero(S0.rows(), S0.cols());
    for (int i = 0;i < S0.rows(); i++)
    {
        Si.block(0, 0, count1, S.cols()) = S.block(0, 0, count1, S.cols());
        Si.block(count1, 0, 1, S0.cols()) = S0(i, all);
        FullPivLU<Matrix<double50, -1, -1>> Lu1(Si);
        FullPivLU<Matrix<double50, -1, -1>> Lu2(S);
        if (Lu1.rank() > Lu2.rank())
        {
            S.block(0, 0, count1 + 1, S0.cols()) = Si.block(0, 0, count1 + 1, S0.cols());
            count1++;
        }
    }

    // ######## Linear Warping ######## \\

    std::tie(Phi_W, Pw, x_w, c) = probtune(Phi_W, S, Le, coef_u[0].cols());

    //phi_w = MatrixCast(Phi_W, NC, (2 * NMu + NMs) * NE);

    //std::cout << "No. of LW Modes: " << x_w.cols() << endl << endl;
    stw(0) += x_w.cols();
    Matrix<double50, -1, -1> sfw = ScaleFactor(Le, x_w, Pw);

    Matrix<double50, -1, -1> coef_zi;
    coef_zi.setZero(NC, x_w.cols() * NE);
    for (int i = 0;i < NE;i++)
    {
        coef_z[i].conservativeResize(NC, NMz + x_w.cols());
        coef_z[i].block(0, NMz, NC, x_w.cols()) = Phi_W[i] * x_w * sfw;
    }

    // ######## In-Extensional Distortion ######## \\

    Array<double50, -1, 1> wi, di;
    for (int i = 0; i < x_w.cols(); i++)
    {
        wi = x_w(all, i);
        di = x_w(c, i);
        if (wi.abs().maxCoeff() < di.abs().maxCoeff() * 1e5)
            d.push_back(i);
    }

    if (d.size() > 0)
    {
        Matrix<double50, -1, -1> rev_ed1 = x_w(c, d);
        Matrix<double50, -1, -1> Ped1 = Inertia(coef_u, coef_u, Le);
        Matrix<double50, -1, -1> sfed1 = ScaleFactor(Le, rev_ed1, Ped1);
        vector<Matrix<double50, -1, -1>> coef_si; coef_si.resize(NE);
        for (int i = 0;i < NE;i++)
        {
            coef_si[i].setZero(NC, x_w.cols());
            coef_si[i] = coef_u[i] * rev_ed1 * sfed1;
        }

        Matrix<double50, -1, -1> Ped2 = Inertia(coef_si, coef_si, Le);
        Matrix<double50, -1, -1> Semp;
        Matrix<double50, -1, -1> evec_ed = eig(Ped2, Semp);

        Matrix<double50, -1, -1> sfed2 = ScaleFactor(Le, evec_ed, Ped2);

        for (int i = 0;i < NE;i++)
        {
            coef_s[i].conservativeResize(NC, NMs + evec_ed.cols());
            coef_s[i].block(0, NMs, NC, evec_ed.cols()) = coef_si[i] * evec_ed * sfed2;
        }

        //std::cout << "No. of InExtDist Modes: " << evec_ed.cols() << endl << endl;
        stw(1) = evec_ed.cols();

        coef_n = nDirectional(csc, cscc, coef_s, coef_n);
    }
    /*else
        std::cout << "No. of InExtDist Modes: " << 0 << endl;

    std::cout << "Linear Warping and In-Extensional Distortion Modes Derived" << endl << endl;*/

    return make_tuple(coef_s, coef_n, coef_z, stw);
}

tuple<vector<Matrix<double50, -1, -1>>, ArrayXi> ExtDist(Matrix<double50, -1, -1> csc, MatrixXi cscc, vector<Matrix<double50, -1, -1>> coef_s, vector<Matrix<double50, -1, -1>> coef_z, vector<Matrix<double50, -1, -1>> coef_u, Array<double50, -1, 1> Le)
{
    //std::cout << "Deriving Extensional Distortion Modes" << endl << endl;
    int NE = Le.rows();
    Index NMz = coef_z[0].cols(), NMs = coef_s[0].cols(), NMu = coef_u[0].cols(), NC = coef_s[0].rows();
    int count1 = 0;
    Matrix<double50, -1, -1> Rx;
    ArrayXi stx(3); stx.setZero();
    Calculus operate;
    Continuity cont;

    vector<Matrix<double50, -1, -1>> w1 = operate.Intcoef(coef_z), Phi_X;
    Phi_X.resize(NE);
    vector<int> c;

    for (int i = 0;i < NE;i++)
    {
        Phi_X[i].setZero(NC, NMz + NMu);

        Phi_X[i].block(0, NMu, NC, NMz) = w1[i];
        Phi_X[i].block(0, 0, NC, NMu) = coef_u[i];
    }

    Matrix<double50, -1, -1> Px = Inertia(Phi_X, Phi_X, Le);
    Matrix<double50, -1, -1> Qx = Inertia(coef_s, Phi_X, Le);
    Rx = cont.Constraint(csc, cscc, Phi_X, 1);

    Matrix<double50, -1, -1> S0, Si, S;
    S0.setZero(Qx.rows() + Rx.rows(), Qx.cols());
    S0.block(0, 0, Qx.rows(), Qx.cols()) = Qx;
    S0.block(Qx.rows(), 0, Rx.rows(), Rx.cols()) = Rx;

    Si.setZero(S0.rows(), S0.cols());
    S.setZero(S0.rows(), S0.cols());
    for (int i = 0;i < S0.rows(); i++)
    {
        Si.block(0, 0, count1, S.cols()) = S.block(0, 0, count1, S.cols());
        Si.block(count1, 0, 1, S0.cols()) = S0(i, all);
        FullPivLU<Matrix<double50, -1, -1>> Lu1(Si);
        FullPivLU<Matrix<double50, -1, -1>> Lu2(S);
        if (Lu1.rank() > Lu2.rank())
        {
            S.block(0, 0, count1 + 1, S0.cols()) = Si.block(0, 0, count1 + 1, S0.cols());
            count1++;
        }
    }

    Matrix<double50, -1, -1> evec;
    std::tie(Phi_X, Px, evec, c) = probtune(Phi_X, S, Le, 0);

    //std::cout << "No. of ExtDist Modes: " << evec.cols() << endl;
    stx(1) = evec.cols();

    Matrix<double50, -1, -1> sfxs = ScaleFactor(Le, evec, Px);

    for (int i = 0;i < NE;i++)
    {
        coef_s[i].conservativeResize(NC, NMs + evec.cols());
        coef_s[i].block(0, NMs, NC, evec.cols()) = Phi_X[i] * evec * sfxs;
    }

    //std::cout << "Extensional Distortion Modes Derived" << endl << endl;
    return std::tie(coef_s, stx);
}

tuple<vector<Matrix<double50, -1, -1>>, ArrayXi> WallBend(Matrix<double50, -1, -1> csc, MatrixXi cscc, vector<Matrix<double50, -1, -1>> coef_nWb, vector<Matrix<double50, -1, -1>> coef_n, vector<Matrix<double50, -1, -1>> coef_u, Array<double50, -1, 1> Le, int m)
{
    //std::cout << "Deriving Wall Bending Modes" << endl << endl;
    int NE = cscc.rows();
    Index NC = coef_n[0].rows(), NMnwb = coef_nWb[0].cols(), NMn = coef_n[0].cols(), NMu = coef_u[0].cols(), NCn = coef_nWb[0].rows();
    int NMtb = NMn + NMnwb;
    int count1 = 0;
    vector<Matrix<double50, -1, -1>> coef0; coef0.resize(NE);
    Calculus integ; Continuity cont;
    ArrayXi ste(3); ste.setZero();

    for (int i = 0;i < NE;i++)
    {
        coef0[i].setZero(NC, NMtb);

        coef0[i].block(0, 0, NCn, NMnwb) = coef_nWb[i];
        coef0[i].block(0, NMnwb, NC, NMn) = coef_n[i];
    }
    vector<Matrix<double50, -1, -1>> coef00 = integ.Intcoef(coef0);
    vector<Matrix<double50, -1, -1>> coef01 = integ.Intcoef(coef00);;
    vector<Matrix<double50, -1, -1>> coef1 = integ.Intcoef(coef_u);
    vector<Matrix<double50, -1, -1>> coef2 = integ.Intcoef(coef1);
    vector<Matrix<double50, -1, -1>> coef3 = integ.Intcoef(coef2);
    vector<Matrix<double50, -1, -1>> coef4 = integ.Intcoef(coef3);
    int NMtmod = NMtb + 5 * NMu;
    vector<Matrix<double50, -1, -1>> Phi_N; Phi_N.resize(NE);

    for (int i = 0;i < NE;i++)
    {
        Phi_N[i].setZero(NC, NMtmod);

        Phi_N[i].block(0, 0, NC, NMtb) = coef01[i];
        Phi_N[i].block(0, NMtb, NC, NMu) = coef_u[i];
        Phi_N[i].block(0, NMtb + NMu, NC, NMu) = coef1[i];
        Phi_N[i].block(0, NMtb + 2 * NMu, NC, NMu) = coef2[i];
        Phi_N[i].block(0, NMtb + 3 * NMu, NC, NMu) = coef3[i];
        Phi_N[i].block(0, NMtb + 4 * NMu, NC, NMu) = coef4[i];
    }

    Matrix<double50, -1, -1> Pn = Inertia(Phi_N, Phi_N, Le);
    Matrix<double50, -1, -1> Qn, Q1;
    vector<int> c;
    Q1 = cont.Constraint(csc, cscc, Phi_N, 2);
    if (NCn == 0)
    {
        Qn = Q1;
    }
    else
    {
        Matrix<double50, -1, -1> Q4 = Inertia(coef_nWb, Phi_N, Le);
        if (cscc.rows() == 2)
        {
            vector<Matrix<double50, -1, -1>> Q2; Q2.resize(NE);
            for (int i = 0;i < NE;i++)
            {
                Q2[i] = coef_n[i](all, 2);
            }
            Matrix<double50, -1, -1> Q3 = Inertia(Q2, Phi_N, Le);
            Qn.block(0, 0, Q1.rows(), Q1.cols()) = Q1;
            Qn.block(Q1.rows(), 0, Q3.rows(), Q3.cols()) = Q3;
            Qn.block(Q1.rows() + Q3.rows(), 0, Q4.rows(), Q4.cols()) = Q4;
        }
        else
        {
            Qn.setZero(Q1.rows() + Q4.rows(), Q1.cols());
            Qn.block(0, 0, Q1.rows(), Q1.cols()) = Q1;
            Qn.block(Q1.rows(), 0, Q4.rows(), Q4.cols()) = Q4;
        }
    }
    Matrix<double50, -1, -1> S0, Si, S;
    S0 = Qn;

    Si.setZero(S0.rows(), S0.cols());
    S.setZero(S0.rows(), S0.cols());
    for (int i = 0;i < S0.rows(); i++)
    {
        Si.block(0, 0, count1, S.cols()) = S.block(0, 0, count1, S.cols());
        Si.block(count1, 0, 1, S0.cols()) = S0(i, all);
        FullPivLU<Matrix<double50, -1, -1>> Lu1(Si);
        FullPivLU<Matrix<double50, -1, -1>> Lu2(S);
        if (Lu1.rank() > Lu2.rank())
        {
            S.block(0, 0, count1 + 1, S0.cols()) = Si.block(0, 0, count1 + 1, S0.cols());
            count1++;
        }
    }

    Matrix<double50, -1, -1> evec;
    std::tie(Phi_N, Pn, evec, c) = probtune(Phi_N, S, Le, 0);

    //std::cout << "No. of WB Modes: " << evec.cols() << endl;
    ste(2) = evec.cols();

    Matrix<double50, -1, -1> sfnwb = ScaleFactor(Le, evec, Pn);

    for (int i = 0;i < NE;i++)
    {
        coef_nWb[i].conservativeResize(NC, NMnwb + evec.cols());

        coef_nWb[i].block(0, NMnwb, NC, evec.cols()) = Phi_N[i] * evec * sfnwb;
    }

    //std::cout << "Wall Bending Modes Derived" << endl << endl;
    return make_tuple(coef_nWb, ste);
}

tuple<vector<Matrix<double50, -1, -1>>, ArrayXi> NLW(Matrix<double50, -1, -1> csc, MatrixXi cscc, vector<Matrix<double50, -1, -1>> coef_s, vector<Matrix<double50, -1, -1>> coef_z, vector<Matrix<double50, -1, -1>> coef_u, Array<double50, -1, 1> Le)
{
    //std::cout << "Deriving Non-Linear Warping Modes" << endl << endl;
    int NE = Le.rows();
    Index NC = coef_s[0].rows(), NMs = coef_s[0].cols(), NMz = coef_z[0].cols(), NMu = coef_u[0].cols();
    int NMsu = NMs + NMu, count1 = 0;
    Calculus integ; Continuity cont;

    vector<Matrix<double50, -1, -1>> coef1 = integ.Intcoef(coef_s), Phi_W;
    Matrix<double50, -1, -1> Rw, S0, S{}, Si;
    ArrayXi stw(3); stw.setZero(); Phi_W.resize(NE);

    for (int i = 0;i < NE;i++)
    {
        Phi_W[i].setZero(NC, NMsu);

        Phi_W[i].block(0, 0, NC, NMs) = coef1[i];
        Phi_W[i].block(0, NMs, NC, NMu) = coef_u[i];
    }

    Matrix<double50, -1, -1> Pw = Inertia(Phi_W, Phi_W, Le);
    Matrix<double50, -1, -1> Qw = Inertia(coef_z, Phi_W, Le);

    Rw = cont.Constraint(csc, cscc, Phi_W, 0);

    S0.setZero(Qw.rows() + Rw.rows(), Qw.cols());
    S0.block(0, 0, Qw.rows(), Qw.cols()) = Qw;
    S0.block(Qw.rows(), 0, Rw.rows(), Qw.cols()) = Rw;

    Si.setZero(S0.rows(), S0.cols());
    S.setZero(S0.rows(), S0.cols());
    for (int i = 0;i < S0.rows(); i++)
    {
        Si.block(0, 0, count1, S.cols()) = S.block(0, 0, count1, S.cols());
        Si.block(count1, 0, 1, S0.cols()) = S0(i, all);
        FullPivLU<Matrix<double50, -1, -1>> Lu1(Si);
        FullPivLU<Matrix<double50, -1, -1>> Lu2(S);
        if (Lu1.rank() > Lu2.rank())
        {
            S.block(0, 0, count1 + 1, S0.cols()) = Si.block(0, 0, count1 + 1, S0.cols());
            count1++;
        }
    }

    Matrix<double50, -1, -1> evec;
    vector<int> c;
    std::tie(Phi_W, Pw, evec, c) = probtune(Phi_W, S, Le, 0);

    //std::cout << "No. of NLW Modes: " << evec.cols() << endl;
    stw(0) = evec.cols();

    Matrix<double50, -1, -1> sfw = ScaleFactor(Le, evec, Pw);

    for (int i = 0;i < NE;i++)
    {
        coef_z[i].conservativeResize(NC, NMz + evec.cols());
        coef_z[i].block(0, NMz, NC, evec.cols()) = Phi_W[i] * evec * sfw;
    }

    //std::cout << "Non-Linear Warping Modes Derived" << endl << endl;
    return std::tie(coef_z, stw);

}
