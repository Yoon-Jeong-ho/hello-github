#pragma once

#include "HigherOrderBeamLib.h"
#include "Calculus.h"

class Edge
{
private:
    Array<double50, -1, 1> var1;
    ArrayXd var2;
    Calculus operate;
public:
	Array<double50, -1, -1> Anglemp(Matrix<double50, -1, -1> csc, MatrixXi cscc)
	{
        var1.setZero(cscc.rows(), 1);
        ArrayXi e;
        for (int i = 0; i < cscc.rows(); i++)
        {
            e = cscc(i, all);
            var1(i) = atan2(csc(e(1), 1) - csc(e(0), 1), csc(e(1), 0) - csc(e(0), 0));
        }
        return var1;
	}
    ArrayXd Angle(MatrixXd csc, MatrixXi cscc)
    {
        var2.setZero(cscc.rows(), 1);
        ArrayXi e;
        for (int i = 0; i < cscc.rows(); i++)
        {
            e = cscc(i, all);
            var2(i) = atan2(csc(e(1), 1) - csc(e(0), 1), csc(e(1), 0) - csc(e(0), 0));
        }
        return var2;
    }
    Array<double50, -1, -1> Lengthmp(Matrix<double50, -1, -1> csc, MatrixXi cscc)
    {
        var1.setZero(cscc.rows(), 1);
        ArrayXi e;
        for (int i = 0; i < cscc.rows(); i++)
        {
            e = cscc(i, all);
            var1(i) = sqrt(pow((csc(e(1), 1) - csc(e(0), 1)), 2) + pow((csc(e(1), 0) - csc(e(0), 0)), 2));
        }
        return var1;
    }
    ArrayXd Length(MatrixXd csc, MatrixXi cscc)
    {
        var2.setZero(cscc.rows(), 1);
        ArrayXi e;
        for (int i = 0; i < cscc.rows(); i++)
        {
            e = cscc(i, all);
            var2(i) = sqrt(pow((csc(e(1), 1) - csc(e(0), 1)), 2) + pow((csc(e(1), 0) - csc(e(0), 0)), 2));
        }
        return var2;
    }
    Matrix<double50, -1, -1> Continuity(vector<Matrix<double50, -1, -1>> coef, int e, double50 s, int d)
    {
        Index NC = coef[0].rows(), NM = coef[0].cols();
        int NE = coef.size();
        vector<Matrix<double50, -1, -1>> dcoef = operate.diffcoefmp(coef, d);
        ArrayXi spow;
        Array<double50, -1, 1> ss;
        Matrix<double50, -1, -1> sp, fval;
        spow.setLinSpaced(NC, 0, NC - 1), ss.setOnes(NC, 1), sp.setZero(1, NC);
        ss *= s;
        for (int i = 0;i < NC;i++)
        {
            sp(i) = pow(ss(i), spow(i));
        }
        fval = sp * dcoef[e].block(0, 0, NC, NM);

        return fval;
    }
    MatrixXd Continuity2(vector<MatrixXd> coef, int e, ArrayXd s, int d)
    {
        int NC = coef[0].rows(), NM = coef[0].cols(), NE = coef.size(), NS = s.rows();
        vector<MatrixXd> dcoef = operate.diffcoef(coef, d);
        ArrayXi spow;
        MatrixXd fval, si(NS, NC);
        spow.setLinSpaced(NC, 0, NC - 1), si.setZero();
        for (int i = 0;i < NC;i++)
        {
            for (int j = 0; j < NS; j++)
                si(j, i) = pow(s(j), spow(i));
        }
        fval = si * dcoef[e].block(0, 0, NC, NM);
        return fval;
    }
};

