#pragma once

#include "HigherOrderBeamLib.h"

class Center
{
private:
    class Rotation
    {
    public:
        Array<double50, -1, 1> Coordinates(Matrix<double50, -1, -1> csc, MatrixXi cscc, Array<double50, -1, 1> alpha, Array<double50, -1, 1> Le)
        {
            Array<double50, -1, 1> X = csc(cscc(all, 0), 0), Y = csc(cscc(all, 0), 1), cR(2, 1);
            cR(0) = (X.cwiseProduct(Le).sum() + Le.square().cwiseProduct(0.5 * cos(alpha)).sum()) / Le.sum();
            cR(1) = (Y.cwiseProduct(Le).sum() + Le.square().cwiseProduct(0.5 * sin(alpha)).sum()) / Le.sum();

            return cR;
        }
        double50 Angle(Matrix<double50, -1, -1> csc, MatrixXi cscc, Array<double50, -1, 1> alpha, Array<double50, -1, 1> Le, Array<double50, -1, 1> cR)
        {
            Array<double50, -1, 1> X = csc(cscc(all, 0), 0), Y = csc(cscc(all, 0), 1);
            Array<double50, -1, 1> xc = X - cR(0), yc = Y - cR(1);
            double50 a, b, c;
            a = yc.square().cwiseProduct(Le).sum() + Le.square().cwiseProduct(sin(alpha)).cwiseProduct(yc).sum() + Le.cube().cwiseProduct(sin(alpha).square()).sum() * 1 / 3;
            b = xc.square().cwiseProduct(Le).sum() + Le.square().cwiseProduct(cos(alpha).cwiseProduct(xc)).sum() + Le.cube().cwiseProduct(cos(alpha).square()).sum() * 1 / 3;
            c = xc.cwiseProduct(yc.cwiseProduct(Le)).sum() + Le.square().cwiseProduct(xc.cwiseProduct(sin(alpha))).sum() * 1 / 2 + Le.square().cwiseProduct(yc.cwiseProduct(cos(alpha))).sum() * 1 / 2 + Le.cube().cwiseProduct(sin(alpha).cwiseProduct(cos(alpha))).sum() * 1 / 3;

            return 0.5 * atan(2 * c / (b - a));
        }
        ArrayXd gcenter(MatrixXd csc, MatrixXi cscc, ArrayXd alpha, ArrayXd Le)
        {
            ArrayXd X = csc(cscc(all, 0), 0), Y = csc(cscc(all, 0), 1), cR(2, 1);
            cR(0) = (X.cwiseProduct(Le).sum() + Le.square().cwiseProduct(0.5 * cos(alpha)).sum()) / Le.sum();
            cR(1) = (Y.cwiseProduct(Le).sum() + Le.square().cwiseProduct(0.5 * sin(alpha)).sum()) / Le.sum();

            return cR;
        }
    };
    class Torsion
    {
    public:
        double50 Angle(Matrix<double50, -1, -1> csc, MatrixXi cscc, Array<double50, -1, 1> alpha, Array<double50, -1, 1> Le, double50 betaR)
        {
            double50 beta = 0.5 * atan2(Le.cwiseProduct(sin(2 * alpha)).sum(), Le.cwiseProduct(cos(2 * alpha)).sum());
            if (abs(cos(alpha - 0).cwiseProduct(sin(alpha - 0)).cwiseProduct(Le).sum()) < Le.sum() * 1e-14)
                beta = 0;
            else if (abs(cos(alpha - betaR).cwiseProduct(sin(alpha - betaR)).cwiseProduct(Le).sum()) < Le.sum() * 1e-14)
                beta = betaR;

            return beta;
        }
        Array<double50, -1, 1> Coordinates(Matrix<double50, -1, -1> csc, MatrixXi cscc, Array<double50, -1, 1> alpha, Array<double50, -1, 1> Le, double50 beta)
        {
            Array<double50, -1, 1> Xe = csc(cscc(all, 0), 0), Ye = csc(cscc(all, 0), 1);
            Matrix<double50, 2, 2> A1, a;
            Matrix<double50, -1, -1> A2(2, 1), b(2, 1);
            A1.setZero(), A2.setZero();
            for (int i = 0; i < cscc.rows(); i++)
            {
                a(0, 0) = Le(i) * sin(alpha(i)) * cos(alpha(i) - beta);
                a(0, 1) = -Le(i) * cos(alpha(i)) * cos(alpha(i) - beta);
                a(1, 0) = Le(i) * sin(alpha(i)) * sin(alpha(i) - beta);
                a(1, 1) = -Le(i) * cos(alpha(i)) * sin(alpha(i) - beta);

                b(0) = Le(i) * Xe(i) * sin(alpha(i)) * cos(alpha(i) - beta) - (Le(i) * Ye(i) * cos(alpha(i)) * cos(alpha(i) - beta));
                b(1) = Le(i) * Xe(i) * sin(alpha(i)) * sin(alpha(i) - beta) - (Le(i) * Ye(i) * cos(alpha(i)) * sin(alpha(i) - beta));

                A1 += a, A2 += b;
            }

            return A1.inverse() * A2;
        }
    };
public:
    Rotation Rot;
    Torsion Tor;
};

