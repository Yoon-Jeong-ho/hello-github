#pragma once

#include "HigherOrderBeamLib.h"

class Hermite
{
private:
    double N1(double k) { return (pow(1 - k, 2) * (2 + k)) / 4; }
    double N2(double L, double k) { return (L / 8) * pow(1 - k, 2) * (1 + k); }
    double N3(double k) { return (pow(1 + k, 2) * (2 - k)) / 4; }
    double N4(double L, double k) { return (L / 8) * pow(1 + k, 2) * (k - 1); }
    double dN1(double L, double k) { return (pow(k - 1, 2) / 4 + ((2 * k - 2) * (k + 2)) / 4) * 2 / L; }
    double dN2(double L, double k) { return ((L * pow(k - 1, 2)) / 8 + (L * (2 * k - 2) * (k + 1)) / 8) * 2 / L; }
    double dN3(double L, double k) { return (-pow(k + 1, 2) / 4 - ((2 * k + 2) * (k - 2)) / 4) * 2 / L; }
    double dN4(double L, double k) { return ((L * pow(k + 1, 2)) / 8 + (L * (2 * k + 2) * (k - 1)) / 8) * 2 / L; }
    double ddN1(double L, double k) { return (3 * k) / 2 * pow(2 / L, 2); }
    double ddN2(double L, double k) { return ((L * (k + 1)) / 4 + (L * (2 * k - 2)) / 4) * pow(2 / L, 2); }
    double ddN3(double L, double k) { return -(3 * k) / 2 * pow(2 / L, 2); }
    double ddN4(double L, double k) { return ((L * (k - 1)) / 4 + (L * (2 * k + 2)) / 4) * pow(2 / L, 2); }

public:
    MatrixXd N(double L, double k, int Nmode)
    {
        MatrixXd eye(Nmode, Nmode), Ni(Nmode, 4 * Nmode);
        eye.setIdentity(), Ni.setZero();
        Ni.block(0, 0, Nmode, Nmode) = eye * N1(k);
        Ni.block(0, Nmode, Nmode, Nmode) = eye * N2(L, k);
        Ni.block(0, 2 * Nmode, Nmode, Nmode) = eye * N3(k);
        Ni.block(0, 3 * Nmode, Nmode, Nmode) = eye * N4(L, k);

        return Ni;
    }
    MatrixXd dN(double L, double k, int Nmode)
    {
        MatrixXd eye(Nmode, Nmode), Ni(Nmode, 4 * Nmode);
        eye.setIdentity(), Ni.setZero();
        Ni.block(0, 0, Nmode, Nmode) = eye * dN1(L, k);
        Ni.block(0, Nmode, Nmode, Nmode) = eye * dN2(L, k);
        Ni.block(0, 2 * Nmode, Nmode, Nmode) = eye * dN3(L, k);
        Ni.block(0, 3 * Nmode, Nmode, Nmode) = eye * dN4(L, k);

        return Ni;
    }
    MatrixXd ddN(double L, double k, int Nmode)
    {
        MatrixXd eye(Nmode, Nmode), Ni(Nmode, 4 * Nmode);
        eye.setIdentity(), Ni.setZero();
        Ni.block(0, 0, Nmode, Nmode) = eye * ddN1(L, k);
        Ni.block(0, Nmode, Nmode, Nmode) = eye * ddN2(L, k);
        Ni.block(0, 2 * Nmode, Nmode, Nmode) = eye * ddN3(L, k);
        Ni.block(0, 3 * Nmode, Nmode, Nmode) = eye * ddN4(L, k);

        return Ni;
    }
};