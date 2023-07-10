#pragma once

#include "HigherOrderBeamLib.h"

VectorXi noddof3(VectorXi Nel, VectorXi Nmode, int b, int n)
{
    int pdof = 0, Neli, Ndofi;
    for (int i = 0; i < b - 1; i++)
    {
        Neli = Nel(i);
        Ndofi = Nmode(i) * 2;
        pdof += (Neli + 1) * Ndofi;
    }

    int Ndofc = Nmode(b - 1) * 2;
    VectorXi dof;
    if (n != INT_MAX)
    {
        dof.setLinSpaced(n * Ndofc - (n - 1) * Ndofc, (n - 1) * Ndofc + 1, n * Ndofc);
        dof += pdof * VectorXi::Ones(dof.size());
    }
    else
    {
        dof.setLinSpaced((Nel(b - 1) + 1) * Ndofc, 0, (Nel(b - 1) + 1) * Ndofc - 1);
        dof += pdof * VectorXi::Ones(dof.size());
    }

    return dof;
}