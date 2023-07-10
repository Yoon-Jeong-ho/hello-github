#pragma once

#include "HigherOrderBeamLib.h"

VectorXi BeamNodeDof(VectorXi Nmode, VectorXi Nels, VectorXi sbdof, int bi, int ni)
{
    VectorXi fix1;

    if (bi == 0)
    {
        if (ni == 0)
            fix1.setLinSpaced(Nmode(bi), 0, Nmode(bi) - 1);
        else if (ni == Nels(bi))
            fix1.setLinSpaced(Nmode(bi), sbdof(0) - 2 * Nmode(bi), sbdof(0) - Nmode(bi) - 1);
        else
            fix1.setLinSpaced(Nmode(bi), 2 * ni * Nmode(bi), (2 * ni + 1) * Nmode(bi) - 1);
    }
    else
    {
        if (ni == 0)
            fix1.setLinSpaced(Nmode(bi), sbdof(bi - 1), sbdof(bi - 1) + Nmode(bi) - 1);
        else if (ni == Nels(bi))
            fix1.setLinSpaced(Nmode(bi), sbdof(bi) - 2 * Nmode(bi), sbdof(bi) - Nmode(bi) - 1);
        else
            fix1.setLinSpaced(Nmode(bi), sbdof(bi - 1) + 2 * ni * Nmode(bi), sbdof(bi - 1) + (2 * ni + 1) * Nmode(bi) - 1);
    }

    return fix1;
}