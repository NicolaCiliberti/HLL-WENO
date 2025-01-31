//
// Created by UTENTE on 19/12/2024.
//

#ifndef PROJET_TIME_H
#define PROJET_TIME_H


#include <vector>
#include "Flux.h"

class Time {
private:
    double dt;
    std::vector<double> conserved;
public:
    Time(double dt, std::vector<double> cons): dt(dt), conserved(cons){};
    Grid& advance(Flux & F, Grid & G);
    Grid& advanceWENO(Grid &grid, Flux &flux);
    Grid1D& advance1D(Flux &flux, Grid1D &grid);
    void isConserved(std::vector<double> newConserved, const std::vector<double> & residue0);
    Grid1D& advance1DWENO(Grid1D &grid, Flux &flux);
    Grid1D& evolveRK3(Grid1D &grid, Flux &flux);
    Grid1D& advance1DLF(Flux &flux, Grid1D &grid);
    Grid1D& evolveRK3WENO(Grid1D &grid, Flux &flux);
    Grid& evolveRK3WENO(Grid &grid, Flux &flux);
    Grid& evolveRK3WENOHLLC(Grid &grid, Flux &flux);
    Grid1D& advance1DHLLC(Flux &flux, Grid1D &grid);
    Grid& advanceHLLC(Flux &F, Grid &grid);
    Grid& evolveRK3(Grid &grid, Flux &flux);
    Grid& evolveRK3WENOHLLCV2(Grid &grid, Flux &flux);
    Grid1D& evolveRK3Transp(Grid1D &grid, Flux &flux);
    Grid1D& testRK3(Grid1D &grid, Flux &flux);
};


#endif //PROJET_TIME_H
