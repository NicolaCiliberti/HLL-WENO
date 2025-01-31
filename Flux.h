//
// Created by UTENTE on 18/12/2024.
//

#ifndef PROJET_FLUX_H
#define PROJET_FLUX_H

#include "Grid.h"
#include "Grid1D.h"


class Flux {
private:
    double dx,dy,dt,sxmax,smax,symax;
    int Nx,Ny;
public:
    Flux(double dx, int Nx,double dt): dx(dx),Nx(Nx),dt(dt){
        dy=0;
        Ny=0;

    }

    void setSymax(double symax);

    void setSxmax(double sxmax);

    void setSmax(double smax);;
    Flux(double dx,double dy, int Nx,int Ny,double dt): dx(dx),Nx(Nx),dy(dy),Ny(Ny),dt(dt){
    };
    std::vector<double> evaluateFFlux(Grid & grid, int i, int j);
    std::vector<double> evaluateFFlux(std::vector<double> & ul, std::vector<double> & ur,double gamma);
    std::vector<double> evaluateGFlux(std::vector<double> & ul, std::vector<double> & ur,double gamma);
    std::vector<double> evaluateGFlux(Grid & grid, int i, int j);
    std::vector<double> evaluateF(double rho, double u, double p, double E);
    std::vector<double> evaluateF(double rho, double u, double v, double p, double E);
    std::vector<double> evaluateG(double rho, double u, double v, double p, double E);
    std::vector<double> evaluateFlux(Grid1D &grid, int i);
    std::vector<double> evaluateF(std::vector<double> & coords, double gamma);
    std::vector<double> evaluateG(std::vector<double> & coords, double gamma);
    std::pair<std::vector<double>,std::vector<double>> WENO1D(Grid1D & grid, int j);
    std::vector<double> WENO5(std::vector<std::vector<double>> & fi);
    std::vector<std::vector<double>>  WENO(Grid &grid, int i, int j, int dir);
    std::vector<double> evaluateLF(std::vector<double> & ul, std::vector<double> & ur,double gamma);
    std::vector<double> evaluateFlux(std::vector<double> & ul, std::vector<double> & ur);
    std::vector<double> evaluateLFx(std::vector<double> & ul, std::vector<double> & ur, double gamma);
    std::vector<double> evaluateLFy(std::vector<double> & ul, std::vector<double> & ur, double gamma);
    std::vector<double> evaluateHLLCF(std::vector<double> & ul, std::vector<double> & ur, double gamma);
    std::vector<double> evaluateHLLCG(std::vector<double> & ub, std::vector<double> & ut, double gamma);
    std::vector<double> evaluateFContact(std::vector<double> & ul, std::vector<double> & ur, double sl, double sr, double sc, double gamma);
    std::vector<double> evaluateGContact(std::vector<double> & ub, std::vector<double> & ut, double sb, double st, double sc, double gamma);
    void updateSmax(double s);
};


#endif //PROJET_FLUX_H
