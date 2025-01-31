// grid.h
#ifndef GRID1D_H
#define GRID1D_H

#include <vector>
#include <iostream>
#include "helpers.h"


class Grid1D {
public:
    // Costruttore
    Grid1D(int nx, double dx, std::vector<std::vector<double>> bound);

    // Inizializza i valori iniziali
    void initialize(const std::string &cas);

    const std::vector<std::vector<double>> &getBound() const;

    // Getter per le dimensioni
    int getNx() const { return nx; }

    double getDx() const;

    // Accesso ai dati
    double& rho(int i);
    double& u(int i);
    double& p(int i);
    double E(int i);
    double H(int i);
    std::vector<double> asCoords(int i);
    void updateRho(std::vector<double> & rhonew);
    void updateU(std::vector<double> & unew);
    void updateEP(std::vector<double> & Enew);
    void update(std::vector<std::vector<double>> & unew);
    double& ghostCell(int i, int var);
    void applyBound();
    std::vector<double> checkSmooth();
    std::vector<double> integrate();
    double sMax();
private:
    int nx; // Numero di celle
    double dx; // Spaziatura del reticolo

    // Variabili fisiche
    std::vector<double> rho_; // Densità
    std::vector<double> u_;   // Velocità x
    std::vector<double> p_;   // Pressione
    std::vector<double> E_;
    std::vector<std::vector<double>> bound;
    std::vector<double> intConserved;
    double gamma=1.4;
    std::string casTest;
public:
    const std::vector<double> &getIntConserved() const;
};


#endif // GRID_H
