//
// Created by UTENTE on 20/12/2024.
//

#include "Print.h"
#include <fstream>
#include <iomanip>

void Print::writeVariableToFile(Grid & grid, const std::string &filename, const std::string &variable) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << " per scrivere i dati." << std::endl;
        return;
    }

    int nx = grid.getNx();
    int ny = grid.getNy();

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double value = 0.0;
            if (variable == "rho") {
                value = grid.rho(i, j);
            } else if (variable == "u") {
                value = grid.u(i, j);
            } else if (variable == "v") {
                value = grid.v(i, j);
            } else if (variable == "p") {
                value = grid.p(i, j);
            } else if (variable == "E") {
                value = grid.E(i, j);
            } else {
                std::cerr << "Errore: variabile " << variable << " non riconosciuta." << std::endl;
                return;
            }
            file << "  " << value;
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Dati di " << variable << " scritti con successo su " << filename << std::endl;
}

void Print::printGridToFiles(Grid & grid) {
    writeVariableToFile(grid,"rho.txt", "rho");
    writeVariableToFile(grid,"u.txt", "u");
    writeVariableToFile(grid,"v.txt", "v");
    writeVariableToFile(grid,"p.txt", "p");
    writeVariableToFile(grid,"E.txt", "E");
}

void Print::printGridToFiles1D(Grid1D & grid) {
    writeVariableToFile1D(grid,"rho.txt", "rho");
    writeVariableToFile1D(grid,"u.txt", "u");
    writeVariableToFile1D(grid,"p.txt", "p");
    writeVariableToFile1D(grid,"E.txt", "E");
}

void Print::exportToVTK(Grid & grid, const std::string& filename) {
    std::ofstream file(filename);
    double gamma = grid.getGamma();

    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << " per scrivere i dati." << std::endl;
        return;
    }

    int nx = grid.getNx();
    int ny = grid.getNy();
    double dx = grid.getDx();
    double dy = grid.getDy();

    // Header del file VTK
    file << "# vtk DataFile Version 3.0\n";
    file << "Grid data\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING " << dx << " " << dy << " 1\n";
    file << "POINT_DATA " << nx * ny << "\n";

    // Scrittura di rho
    file << "SCALARS rho double\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << grid.rho(i, j) << "\n";
        }
    }

    // Scrittura di u
    file << "SCALARS u double\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << grid.u(i, j) << "\n";
        }
    }

    // Scrittura di v
    file << "SCALARS v double\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << grid.v(i, j) << "\n";
        }
    }

    // Scrittura di p
    file << "SCALARS p double\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << grid.p(i, j) << "\n";
        }
    }

    // Scrittura di E
    file << "SCALARS E double\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << grid.E(i, j) << "\n";
        }
    }

    file.close();
    std::cout << "File VTK scritto con successo: " << filename << std::endl;
}

void Print::writeVariableToFile1D(Grid1D & grid1, const std::string &filename, const std::string &variable) {
    double gamma = 1.4;
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << " per scrivere i dati." << std::endl;
        return;
    }

    int nx = grid1.getNx();

    for (int i = 0; i < nx; ++i) {
        double value = 0;
        if (variable == "rho") {
            value = grid1.rho(i);
        } else if (variable == "u") {
            value = grid1.u(i);
        } else if (variable == "p") {
            value = grid1.p(i);
        } else if (variable == "E") {
            value = grid1.p(i)/((gamma-1)*grid1.rho(i));
        } else {
            std::cerr << "Errore: variabile " << variable << " non riconosciuta." << std::endl;
            return;
        }
        file << "  " << value;
    }
    file << std::endl;
    file.close();
    std::cout << "Dati di " << variable << " scritti con successo su " << filename << std::endl;
}

void Print::printSmoothGrid(int n) {
    std::ofstream file("smooth/" + std::to_string(n)+"x"+std::to_string(n));
    double dx = 1.0/n;
    double value;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            value = mediaIntegrale(i*dx,(i+1)*dx,j*dx,(j+1)*dx,0.5);
            file << "  " << value;
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Dati smooth scritti con successo" << std::endl;
}

