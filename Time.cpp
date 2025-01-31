//
// Created by UTENTE on 19/12/2024.
//

#include "Time.h"

//Avancement en temps 2D avec Euler explicite HLL
Grid& Time::advance(Flux &F, Grid &grid) {
    double gamma = grid.getGamma();
    int Nx = grid.getNx();
    int Ny = grid.getNy();
    double dx = grid.getDx();
    double dy = grid.getDy();
    std::vector<double> tempCoords(4);
    std::vector<std::vector<double>> unp1(4, std::vector<double>(Nx*Ny, 0));
    std::vector<double> source(4,0);
    std::vector<double> ut,ub,u,ul,ur;
    F.setSxmax(grid.getSxmax());
    F.setSymax(grid.getSymax());
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if(grid.getProb()=="RIT"){
                source[2] = -grid.rho(i,j)*9.81;
            }
            u = grid.asCoords(i,j);
            ut = grid.asCoords(i,j+1);
            ub = grid.asCoords(i,j-1);
            ul = grid.asCoords(i-1,j);
            ur = grid.asCoords(i+1,j);
            tempCoords = grid.asCoords(i,j) - dt/dx * (F.evaluateFFlux(u,ur,gamma) - F.evaluateFFlux(ul,u,gamma))
                    - dt/dy * (F.evaluateGFlux(u,ut,gamma) - F.evaluateGFlux(ub,u,gamma)); //+ dt * source; //Manca sorgente
            for (int k=0; k<4; k++){
                unp1[k][i*Ny + j] = tempCoords[k];
            }
        }
    }
    grid.applyBound();
    grid.update(unp1);
    return grid;
}

//Avancement en temps 2D avec WENO5 + Euler explicite
Grid& Time::advanceWENO(Grid &grid, Flux &flux) {
    double gamma = grid.getGamma();
    int Nx = grid.getNx();
    int Ny = grid.getNy();
    double dx = grid.getDx();
    double dy = grid.getDy();
    std::vector<double> tempCoords(4);
    std::vector<std::vector<double>> unp1(4, std::vector<double>(Nx*Ny, 0));
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            auto fluxes = flux.WENO(grid,i,j,2);
            auto fluxes2 = flux.WENO(grid,i-1,j,0);
            auto fluxes3 = flux.WENO(grid,i,j-1,1);
            std::vector<double> source(4,0); //{0,0,grid.rho(i,j)*3,grid.rho(i,j)*grid.v(i,j)*3};
            tempCoords = grid.asCoords(i,j) - dt / dx *
                                              (flux.evaluateFFlux(fluxes[0],fluxes[1],gamma)
                                               - flux.evaluateFFlux(fluxes2[0],fluxes2[1],gamma))
                                              - dt / dy *
                                              (flux.evaluateGFlux(fluxes[2],fluxes[3],gamma)
                                               - flux.evaluateGFlux(fluxes3[0],fluxes3[1],gamma)); //Manca sorgente
            for (int k=0; k<4; k++){
                unp1[k][i*Ny + j] = tempCoords[k];
            }
        }
    }
    grid.applyBound();
    grid.update(unp1);
    return grid;
}

//Avancement en temps 2D avec Euler explicite HLLC
Grid& Time::advanceHLLC(Flux &F, Grid &grid) {
    double gamma = grid.getGamma();
    int Nx = grid.getNx();
    int Ny = grid.getNy();
    double dx = grid.getDx();
    double dy = grid.getDy();
    std::vector<double> tempCoords(4);
    std::vector<std::vector<double>> unp1(4, std::vector<double>(Nx*Ny, 0));
    std::vector<double> source(4,0);
    std::vector<double> ut,ub,u,ul,ur;
    F.setSxmax(grid.getSxmax());
    F.setSymax(grid.getSymax());
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            if(grid.getProb()=="RIT"){
                source[2] = -grid.rho(i,j)*9.81;
            }
            u = grid.asCoords(i,j);
            ut = grid.asCoords(i,j+1);
            ub = grid.asCoords(i,j-1);
            ul = grid.asCoords(i-1,j);
            ur = grid.asCoords(i+1,j);
            tempCoords = grid.asCoords(i,j) - dt/dx * (F.evaluateHLLCF(u,ur,gamma) - F.evaluateHLLCF(ul,u,gamma))
                         - dt/dy * (F.evaluateHLLCG(u,ut,gamma) - F.evaluateHLLCG(ub,u,gamma)); //+ dt * source; //Manca sorgente
            for (int k=0; k<4; k++){
                unp1[k][i*Ny + j] = tempCoords[k];
            }
        }
    }
    grid.applyBound();
    grid.update(unp1);
    return grid;
}

//Controle des quantitées conservatifs
void Time::isConserved(std::vector<double> newConserved, const std::vector<double>& conserved0) {
    int dim = newConserved.size();
    std::vector<double> residue(dim);
    std::vector<double> alt_residue(dim);
    std::cout << "Previous conserved" << std::endl;
    print(conserved);
    std::cout << "New conserved" << std::endl;
    print(newConserved);
    conserved = newConserved;
    for (int i = 0; i<dim; i++){
        residue[i] = (conserved[i]/newConserved[i]) - 1;
        alt_residue[i] = (conserved0[i]/newConserved[i]) - 1;
    }
    std::cout << "Residue" << std::endl;
    //print(residue);
    //std::cout << "Residue" << std::endl;
    print(alt_residue);
}

//Avancement en temps 1D avec Euler Explicite HLL
Grid1D& Time::advance1D(Flux &flux, Grid1D &grid) {
    int Nx = grid.getNx();
    double dx = grid.getDx();
    std::vector<std::vector<double>> unp1(3, std::vector<double>(Nx, 0));
    std::vector<double> tempCoords(3);
    std::vector<double> fluxLeft,fluxRight;
    for (int i = 0; i < Nx; i++) {
        fluxLeft = flux.evaluateFlux(grid, i - 1);
        fluxRight = flux.evaluateFlux(grid, i);
        tempCoords = grid.asCoords(i) - dt / dx * (fluxRight - fluxLeft);
        for (int k = 0; k < 3; ++k) {
            unp1[k][i] = tempCoords[k];
        }
    }
    grid.applyBound();// Aggiorna la griglia con i nuovi valori
    grid.update(unp1);
    return grid;
}

//Avancement en temps 1D avec Euler Explicite Rusanov
Grid1D& Time::advance1DLF(Flux &flux, Grid1D &grid) {
    double gamma = 1.4;
    int Nx = grid.getNx();
    double dx = grid.getDx();
    std::vector<std::vector<double>> unp1(3, std::vector<double>(Nx, 0));
    std::vector<double> tempCoords(3);
    std::vector<double> fluxLeft,fluxRight,ul,ur;
    flux.setSmax(grid.sMax());
    for (int i = 0; i < Nx; i++) {
        ul = grid.asCoords(i-1);
        ur = grid.asCoords(i);
        fluxLeft = flux.evaluateLF(ul,ur,gamma);
        ul = grid.asCoords(i);
        ur = grid.asCoords(i+1);
        fluxRight = flux.evaluateLF(ul,ur,gamma);
        tempCoords = grid.asCoords(i) - dt / dx * (fluxRight - fluxLeft);
        for (int k = 0; k < 3; ++k) {
            unp1[k][i] = tempCoords[k];
        }
    }
    grid.applyBound();// Aggiorna la griglia con i nuovi valori
    grid.update(unp1);
    flux.updateSmax(grid.sMax());
    return grid;
}

//Avancement en temps 1D avec Euler Explicite HLLC
Grid1D& Time::advance1DHLLC(Flux &flux, Grid1D &grid) {
    double gamma = 1.4;
    int Nx = grid.getNx();
    double dx = grid.getDx();
    std::vector<std::vector<double>> unp1(3, std::vector<double>(Nx, 0));
    std::vector<double> tempCoords(3);
    std::vector<double> fluxLeft,fluxRight,ul,ur;
    flux.setSmax(grid.sMax());
    for (int i = 0; i < Nx; i++) {
        ul = grid.asCoords(i-1);
        ur = grid.asCoords(i);
        fluxLeft = flux.evaluateHLLCF(ul,ur,gamma);
        ul = grid.asCoords(i);
        ur = grid.asCoords(i+1);
        fluxRight = flux.evaluateHLLCF(ul,ur,gamma);
        tempCoords = grid.asCoords(i) - dt / dx * (fluxRight - fluxLeft);
        for (int k = 0; k < 3; ++k) {
            unp1[k][i] = tempCoords[k];
        }
    }
    grid.applyBound();// Aggiorna la griglia con i nuovi valori
    grid.update(unp1);
    flux.updateSmax(grid.sMax());
    return grid;
}

//Avancement en temps 1D avec WENO5 + Euler Explicite
Grid1D& Time::advance1DWENO(Grid1D &grid, Flux &flux) {
    int Nx = grid.getNx();
    double dx = grid.getDx();
    std::vector<std::vector<double>> unp1(3, std::vector<double>(Nx, 0));
    std::vector<double> tempCoords(3);
    std::vector<double> fluxLeft,fluxRight;
    std::pair<std::vector<double>,std::vector<double>> fPM;
    fluxLeft = flux.evaluateFlux(grid, -1);
    fluxRight = flux.evaluateFlux(grid, 0);
    tempCoords = grid.asCoords(0) - dt / dx * (fluxRight - fluxLeft);
    for (int k = 0; k < 3; ++k) {
        unp1[k][0] = tempCoords[k];
    }
    for (int i = 1; i < Nx-1; i++) {
        fPM = flux.WENO1D(grid,i-1);
        fluxLeft = flux.evaluateFlux(fPM.first,fPM.second);
        fPM = flux.WENO1D(grid,i);
        fluxRight = flux.evaluateFlux(fPM.first,fPM.second);
        tempCoords = grid.asCoords(i) - dt / dx * (fluxRight - fluxLeft);
        //print(tempCoords);
        for (int k = 0; k < 3; ++k) {
            unp1[k][i] = tempCoords[k];
        }
    }
    fluxLeft = flux.evaluateFlux(grid, Nx-1);
    fluxRight = flux.evaluateFlux(grid, Nx);
    tempCoords = grid.asCoords(Nx-1) - dt / dx * (fluxRight - fluxLeft);
    for (int k = 0; k < 3; ++k) {
        unp1[k][Nx-1] = tempCoords[k];
    }
    grid.applyBound();// Aggiorna la griglia con i nuovi valori
    grid.update(unp1);
    return grid;
}

//Avancement en temps 2D WENO5 + SSPRK3
Grid1D& Time::evolveRK3WENO(Grid1D &grid, Flux &flux) {
    double gamma = 1.4;
    int dx = grid.getDx();
    Grid1D gridTemp = grid;
    flux.updateSmax(grid.sMax());
    //Grid1D grid2(grid.getDx(), grid.getNx(),grid.getBound());
    std::vector<std::vector<double>> u1 (3,std::vector<double> (grid.getNx()));
    std::vector<double> tempCoords(3);
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(grid, j);
        auto fluxes2 = flux.WENO1D(grid, j-1);
        tempCoords = grid.asCoords(j) - dt / grid.getDx() *
                (flux.evaluateFlux(fluxes.first,fluxes.second)
                - flux.evaluateFlux(fluxes2.first,fluxes2.second));
        //print(tempCoords);
        for (int k = 0; k < 3; ++k) {
            u1[k][j] = tempCoords[k];
        }
    }
    //printf("Updating GridTemp 1 \n");
    gridTemp.update(u1);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> u2 (3,std::vector<double> (grid.getNx()));
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(gridTemp, j);
        auto fluxes2 = flux.WENO1D(gridTemp, j-1);
        tempCoords = 0.75 * grid.asCoords(j) + 0.25 * (gridTemp.asCoords(j) - dt/grid.getDx()
                * (flux.evaluateFlux(fluxes.first,fluxes.second)
                   - flux.evaluateFlux(fluxes2.first,fluxes2.second)));
        //print(tempCoords);
        for (int k = 0; k < 3; ++k) {
            u2[k][j] = tempCoords[k];
        }
    }

    //printf("Updating GridTemp 2 \n");
    gridTemp.update(u2);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> unp1 (3,std::vector<double> (grid.getNx()));
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(gridTemp, j);
        auto fluxes2 = flux.WENO1D(gridTemp, j-1);
        tempCoords = 1.0/3 * grid.asCoords(j) + 2.0/3 * (gridTemp.asCoords(j) - dt/grid.getDx()
                * (flux.evaluateFlux(fluxes.first,fluxes.second)
                   - flux.evaluateFlux(fluxes2.first,fluxes2.second)));
        for (int k = 0; k < 3; ++k) {
            unp1[k][j] = tempCoords[k];
        }
    }
    //printf("Updating Grid \n");
    grid.update(unp1);
    flux.updateSmax(grid.sMax());
    grid.applyBound();
    return grid;
}

//Avancement en temps 1D SSPRK3 HLL
Grid1D& Time::evolveRK3(Grid1D &grid, Flux &flux) {
    Grid1D gridTemp = grid;
    //Grid1D grid2(grid.getDx(), grid.getNx(),grid.getBound());
    std::vector<std::vector<double>> u1 (3,std::vector<double> (grid.getNx()));
    std::vector<double> tempCoords(3);
    for (int j = 0; j < grid.getNx(); j++) {
        tempCoords = grid.asCoords(j)
                - dt / grid.getDx() * (flux.evaluateFlux(grid,j) - flux.evaluateFlux(grid,j-1));
        //print(tempCoords);
        for (int k = 0; k < 3; ++k) {
            u1[k][j] = tempCoords[k];
        }
    }
    //printf("Updating GridTemp 1 \n");
    gridTemp.update(u1);
    gridTemp.applyBound();

    std::vector<std::vector<double>> u2 (3,std::vector<double> (grid.getNx()));
    for (int j = 0; j < grid.getNx(); j++) {
        tempCoords = 0.75 * grid.asCoords(j) + 0.25 * (gridTemp.asCoords(j)
                - dt / grid.getDx() * (flux.evaluateFlux(gridTemp,j) - flux.evaluateFlux(gridTemp,j-1)));
        //print(tempCoords);
        for (int k = 0; k < 3; ++k) {
            u2[k][j] = tempCoords[k];
        }
    }

    //printf("Updating GridTemp 2 \n");
    gridTemp.update(u2);
    gridTemp.applyBound();

    std::vector<std::vector<double>> unp1 (3,std::vector<double> (grid.getNx()));
    for (int j = 0; j < grid.getNx(); j++) {
        tempCoords = 1.0/3 * grid.asCoords(j) + 2.0/3 * (gridTemp.asCoords(j)
                - dt / grid.getDx() * (flux.evaluateFlux(gridTemp,j) - flux.evaluateFlux(gridTemp,j-1)));
        for (int k = 0; k < 3; ++k) {
            unp1[k][j] = tempCoords[k];
        }
    }
    //printf("Updating Grid \n");
    grid.update(unp1);
    grid.applyBound();
    return grid;
}

//Avancement en temps 2D SSPRK3 HLL
Grid& Time::evolveRK3(Grid &grid, Flux &flux) {
    double gamma = grid.getGamma();
    double g = 1;
    Grid gridTemp = grid;
    flux.updateSmax(grid.sMax());
    std::vector<std::vector<double>> u1 (4,std::vector<double> (grid.getNx()*grid.getNy())),
            fluxes,fluxes2,fluxes3;
    std::vector<double> tempCoords(4), ur, ul, ub, ut, u;
    std::vector<double> source(4,0);
    for(int i = 0; i < grid.getNx(); i++) {
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -grid.rho(i,j)*g;
                source[3] = -grid.v(i,j)*grid.rho(i,j)*g;
            }
            u = grid.asCoords(i,j);
            ut = grid.asCoords(i,j+1);
            ur = grid.asCoords(i+1,j);
            ul = grid.asCoords(i-1,j);
            ub = grid.asCoords(i,j-1);
            tempCoords = grid.asCoords(i,j) - dt / grid.getDx() *
                                            (flux.evaluateFFlux(u,ur,gamma) - flux.evaluateFFlux(ul,u,gamma))
                                            - dt / grid.getDy() *
                                            (flux.evaluateGFlux(u,ut,gamma) - flux.evaluateGFlux(ub,u,gamma))
                                            + dt * source;
            //print(tempCoords);
            for (int k = 0; k < 4; ++k) {
                u1[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }
    //printf("Updating GridTemp 1 \n");
    gridTemp.update(u1);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> u2 (4,std::vector<double> (grid.getNx()*grid.getNy()));
    for(int i=0; i<grid.getNx(); i++){
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -gridTemp.rho(i,j)*g;
                source[3] = -gridTemp.v(i,j)*grid.rho(i,j)*g;
            }
            u = gridTemp.asCoords(i,j);
            ut = gridTemp.asCoords(i,j+1);
            ur = gridTemp.asCoords(i+1,j);
            ul = gridTemp.asCoords(i-1,j);
            ub = gridTemp.asCoords(i,j-1);
            tempCoords = 0.75 * grid.asCoords(i,j) + 0.25 * (gridTemp.asCoords(i,j)
                                                             - dt / grid.getDx() *
                                                               (flux.evaluateFFlux(u,ur,gamma) - flux.evaluateFFlux(ul,u,gamma))
                                                             - dt / grid.getDy() *
                                                               (flux.evaluateGFlux(u,ut,gamma) - flux.evaluateGFlux(ub,u,gamma))
                                                             + dt * source);
            for (int k = 0; k < 4; ++k) {
                u2[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }

    //printf("Updating GridTemp 2 \n");
    gridTemp.update(u2);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> unp1 (4,std::vector<double> (grid.getNx()*grid.getNy()));
    for(int i = 0; i < grid.getNx(); i++){
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -gridTemp.rho(i,j)*g;
                source[3] = -gridTemp.v(i,j)*grid.rho(i,j)*g;
            }
            u = gridTemp.asCoords(i,j);
            ut = gridTemp.asCoords(i,j+1);
            ur = gridTemp.asCoords(i+1,j);
            ul = gridTemp.asCoords(i-1,j);
            ub = gridTemp.asCoords(i,j-1);
            tempCoords = 1.0/3 * grid.asCoords(i,j) + 2.0/3 * (gridTemp.asCoords(i,j)
                                                             - dt / grid.getDx() *
                                                               (flux.evaluateFFlux(u,ur,gamma) - flux.evaluateFFlux(ul,u,gamma))
                                                             - dt / grid.getDy() *
                                                               (flux.evaluateGFlux(u,ut,gamma) - flux.evaluateGFlux(ub,u,gamma))
                                                             + dt * source);
            for (int k = 0; k < 4; ++k) {
                unp1[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }
    //printf("Updating Grid \n");
    grid.update(unp1);
    flux.updateSmax(grid.sMax());
    grid.applyBound();
    return grid;
}

//Avancement en temps 2D SSPRK3 WENO HLL
Grid& Time::evolveRK3WENO(Grid &grid, Flux &flux) {
    double gamma = grid.getGamma();
    double g = 1;
    Grid gridTemp = grid;
    flux.updateSmax(grid.sMax());
    std::vector<std::vector<double>> u1 (4,std::vector<double> (grid.getNx()*grid.getNy())),
    fluxes,fluxes2,fluxes3;
    std::vector<double> tempCoords(4);
    std::vector<double> source(4,0);
    for(int i = 0; i < grid.getNx(); i++) {
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -grid.rho(i,j)*g;
                source[3] = -grid.v(i,j)*grid.rho(i,j)*g;
            }
            fluxes = flux.WENO(grid,i,j,2);
            fluxes2 = flux.WENO(grid,i-1,j,0);
            fluxes3 = flux.WENO(grid,i,j-1,1);

            //Prémiere ètape RK3
            tempCoords = grid.asCoords(i,j) - dt / grid.getDx() *
                                            (flux.evaluateFFlux(fluxes[0],fluxes[1],gamma)
                                             - flux.evaluateFFlux(fluxes2[0],fluxes2[1],gamma))
                                             - dt / grid.getDy() *
                                             (flux.evaluateGFlux(fluxes[2],fluxes[3],gamma)
                                              - flux.evaluateGFlux(fluxes3[0],fluxes3[1],gamma))
                                              + dt * source;
            //print(tempCoords);
            for (int k = 0; k < 4; ++k) {
                u1[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }
    //Mise à jour du maillage
    gridTemp.update(u1);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> u2 (4,std::vector<double> (grid.getNx()*grid.getNy()));
    for(int i=0; i<grid.getNx(); i++){
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -gridTemp.rho(i,j)*g;
                source[3] = -gridTemp.v(i,j)*grid.rho(i,j)*g;
            }
            fluxes = flux.WENO(gridTemp,i,j,2);
            fluxes2 = flux.WENO(gridTemp,i-1,j,0);
            fluxes3 = flux.WENO(gridTemp,i,j-1,1);
            //Deuxiéme ètape RK3
            tempCoords = 0.75 * grid.asCoords(i,j) + 0.25 * (gridTemp.asCoords(i,j)
                                                   - dt / grid.getDx() *
                                                     (flux.evaluateFFlux(fluxes[0],fluxes[1],gamma)
                                                      - flux.evaluateFFlux(fluxes2[0],fluxes2[1],gamma))
                                                   - dt / grid.getDy() *
                                                     (flux.evaluateGFlux(fluxes[2],fluxes[3],gamma)
                                                      - flux.evaluateGFlux(fluxes3[0],fluxes3[1],gamma))
                                                  + dt * source);
            for (int k = 0; k < 4; ++k) {
                u2[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }

    //Mise à jour du maillage temporanée
    gridTemp.update(u2);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> unp1 (4,std::vector<double> (grid.getNx()*grid.getNy()));
    for(int i = 0; i < grid.getNx(); i++){
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -gridTemp.rho(i,j)*g;
                source[3] = -gridTemp.v(i,j)*grid.rho(i,j)*g;
            }
            fluxes = flux.WENO(gridTemp,i,j,2);
            fluxes2 = flux.WENO(gridTemp,i-1,j,0);
            fluxes3 = flux.WENO(gridTemp,i,j-1,1);

            //Troisiéme étape RK3
            tempCoords = 1.0/3 * grid.asCoords(i,j) + 2.0/3 * (gridTemp.asCoords(i,j)
                                                             - dt / grid.getDx() *
                                                             (flux.evaluateFFlux(fluxes[0],fluxes[1],gamma)
                                                             - flux.evaluateFFlux(fluxes2[0],fluxes2[1],gamma))
                                                             - dt/grid.getDy() *
                                                             (flux.evaluateGFlux(fluxes[2],fluxes[3],gamma)
                                                             - flux.evaluateGFlux(fluxes3[0],fluxes3[1],gamma))
                                                             + dt * source);
            for (int k = 0; k < 4; ++k) {
                unp1[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }

    //Mise à jour finale
    grid.update(unp1);
    flux.updateSmax(grid.sMax());
    grid.applyBound();
    return grid;
}

//Avancement en temps 2D SSPRK3 WENO HLLC
Grid& Time::evolveRK3WENOHLLC(Grid &grid, Flux &flux) {
    double gamma = grid.getGamma();
    double g = 1;
    Grid gridTemp = grid;
    flux.updateSmax(grid.sMax());
    std::vector<std::vector<double>> u1 (4,std::vector<double> (grid.getNx()*grid.getNy())),
            fluxes,fluxes2,fluxes3;
    std::vector<double> tempCoords(4);
    std::vector<double> source(4,0);
    for(int i = 0; i < grid.getNx(); i++) {
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -grid.rho(i,j)*g;
                source[3] = -grid.v(i,j)*grid.rho(i,j)*g;
            }
            fluxes = flux.WENO(grid,i,j,2); //fl[0] = u-(i+1/2)  fl[1] = u+(i+1/2)
            fluxes2 = flux.WENO(grid,i-1,j,0);
            fluxes3 = flux.WENO(grid,i,j-1,1);
            tempCoords = grid.asCoords(i,j) - dt / grid.getDx() *
                                                (flux.evaluateHLLCF(fluxes[0],fluxes[1],gamma) //F(i+1)
                                                - flux.evaluateHLLCF(fluxes2[0],fluxes2[1],gamma))
                                                - dt / grid.getDy() *
                                                (flux.evaluateHLLCG(fluxes[2],fluxes[3],gamma)
                                                - flux.evaluateHLLCG(fluxes3[0],fluxes3[1],gamma))
                                                + dt * source;
            //print(tempCoords);
            for (int k = 0; k < 4; ++k) {
                u1[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }
    //printf("Updating GridTemp 1 \n");
    gridTemp.update(u1);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> u2 (4,std::vector<double> (grid.getNx()*grid.getNy()));
    for(int i=0; i<grid.getNx(); i++){
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -gridTemp.rho(i,j)*g;
                source[3] = -gridTemp.v(i,j)*grid.rho(i,j)*g;
            }
            fluxes = flux.WENO(gridTemp,i,j,2);
            fluxes2 = flux.WENO(gridTemp,i-1,j,0);
            fluxes3 = flux.WENO(gridTemp,i,j-1,1);
            tempCoords = 0.75 * grid.asCoords(i,j) + 0.25 * (gridTemp.asCoords(i,j)
                                                            - dt / grid.getDx() *
                                                            (flux.evaluateHLLCF(fluxes[0],fluxes[1],gamma)
                                                            - flux.evaluateHLLCF(fluxes2[0],fluxes2[1],gamma))
                                                            - dt / grid.getDy() *
                                                            (flux.evaluateHLLCG(fluxes[2],fluxes[3],gamma)
                                                            - flux.evaluateHLLCG(fluxes3[0],fluxes3[1],gamma))
                                                            + dt * source);
            for (int k = 0; k < 4; ++k) {
                u2[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }

    //printf("Updating GridTemp 2 \n");
    gridTemp.update(u2);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> unp1 (4,std::vector<double> (grid.getNx()*grid.getNy()));
    for(int i = 0; i < grid.getNx(); i++){
        for (int j = 0; j < grid.getNy(); j++) {
            if(grid.getProb()=="RTI"){
                source[2] = -gridTemp.rho(i,j)*g;
                source[3] = -gridTemp.v(i,j)*grid.rho(i,j)*g;
            }
            fluxes = flux.WENO(gridTemp,i,j,2);
            fluxes2 = flux.WENO(gridTemp,i-1,j,0);
            fluxes3 = flux.WENO(gridTemp,i,j-1,1);
            tempCoords = 1.0/3 * grid.asCoords(i,j) + 2.0/3 * (gridTemp.asCoords(i,j)
                                                            - dt / grid.getDx() *
                                                            (flux.evaluateHLLCF(fluxes[0],fluxes[1],gamma)
                                                            - flux.evaluateHLLCF(fluxes2[0],fluxes2[1],gamma))
                                                            - dt / grid.getDy() *
                                                            (flux.evaluateHLLCG(fluxes[2],fluxes[3],gamma)
                                                            - flux.evaluateHLLCG(fluxes3[0],fluxes3[1],gamma))
                                                            + dt * source);
            for (int k = 0; k < 4; ++k) {
                unp1[k][i*grid.getNy() + j] = tempCoords[k];
            }
        }
    }
    //printf("Updating Grid \n");
    grid.update(unp1);
    flux.updateSmax(grid.sMax());
    grid.applyBound();
    return grid;
}

//Avancement en temps 2D SSPRK3 WENO Equation Transport
Grid1D& Time::evolveRK3Transp(Grid1D &grid, Flux &flux) {
    int dx = grid.getDx();
    Grid1D gridTemp = grid;
    flux.updateSmax(grid.sMax());
    //Grid1D grid2(grid.getDx(), grid.getNx(),grid.getBound());
    std::vector<std::vector<double>> u1 (3,std::vector<double> (grid.getNx(),0));
    std::vector<double> tempCoords(3);
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(grid, j);
        auto fluxes2 = flux.WENO1D(grid, j-1);
        tempCoords = grid.asCoords(j)
                - dt / grid.getDx() * (fluxes.first - fluxes2.first);
        //print(tempCoords);
        for (int k = 0; k < 1; ++k) {
            u1[k][j] = tempCoords[k];
        }
    }
    //printf("Updating GridTemp 1 \n");
    gridTemp.update(u1);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> u2 (3,std::vector<double> (grid.getNx(),0));
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(gridTemp, j);
        auto fluxes2 = flux.WENO1D(gridTemp, j-1);
        tempCoords = 0.75 * grid.asCoords(j) + 0.25 * (gridTemp.asCoords(j)
                - dt/grid.getDx() * (fluxes.first - fluxes2.first));
        //print(tempCoords);
        for (int k = 0; k < 1; ++k) {
            u2[k][j] = tempCoords[k];
        }
    }

    //printf("Updating GridTemp 2 \n");
    gridTemp.update(u2);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> unp1 (3,std::vector<double> (grid.getNx(),0));
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(gridTemp, j);
        auto fluxes2 = flux.WENO1D(gridTemp, j-1);
        tempCoords = 1.0/3 * grid.asCoords(j) + 2.0/3 * (gridTemp.asCoords(j)
                - dt/grid.getDx() * (fluxes.first  - fluxes2.first ));
        for (int k = 0; k < 1; ++k) {
            unp1[k][j] = tempCoords[k];
        }
    }
    //printf("Updating Grid \n");
    grid.update(unp1);
    flux.updateSmax(grid.sMax());
    grid.applyBound();
    return grid;
}

//Test pour verifier ordre SSPRK3
Grid1D &Time::testRK3(Grid1D &grid, Flux &flux) {
    int dx = grid.getDx();
    Grid1D gridTemp = grid;
    flux.updateSmax(grid.sMax());
    //Grid1D grid2(grid.getDx(), grid.getNx(),grid.getBound());
    std::vector<std::vector<double>> u1 (3,std::vector<double> (grid.getNx(),0));
    std::vector<double> tempCoords(3);
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(grid, j);
        auto fluxes2 = flux.WENO1D(grid, j-1);
        tempCoords = grid.asCoords(j)
                     - dt / grid.getDx() * (fluxes.first - fluxes2.first);
        //print(tempCoords);
        for (int k = 0; k < 1; ++k) {
            u1[k][j] = tempCoords[k];
        }
    }
    //printf("Updating GridTemp 1 \n");
    gridTemp.update(u1);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> u2 (3,std::vector<double> (grid.getNx(),0));
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(gridTemp, j);
        auto fluxes2 = flux.WENO1D(gridTemp, j-1);
        tempCoords = 0.75 * grid.asCoords(j) + 0.25 * (gridTemp.asCoords(j)
                                                       - dt/grid.getDx() * (fluxes.first - fluxes2.first));
        //print(tempCoords);
        for (int k = 0; k < 1; ++k) {
            u2[k][j] = tempCoords[k];
        }
    }

    //printf("Updating GridTemp 2 \n");
    gridTemp.update(u2);
    flux.updateSmax(gridTemp.sMax());
    gridTemp.applyBound();

    std::vector<std::vector<double>> unp1 (3,std::vector<double> (grid.getNx(),0));
    for (int j = 0; j < grid.getNx(); j++) {
        auto fluxes = flux.WENO1D(gridTemp, j);
        auto fluxes2 = flux.WENO1D(gridTemp, j-1);
        tempCoords = 1.0/3 * grid.asCoords(j) + 2.0/3 * (gridTemp.asCoords(j)
                                                         - dt/grid.getDx() * (fluxes.first  - fluxes2.first ));
        for (int k = 0; k < 1; ++k) {
            unp1[k][j] = tempCoords[k];
        }
    }
    //printf("Updating Grid \n");
    grid.update(unp1);
    flux.updateSmax(grid.sMax());
    grid.applyBound();
    return grid;
}
