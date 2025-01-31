//
// Created by UTENTE on 09/12/2024.
//

#include <valarray>
#include "Flux.h"
#include "Time.h"
#include "Print.h"
#include "Grid1D.h"
#include "helpers.h"

#ifndef PROJET_TEST_H
#define PROJET_TEST_H

#endif //PROJET_TEST_H

void test1(){
    int nx = 50; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.002; // Passo temporale
    double T = 0.25;
    Grid1D grid(nx, dx, {{1,0,1},{0.125,0,0.1}});
    grid.initialize("Sod1"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1D(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test1_1(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.002; // Passo temporale
    double T = 0.25;
    Grid1D grid(nx, dx, {{1,0,1},{0.125,0,0.1}});
    grid.initialize("Sod1"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1DLF(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test1_2(){
    int nx = 50; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.002; // Passo temporale
    double T = 0.25;
    Grid1D grid(nx, dx, {{1,0,1},{0.125,0,0.1}});
    grid.initialize("Sod1"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1DHLLC(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test2(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.2;
    Grid1D grid(nx, dx, {{1,0.75,1},{0.125,0,0.1}});
    grid.initialize("Sod2"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1D(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test2_1(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.2;
    Grid1D grid(nx, dx, {{1,0.75,1},{0.125,0,0.1}});
    grid.initialize("Sod2"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1DLF(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test2_2(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.2;
    Grid1D grid(nx, dx, {{1,0.75,1},{0.125,0,0.1}});
    grid.initialize("Sod2"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1DHLLC(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test3(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.15;
    Grid1D grid(nx, dx, {{1,-2,0.4},{1,2,.4}});
    grid.initialize("123"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1D(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test3_1(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.15;
    Grid1D grid(nx, dx, {{1,-2,0.4},{1,2,.4}});
    grid.initialize("123"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1DLF(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test3_2(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.15;
    Grid1D grid(nx, dx, {{1,-2,0.4},{1,2,.4}});
    grid.initialize("123"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1DHLLC(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
}

void test4(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dt = 0.001; // Passo temporale
    double T = 0.25;
    Grid1D grid(nx, dx, {{1,0,1},{0.125,0,0.1}});
    grid.initialize("Sod1");
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
}

void test5(){
    int nx = 200; // Nombre de celles en x
    double dx = 1.0/nx; // Passé spatial en x
    double dt = 0.001; // Passé temporal
    double T = 0.2;
    Grid1D grid(nx, dx, {{1,0.75,1},{0.125,0,0.1}});
    grid.initialize("Sod2"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test6(){
    int nx = 100; // Numero di celle in x
    double dx = 1.0/nx; // Spaziatura in x
    double dt = 0.001; // Passo temporale
    double T = 0.15;
    Grid1D grid(nx, dx, {{1,-2,0.4},{1,2,0.4}});
    grid.initialize("123"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<10; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1D(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    for (int i = 10; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test6_1(){
    int nx = 100; // Numero di celle in x
    double dx = 2.0/nx; // Spaziatura in x
    double dt = 0.01; // Passo temporale
    double T = 2;
    Grid1D grid(nx, dx, {{1,-2,0.4},{1,2,0.4}});
    grid.initialize("Smooth"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1D(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test7(){
    int nx = 50; // Numero di celle in x
    int ny = 50; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit = 300;
    std::string prob = "Case3";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance(flux,grid);
        time.isConserved(grid.integrate(),residue0);
        //print.exportToVTK(grid, "vtk/results" + std::to_string(i) + ".vtk");
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test7_1(){
    int nx = 200; // Numero di celle in x
    int ny = 200; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit = 300;
    std::string prob = "Case3";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
        //print.exportToVTK(grid, "vtk/results" + std::to_string(i) + ".vtk");
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}


void test8(){
    int nx = 100; // Numero di celle in x
    int ny = 100; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit =  300;
    std::string prob = "Case12";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance(flux,grid);
        time.isConserved(grid.integrate(),residue0);
        //print.exportToVTK(grid, "vtk/results" + std::to_string(i) + ".vtk");
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test8_1(){
    int nx = 200; // Numero di celle in x
    int ny = 200; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit =  250;
    std::string prob = "Case12";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
        //print.exportToVTK(grid, "vtk/results" + std::to_string(i) + ".vtk");
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test9_1(){
    int nx = 200; // Numero di celle in x
    int ny = 200; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit =  200;
    std::string prob = "Case15";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
        //print.exportToVTK(grid, "vtk/results" + std::to_string(i) + ".vtk");
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test9(){
    int nx = 50; // Numero di celle in x
    int ny = 50; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit = 200;
    std::string prob = "Case15";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,ny,nx,dt);
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        grid = time.advance(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test10(){
    int nx = 30; // Numero di celle in x
    int ny = 10*nx; // Numero di celle in y
    double dx = 0.25/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit = 2000;
    std::string prob = "RTI";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENOHLLC(grid,flux);
        time.isConserved(grid.integrate(),residue0);

        if(std::fmod(i,100)==0){
            print.exportToVTK(grid, "vtk/results" + std::to_string(i/100) + ".vtk");
        }
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test11(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.002; // Passo temporale
    int nit = 125;
    Grid1D grid(nx, dx, {{1,0,1},{0.125,0,0.1}});
    grid.initialize("Sod1"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    std::pair<std::vector<double>,std::vector<double>> us;
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<nit/*T/dt*/; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        grid = time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    std::cout << "Final t: " << nit*dt <<std::endl;
    //us = flux.WENO1D(grid,5);
    //print(us.first);
    //print(us.second);
    Print print;
    print.printGridToFiles1D(grid);
}

void test12(){
    int nx = 40; // Numero di celle in x
    int ny = 40; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.005; // Passo temporale
    int nit = 100;
    std::string prob = "Smooth";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print printS;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i + 1 << std::endl;
        grid.setCurrdT((i+1)*dt);
        //time.evolveRK3WENO(grid,flux);
        time.advance(flux,grid);
        time.isConserved(grid.integrate(),residue0);
        //print.exportToVTK(grid, "vtk/results" + std::to_string(i) + ".vtk");
    }
    std::string probT = "SmoothT";
    print(grid.checkSmooth());
    printS.printGridToFiles(grid);
    printS.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test13(){
    int nx = 20; // Numero di celle in x
    int ny = 20; // Numero di celle in y
    double dx = 2.0/nx; // Spaziatura in x
    double dy = 2.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit = 100;
    std::string prob = "Smooth";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print printS;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i + 1 << std::endl;
        grid.setCurrdT((i+1)*dt);
        time.evolveRK3WENOHLLC(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    print(grid.checkSmooth());
    printS.printGridToFiles(grid);
    printS.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test14(){
    int nx = 20; // Numero di celle in x
    int ny = 20; // Numero di celle in y
    double dx = 2.0/nx; // Spaziatura in x
    double dy = 2.0/ny; // Spaziatura in y
    std::string prob = "SmoothT";
    Grid grid(nx, ny, dx, dy,prob);
    Print printS;
    grid.setCurrdT(0.1);
    grid.initialize();
    printS.exportToVTK(grid, "resultsExact.vtk");
}



void test15(){
    int n = 80;
    double dx = 1.0/n;
    Flux flux(1,1,1,1,1);
    std::vector<std::vector<double>> u(5,std::vector<double> (1));
    std::vector<double> uweno;
    for(int i=0; i<5; i++){
        //u[i][0] = n*(polynomial((i-2+0.5)*dx + 1,1) - polynomial((i-2-0.5)*dx + 1,1));
        u[i][0] = -n*(std::cos((i-2+0.5)*dx) - std::cos((i-2-0.5)*dx));
        std::cout << "Mean value: " << u[i][0] << std::endl;
    }
    uweno = flux.WENO5(u);
    std::cout << "Reconstructed: " << uweno[0] << std::endl;
    //std::cout << "Exact: "<< polynomial(dx/2 + 1,0) << std::endl;
    std::cout << "Error: "<< sin(dx/2) - uweno[0] << std::endl;

}

void test16(){
    int nx = 20; // Numero di celle in x
    int ny = 20; // Numero di celle in y
    double dx = 2.0/nx; // Spaziatura in x
    double dy = 2.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit = 100;
    std::string prob = "Smooth";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print printS;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i + 1 << std::endl;
        grid.setCurrdT((i+1)*dt);
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    print(grid.checkSmooth());
    printS.printGridToFiles(grid);
    printS.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test17(){
    int nx = 1024; // Numero di celle in x
    // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dt = 0.1*pow(dx,5.0/3); // Passo temporale
    //int nit = 1000;
    double t = 0.1;
    std::string prob = "Smooth";
    Grid1D grid(nx, dx,{{0,0,0},{0,0,0}});
    grid.initialize(prob); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    Print printS;
    for (int i = 0; i<t/dt ; i++) {
        std::cout << "Iteration: " << i + 1 << " Time: " << dt*(i) << std::endl;
        //grid.setCurrdT((i+1)*dt);
        time.evolveRK3Transp(grid,flux);
        //time.isConserved(grid.integrate(),residue0);
    }
    time.isConserved(grid.integrate(),residue0);
    //print(grid.checkSmooth());
    printS.printGridToFiles1D(grid);
    //printS.exportToVTK(grid, "results.vtk");
    std::cout<< "dt: " << dt << std::endl;
}
