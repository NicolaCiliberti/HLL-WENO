//
// Created by UTENTE on 24/12/2024.
//

#include "Grid1D.h"

#include <iostream>
#include <valarray>



Grid1D::Grid1D(int nx, double dx, std::vector<std::vector<double>> bound)
        : nx(nx), dx(dx), bound(bound) {
    rho_.resize(nx, 0.0);
    u_.resize(nx , 0.0);
    p_.resize(nx, 0.0);
}

//Initialiser le maillage
void Grid1D::initialize(const std::string &cas) {
    casTest = cas;
    double x,rl,pl,ul,rr,pr,ur,xs;
    if(cas == "Sod1"){ //Sod classique
        for (int i = 0; i < nx; ++i) {
            x= dx*i;
            if (x < 0.5) {
                rho(i) = 1.0;
                p(i) = 1.0;
                u(i) = 0.0;
            } else {
                rho(i) = 0.125;
                p(i) = 0.1;
                u(i) = 0.0;
            }
        }
    } else if (cas == "Sod2") { //Sod modifié

        for (int i = 0; i < nx; ++i) {
            x= dx*i;
            if (x < 0.3) {
                rho(i) = 1.0;
                p(i) = 1.0;
                u(i) = 0.75;
            } else {
                rho(i) = 0.125;
                p(i) = 0.1;
                u(i) = 0.0;
            }
        }
    } else if (cas == "123") { //Noh
        for (int i = 0; i < nx; ++i) {
            x= dx*i;
            if (x < 0.5) {
                rho(i) = 1.0;
                p(i) = 0.4;
                u(i) = -2;
            } else {
                rho(i) = 1.0;
                p(i) = 0.4;
                u(i) = 2;
            }
        }
    } else if (cas == "Smooth"){ //Accuracy test
        casTest = cas;
        for (int i = 0; i < nx; ++i) {
            x= dx*i;
            rho(i) = std::exp(-300*pow(x-0.5,2));
            p(i) = 0;
            u(i) = 0;
        }
    }
    applyBound();
    intConserved = integrate();
    std::cout << "Initial conserved: " << std::endl;
    print(intConserved);
}


double Grid1D::getDx() const {
    return dx;
}

//Retourner les variables conservatifs
std::vector<double> Grid1D::asCoords(int i) {
    std::vector<double> coords(3);
    coords[0] = rho(i);
    coords[1] = rho(i)*u(i);
    coords[2] = E(i);
    return coords;
}

double Grid1D::E(int i) {
    return p(i)/(gamma-1) + 0.5*rho(i)*(u(i)*u(i));
}

double Grid1D::H(int i) {
    return (E(i) + p(i))/rho(i);
}

void Grid1D::updateRho(std::vector<double> &rhonew) {
    rho_ = rhonew;
}

void Grid1D::updateU(std::vector<double> & unew) {
    for(int i=0; i<nx; i++){
        //printf("Old u: %f\n",u_[i]);
        u_[i] = (1/rho(i))*unew[i];
        //printf("New u: %f\n",u_[i]);
    }
}

void Grid1D::updateEP(std::vector<double> &Enew) {
    //printf("Old E: \n");
    //print(E_);
    E_ = Enew;
    //printf("New E: \n");
    //print(Enew);
    double tempP;
    for (int i=0; i<nx; i++){
        //printf("New p: %f\n",(gamma-1)*(E_[i] - 0.5*rho(i)*(u(i)*u(i))));
        tempP = (gamma-1)*(E_[i] - 0.5*rho(i)*(u(i)*u(i)));
        p_[i] = max(tempP,0);
    }
}

void Grid1D::update(std::vector<std::vector<double>> &unew) {
    updateRho(unew[0]);
    updateU(unew[1]);
    updateEP(unew[2]);
}

double &Grid1D::rho(int i) {
    if (i < 0 || i > nx-1) {
        if(casTest=="Smooth"){
            if(i<0){
                return rho(nx+i);
            } else {
                return rho(i-nx);
            }
        }
        return ghostCell(i, 0);
    } else {
        return rho_[i];
    }
}

double &Grid1D::u(int i) {
    if (i < 0 || i > nx-1) {
        if(casTest=="Smooth"){
            return u(0);
        }
        return ghostCell(i, 1);
    } else {
        return u_[i];
    }
}

double &Grid1D::p(int i) {
    if (i < 0 || i > nx-1) {
        if(casTest=="Smooth"){
            return p(0);
        }
        return ghostCell(i,2);
    } else {
        return p_[i];
    }
}

//Integrer les variables conservatifs
std::vector<double> Grid1D::integrate() {
    std::vector<double> sums(3,0);
    for(int i=0; i<nx; i++){
        sums[0] += dx*rho(i);
        sums[1] += dx*rho(i)*u(i);
        sums[2] += dx*E(i);
    }
    return sums;
}

const std::vector<double> &Grid1D::getIntConserved() const {
    return intConserved;
}

double &Grid1D::ghostCell(int i, int var) {
    if(i < 0){
        return bound[0][var];
    } else if (i > nx-1 ) {
        return bound[1][var];
    } else {
        std::cerr << "Chiamata errata Ghost Cell" << std::endl;
        return bound[0][0];
    }
}

//Conditions de bord
void Grid1D::applyBound() {
    if(casTest != "Smooth" ){
        bound[0] = {rho(0),u(0),p(0)};
        bound[1] = {rho(nx-1),u(nx-1),p(nx-1)};
    }
     //bound[0] = {1,0,1};
    //bound[1] = {0.125,0,0.1};
}

const std::vector<std::vector<double>> &Grid1D::getBound() const {
    return bound;
}

double Grid1D::sMax() {
    double max = 0;
    for(int i=0; i<nx; i++){
        if(std::abs(u(i)) + sqrt(gamma*p(i)/rho(i)) > max) {
            max = std::abs(u(i)) + sqrt(gamma*p(i)/rho(i));
        }
    }
    return max;
}

//Comparaison avec solution exacte
std::vector<double> Grid1D::checkSmooth() {
    std::vector<double> out(3,0);
    double x;
    for(int i = 0; i < nx; i++) {
        x = i * (dx) * (3.0/2);
        out[0] += dx*std::abs(rho(i) - exp(-300*pow(x - 0.5 - 0.1,2)));
        std::cout << rho(i) << " "<< exp(-300*pow(x - 0.5 - 0.1,2)) << std::endl;
        out[1] += dx*std::abs(u(i) - 1);
        out[2] += dx*std::abs(p(i) - 1);
    }
    return out;
}





