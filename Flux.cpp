//
// Created by UTENTE on 18/12/2024.
//

#include <valarray>
#include "Flux.h"

//Evaluer Flux HLL horizontal
std::vector<double> Flux::evaluateFFlux(Grid &grid, int i, int j) { //Evaluates flux F(i+1/2,j)
    double sr,sl,al,ar;
    //Vitesses son gauche et droite
    al = sqrt(grid.getGamma()*grid.p(i,j)/grid.rho(i,j));
    ar = sqrt(grid.getGamma()*grid.p(i+1,j)/grid.rho(i+1,j));
    //Vitesses onde gauche et droite
    sl = min(grid.u(i,j)-al,grid.u(i+1,j)-ar);
    sl = min(0,sl);
    sr = max(grid.u(i,j)+al,grid.u(i+1,j)+ar);
    sr = max(0,sr);

    //Calcul du fluxe
    std::vector<double> fl,fr,ul,ur;
    std::vector<double> flux(4);
    ul = grid.asCoords(i, j);
    ur = grid.asCoords(i+1, j);
    if(sl>sr){
        std::cerr << "Error: SL > SR" << std::endl;
    }
    if(sr<=0){ //Flux FR
        return evaluateF(ur, grid.getGamma());
    } else if (sl>=0){ //Flux FL
        return evaluateF(ul, grid.getGamma());
    } else { //Flux FHLL
        fl = evaluateF(ul, grid.getGamma());
        fr = evaluateF(ur, grid.getGamma());
        flux = (1/(sr-sl))*(sr * fl - sl * fr + sl * sr * (ur - ul));

    }
    return flux;
}

//Evaluer Flux HLL horizontal V2 (Sans utiliser Grid)
std::vector<double> Flux::evaluateFFlux(std::vector<double> & ul, std::vector<double> & ur, double gamma) { //Evaluates flux F(i+1/2,j)
    double al, ar, sl, sr, pl, pr;
    std::vector<double> fl,fr;

    pl = max(0, (gamma-1)*(ul[3] - 0.5*(ul[1]*ul[1] + ul[2]*ul[2])/ul[0]));
    pr = max(0, (gamma-1)*(ur[3] - 0.5*(ur[1]*ur[1] + ur[2]*ur[2])/ur[0]));

    al = sqrt(gamma * pl / ul[0]);  // Vitesse son gauche
    ar = sqrt(gamma * pr / ur[0]);  // Vitesse son droite

    sl = std::min(ul[1]/ul[0] - al, ur[1]/(ur[0]) - ar);// Vitesse onde gauche
    sr = std::max(ul[1]/ul[0] + al, ur[1]/(ur[0]) + ar);  // Vitesse onde droite

    std::vector<double> flux(4);
    if(sr < sl){
        std::cerr << "Warning: SR < SL" << std::endl;
    }
    if (sr <= 0) {
        //return evaluateF(grid.rho(i + 1), grid.u(i + 1), grid.p(i + 1), grid.E(i + 1));
        return evaluateF(ur,gamma);
    } else if (sl >= 0) {
        //return evaluateF(grid.rho(i), grid.u(i), grid.p(i), grid.E(i));
        return evaluateF(ul,gamma);
    } else {
        fl = evaluateF(ul,gamma);
        fr = evaluateF(ur,gamma);
        flux = (1 / (sr - sl)) * ((sr * fl) - (sl * fr) + sl * sr * (ur -  ul));
    }
    return flux;
}

//Evaluer Flux HLL vertical
std::vector<double> Flux::evaluateGFlux(std::vector<double> & ub, std::vector<double> & ut, double gamma) { //Evaluates flux F(i+1/2,j)
    double ab, at, sb, st, pb, pt;
    std::vector<double> fb,ft;


    pb = max(0, (gamma-1)*(ub[3] - 0.5*(ub[1]*ub[1] + ub[2]*ub[2])/ub[0]));
    pt = max(0, (gamma-1)*(ut[3] - 0.5*(ut[1]*ut[1] + ut[2]*ut[2])/ut[0]));

    at = sqrt(gamma * pt / ut[0]);  // Vitesse son gauche
    ab = sqrt(gamma * pb / ub[0]);  // Vitesse son droite

    /*
    printf("Rhol: %f, Rhor: %f\n",ub[0], ut[0]);
    printf("Pl: %f, Pr: %f\n",pCoords(ub), pCoords(ut));
    printf("Al: %f, Ar: %f\n",al,ar);
*/
    sb = std::min(ub[2]/(ub[0]) - ab, ut[2]/(ut[0]) - at);  // Vitesse onde gauche
    st = std::max(ub[2]/(ub[0]) + ab, ut[2]/(ut[0]) + at);  // Vitesse onde droite

    std::vector<double> flux(4);
    if(st < sb){
        std::cerr << "Warning: ST < SB" << std::endl;
    }
    if (st <= 0) {
        //retutn evaluateF(grid.rho(i + 1), grid.u(i + 1), grid.p(i + 1), grid.E(i + 1));
        return evaluateG(ut, gamma);
    } else if (sb >= 0) {
        //retuyn evaluateF(grid.rho(i), grid.u(i), grid.p(i), grid.E(i));
        return evaluateG(ub, gamma);
    } else {
        fb = evaluateG(ub, gamma);
        ft = evaluateG(ut, gamma);
        flux = (1 / (st - sb)) * ((st * fb) - (sb * ft) + (st * sb * (ut -  ub)));
    }
    return flux;
}

//Evaluer Flux HLL vertical V2 (Sans utiliser Grid)
std::vector<double> Flux::evaluateGFlux(Grid &grid, int i, int j) { //Evaluates flux G(i,j+1/2)
    double st,sb,at,ab;
    double gamma = grid.getGamma();

    ab = sqrt(gamma*grid.p(i,j)/grid.rho(i,j));
    at = sqrt(gamma*grid.p(i,j+1)/grid.rho(i,j+1));

    sb = min(grid.v(i,j)-ab,grid.v(i,j+1)-at);
    st = max(grid.v(i,j)+ab,grid.v(i,j+1)+at);

    if(sb>st){
        std::cerr << "Error: SB > ST" << std::endl;
    }
    std::vector<double> fb,ft,flux,ub,ut;

    if(st<=0){
        return evaluateG(grid.rho(i,j+1),grid.u(i,j+1),grid.v(i,j+1),grid.p(i,j+1),grid.E(i,j+1));
    } else if (sb>=0){
        return evaluateG(grid.rho(i,j),grid.u(i,j),grid.v(i,j),grid.p(i,j),grid.E(i,j));
    } else {
        fb = evaluateG(grid.rho(i,j),grid.u(i,j),grid.v(i,j),grid.p(i,j),grid.E(i,j));
        ft = evaluateG(grid.rho(i,j+1),grid.u(i,j+1),grid.v(i,j+1),grid.p(i,j+1),grid.E(i,j+1));
        ut = grid.asCoords(i,j+1);
        ub = grid.asCoords(i,j);
        flux = (st * fb - sb * ft + st * sb * (ut-ub))*(1/(st-sb));
    }
    return flux;
}

//Evaluer Flux HLLC horizontal
std::vector<double> Flux::evaluateHLLCF(std::vector<double> & ul, std::vector<double> & ur, double gamma) {
    double al, ar, sl, sr, pl, pr, sc;
    std::vector<double> fl, fr, flux;

    if(ul.size()==4){
        pl = max(0, (gamma-1)*(ul[3] - 0.5*(ul[1]*ul[1] + ul[2]*ul[2])/ul[0]));
        pr = max(0, (gamma-1)*(ur[3] - 0.5*(ur[1]*ur[1] + ur[2]*ur[2])/ur[0]));
    } else {
        pl = max(0, (gamma-1)*(ul[2] - 0.5*(ul[1]*ul[1]/ul[0])));
        pr = max(0, (gamma-1)*(ur[2] - 0.5*(ur[1]*ur[1]/ur[0])));
    }

    al = sqrt(gamma * pl / ul[0]);  // Velocità del suono sinistra
    ar = sqrt(gamma * pr / ur[0]);  // Velocità del suono destra

    sl = std::min(ul[1]/ul[0] - al, ur[1]/ur[0] - ar);  // Velocità dell'onda sinistra
    sr = std::max(ul[1]/ul[0] + al, ur[1]/ur[0] + ar);  // Velocità dell'onda destra

    // Calcolare la velocità di contatto s_c per HLLC
    sc = (pr - pl + ul[0] * sl - ur[0] * sr) / (ul[0] * (sl - ul[1]/ul[0]) - ur[0] * (sr - ur[1]/ur[0]));

    // Gestione delle condizioni al contatto
    if (sr <= 0) {
        return evaluateF(ur, gamma);
    } else if (sl >= 0) {
        return evaluateF(ul, gamma);
    } else if (sl < sc) {
        fl = evaluateF(ul, gamma);
        fr = evaluateF(ur, gamma);
        flux = (1 / (sr - sl)) * ((sr * fl) - (sl * fr) + sl * sr * (ur - ul));
    } else {
        // Calcolare il flusso usando s_c
        flux = evaluateFContact(ul, ur, sl, sr, sc, gamma);
    }

    return flux;
}

//Evaluer Flux HLLC vertical
std::vector<double> Flux::evaluateHLLCG(std::vector<double> & ub, std::vector<double> & ut, double gamma) {
    double ab, at, sb, st, pb, pt, sc;
    std::vector<double> fb, ft, flux;

    pb = max(0, (gamma-1)*(ub[3] - 0.5*(ub[1]*ub[1] + ub[2]*ub[2])/ub[0]));
    pt = max(0, (gamma-1)*(ut[3] - 0.5*(ut[1]*ut[1] + ut[2]*ut[2])/ut[0]));

    ab = sqrt(gamma * pb / ub[0]);  // Velocità del suono sinistra
    at = sqrt(gamma * pt / ut[0]);  // Velocità del suono destra

    sb = std::min(ub[2]/ub[0] - ab, ut[2]/ut[0] - at);  // Velocità dell'onda sinistra
    st = std::max(ub[2]/ub[0] + ab, ut[2]/ut[0] + at);  // Velocità dell'onda destra

    // Calcolare la velocità di contatto s_c per HLLC
    sc = (pt - pb + ub[0] * sb - ut[0] * st) / (ub[0] * (sb - ub[2]/ub[0]) - ut[0] * (st - ut[2]/ut[0]));

    // Gestione delle condizioni al contatto
    if (st <= 0) {
        return evaluateG(ut, gamma);
    } else if (sb >= 0) {
        return evaluateG(ub, gamma);
    } else if (sb < sc) {
        fb = evaluateG(ub, gamma);
        ft = evaluateG(ut, gamma);
        flux = (1 / (st - sb)) * ((st * fb) - (sb * ft) + sb * st * (ut - ub));
    } else {
        // Calcolare il flusso usando s_c
        flux = evaluateGContact(ub, ut, sb, st, sc, gamma);
    }

    return flux;
}

//Evaluer Onde de Contact horizontale pour HLLC
std::vector<double> Flux::evaluateFContact(std::vector<double> & ul, std::vector<double> & ur, double sl, double sr, double sc, double gamma) {
    // Calcolare il flusso alla superficie di contatto
    if(ul.size()==4){
        std::vector<double> flux(4);

        // Utilizzare sl, sr, sc per il calcolo del flusso
        flux[0] = (ur[0] * (sr - sc) - ul[0] * (sl - sc)) / (sr - sl);  // Densità
        flux[1] = (ur[1] * (sr - sc) - ul[1] * (sl - sc)) / (sr - sl);  // Velocità x
        flux[2] = (ur[2] * (sr - sc) - ul[2] * (sl - sc)) / (sr - sl);  // Velocità y
        flux[3] = (ur[3] * (sr - sc) - ul[3] * (sl - sc)) / (sr - sl);  // Energia
        return flux;
    } else {
        std::vector<double> flux(3);

        // Utilizzare sl, sr, sc per il calcolo del flusso
        flux[0] = (ur[0] * (sr - sc) - ul[0] * (sl - sc)) / (sr - sl);  // Densità
        flux[1] = (ur[1] * (sr - sc) - ul[1] * (sl - sc)) / (sr - sl);  // Velocità x
        flux[2] = (ur[2] * (sr - sc) - ul[2] * (sl - sc)) / (sr - sl);  // Energia
        return flux;
    }


}

//Evaluer Onde de Contact verticale pour HLLC
std::vector<double> Flux::evaluateGContact(std::vector<double> & ub, std::vector<double> & ut, double sb, double st, double sc, double gamma) {
    std::vector<double> flux(4);

    flux[0] = (ut[0] * (st - sc) - ub[0] * (sb - sc)) / (st - sb);  // Densità
    flux[1] = (ut[1] * (st - sc) - ub[1] * (sb - sc)) / (st - sb);  // Velocità x
    flux[2] = (ut[2] * (st - sc) - ub[2] * (sb - sc)) / (st - sb);  // Velocità y
    flux[3] = (ut[3] * (st - sc) - ub[3] * (sb - sc)) / (st - sb);  // Energia

    return flux;
}

std::vector<double> Flux::evaluateF(double rho, double u, double v, double p, double E) {
    std::vector<double> f(4);
    f[0] = rho*u;
    f[1] = rho*u*u + p;
    f[2] = rho*u*v;
    f[3] = u*(E + p);
    return f;
}

std::vector<double> Flux::evaluateF(std::vector<double> & coords, double gamma){
    if (coords.size()==3){
        std::vector<double> f(3);
        f[0] = coords[1];
        f[1] = coords[1]*coords[1]/coords[0] + (gamma-1)*(coords[2] - 0.5*( coords[1]*coords[1]/coords[0]));
        f[2] = coords[1]/coords[0]*(coords[2] + (gamma-1)*(coords[2] - 0.5*( coords[1]*coords[1]/coords[0])));
        return f;
    } else {
        std::vector<double> f(4);
        f[0] = coords[1];
        f[1] = coords[1]*coords[1]/coords[0] + (gamma-1)*(coords[3] - 0.5*(coords[1]*coords[1]+ coords[2]*coords[2])/coords[0]);
        f[2] = coords[1]*coords[2]/coords[0];
        f[3] = coords[1]/coords[0]*(coords[3] + (gamma-1)*(coords[3] - 0.5*( coords[1]*coords[1]/coords[0] + coords[2]*coords[2]/coords[0])));
        return f;
    }
}

std::vector<double> Flux::evaluateG(std::vector<double> & coords,double gamma){
        std::vector<double> f(4);
        f[0] = coords[2];
        f[1] = coords[1]*coords[2]/coords[0];
        f[2] = coords[2]*coords[2]/coords[0] + (gamma-1)*(coords[3] - 0.5*(coords[1]*coords[1]+ coords[2]*coords[2])/coords[0]);
        f[3] = coords[2]/coords[0]*(coords[3] + (gamma-1)*(coords[3] - 0.5*( coords[1]*coords[1]/coords[0] + coords[2]*coords[2]/coords[0])));
        return f;
}

std::vector<double> Flux::evaluateF(double rho, double u, double p, double E) {
    std::vector<double> f(3);
    f[0] = rho*u;
    f[1] = rho*u*u + p;
    f[2] = u*(E + p);
    return f;
}

std::vector<double> Flux::evaluateG(double rho, double u, double v, double p, double E) {
    std::vector<double> g(4);
    g[0] = rho*v;
    g[1] = rho*u*v;
    g[2] = rho*v*v + p;
    g[3] = v*(E + p);
    return g;
}

// Evaluate Flux HLL 1D
std::vector<double> Flux::evaluateFlux(Grid1D &grid, int i) {
    double al, ar, sl, sr;
    std::vector<double> ul,ur,fl,fr;

    al = sqrt(1.4 * grid.p(i) / grid.rho(i));  // Vitesse son gauche
    ar = sqrt(1.4 * grid.p(i + 1) / grid.rho(i + 1));  // Vitesse son droite

    sl = std::min(grid.u(i) - al, grid.u(i + 1) - ar);  // Vitesse onde gauche
    sr = std::max(grid.u(i) + al, grid.u(i + 1) + ar);  // Vitesse onde droite

    std::vector<double> flux(3);
    ul = grid.asCoords(i);
    ur = grid.asCoords(i+1);
    if(sr < sl){
        std::cerr << "Warning: SR < SL" << std::endl;
    }
    if (sr <= 0) {
        //return evaluateF(grid.rho(i + 1), grid.u(i + 1), grid.p(i + 1), grid.E(i + 1));
        return evaluateF(ur,1.4);
    } else if (sl >= 0) {
        //return evaluateF(grid.rho(i), grid.u(i), grid.p(i), grid.E(i));
        return evaluateF(ul,1.4);
    } else {
        fl = evaluateF(ul,1.4);
        fr = evaluateF(ur,1.4);
        flux = (1 / (sr - sl)) * (sr * fl - sl * fr + sl * sr * (grid.asCoords(i + 1) - grid.asCoords(i)));
    }
    return flux;
}

// Evaluate Flux HLL 1D V2 (Sans Grid1D)
std::vector<double> Flux::evaluateFlux(std::vector<double> & ul, std::vector<double> & ur) {
    double al, ar, sl, sr, pl, pr;
    std::vector<double> fl,fr;
    double gamma = 1.4;

    pl = (gamma-1)*(ul[2] - 0.5*(ul[1]*ul[1])/ul[0]);
    pr = (gamma-1)*(ur[2] - 0.5*(ur[1]*ur[1])/ur[0]);

    pl = max(0, pl);
    pr = max(0, pr);


    al = sqrt(gamma * pl / ul[0]);  // Vitesse son gauche
    ar = sqrt(gamma * pr / ur[0]);  // Vitesse son droite

    /*
    printf("Rhol: %f, Rhor: %f\n",ul[0], ur[0]);
    printf("Pl: %f, Pr: %f\n",pCoords(ul), pCoords(ur));
    printf("Al: %f, Ar: %f\n",al,ar);
*/
    sl = std::min(ul[1]/ul[0] - al, ur[1]/ur[0] - ar);  // Vitesse onde gauche
    sr = std::max(ul[1]/ul[0] + al, ur[1]/ur[0] + ar);  // Vitesse onde droite

    std::vector<double> flux(3);
    if(sr < sl){
        std::cerr << "Warning: SR < SL" << std::endl;
    }
    if (sr <= 0) {
        //return evaluateF(grid.rho(i + 1), grid.u(i + 1), grid.p(i + 1), grid.E(i + 1));
        return evaluateF(ur,gamma);
    } else if (sl >= 0) {
        //return evaluateF(grid.rho(i), grid.u(i), grid.p(i), grid.E(i));
        return evaluateF(ul,gamma);
    } else {
        fl = evaluateF(ul,gamma);
        fr = evaluateF(ur,gamma);
        flux = (1 / (sr - sl)) * (sr * fl - sl * fr + sl * sr * (ur -  ul));
        //printf("SR - SL: %f\n",sr-sl);
        //print(flux);
    }
    return flux;
}

// Evaluate Flux Lax-Friedrichs 1D
std::vector<double> Flux::evaluateLF(std::vector<double> & ul, std::vector<double> & ur, double gamma) {
    double al, ar, sl, sr,pl,pr,s_star;
    std::vector<double> fl,fr;

    pl = max(0, (gamma-1)*(ul[2] - 0.5*(ul[1]*ul[1])/ul[0]));
    pr = max(0, (gamma-1)*(ur[2] - 0.5*(ur[1]*ur[1])/ur[0]));

    al = sqrt(gamma * pl / ul[0]);  // Vitesse son gauche
    ar = sqrt(gamma * pr / ur[0]);  // Vitesse son droite

    /*
    printf("Rhol: %f, Rhor: %f\n",ul[0], ur[0]);
    printf("Pl: %f, Pr: %f\n",pCoords(ul), pCoords(ur));
    printf("Al: %f, Ar: %f\n",al,ar);
*/
    sl = std::min(ul[1]/ul[0] - al, ur[1]/ur[0] - ar);  // Vitesse onde gauche
    sr = std::max(ul[1]/ul[0] + al, ur[1]/ur[0] + ar);  // Vitesse onde droite

    std::vector<double> flux(3);
    fl = evaluateF(ul,gamma);
    fr = evaluateF(ur,gamma);
    flux = 0.5*(fl+fr) - 0.5 * smax * (ur-ul);
    return flux;
}

// Evaluate Flux Lax-Friedrichs 2D horizontal
std::vector<double> Flux::evaluateLFx(std::vector<double> & ul, std::vector<double> & ur, double gamma) {
    double al, ar, sl, sr,pl,pr,s_star;
    std::vector<double> fl,fr;

    std::vector<double> flux(4);
    fl = evaluateF(ul,gamma);
    fr = evaluateF(ur,gamma);
    flux = 0.5*(fl+fr) - 0.5 * sxmax * (ur-ul);
    return flux;
}

// Evaluate Flux Lax-Friedrichs 2D vertical
std::vector<double> Flux::evaluateLFy(std::vector<double> & ul, std::vector<double> & ur, double gamma) {
    double al, ar, sl, sr,pl,pr,s_star;
    std::vector<double> fl,fr;

    std::vector<double> flux(4);
    fl = evaluateG(ul,gamma);
    fr = evaluateG(ur,gamma);
    flux = 0.5*(fl+fr) - 0.5 * symax * (ur-ul);
    return flux;
}

// Reconstruction WENO 1D donnée la celle
std::pair<std::vector<double>, std::vector<double>> Flux::WENO1D(Grid1D &grid, int j) {

    //printf("i: %d\n",j);
    std::vector<double> ul(3), ur(3);
    std::vector<std::vector<double>> f(5, std::vector<double>(3));

    //UR - stencil {i-1,...,i+3}
    for (int l = 0; l < 5; l++) {
        f[l] = grid.asCoords(j - l + 3);
    }

    ur = WENO5(f);

    //UL - stencil {i-2,...,i+2}

    for (int l = 0; l < 5; l++) {
        f[l] = grid.asCoords(j + l - 2); //fi[0] =  i-2
    }

    ul = WENO5(f);

    return {ul, ur};
}

// Reconstruction WENO generale
std::vector<double> Flux::WENO5(std::vector<std::vector<double>> & fi) {
    const std::vector<double> dweights = {0.6, 0.1, 0.3};
    int dim = fi[0].size();
    std::vector<std::vector<double>> weights(3, std::vector<double>(dim));
    double sumw;
    std::vector<double> f0(dim),f1(dim),f2(dim),u(dim);
    std::vector<std::vector<double>> beta(3, std::vector<double>(dim));

    // Poids non linéaires
    for(int i=0; i<dim; i++) {
        beta[0][i] = (13.0 / 12) * pow((fi[1][i] - 2 * fi[2][i] + fi[3][i]),2)
                + 0.25 * pow((fi[1][i] - fi[3][i]),2);
        beta[1][i] = (13.0 / 12) * pow((fi[0][i] - 2 * fi[1][i] + fi[2][i]),2)
                + 0.25 * pow((3 * fi[2][i] - 4 * fi[1][i] + fi[0][i]),2);
        beta[2][i] = (13.0 / 12) * pow((fi[2][i] - 2 * fi[3][i] + fi[4][i]),2)
                + 0.25 * pow((3 * fi[2][i] - 4 * fi[3][i] + fi[4][i]),2);
    }

    //Calcul des poids
    for (int i = 0; i < dim; i++) {
        for (int k = 0; k < 3; k++) {
            weights[k][i] =  dweights[k]*std::pow( 1/(1e-6 + beta[k][i]), 2);
        }
    }
    //Normalisation des poids
    for (int i = 0; i < dim; i++) {
        sumw = weights[0][i] + weights[1][i] + weights[2][i];
        for (int k = 0; k < 3; k++) {
            weights[k][i] = weights[k][i] / sumw;
        }
    }

    //Polynomes interpolés
    for(int i=0; i<dim; i++){
        f0[i] = -1.0 / 6.0 * fi[1][i] + 5.0 / 6.0 * fi[2][i] + 1.0 / 3.0 * fi[3][i];
        f1[i] = 1.0 / 3.0 * fi[0][i] - 7.0 / 6.0 * fi[1][i] + 11.0 / 6.0 * fi[2][i];
        f2[i] = 1.0 / 3.0 * fi[2][i] + 5.0 / 6.0 * fi[3][i] - 1.0 / 6.0 * fi[4][i];
    }

    //Reconstruction finale
    for (int k = 0; k < dim; k++) {
        u[k] = weights[0][k] * f0[k] +  weights[1][k] * f1[k] +  weights[2][k] * f2[k];
    }
    return u;
}

// Reconstruction WENO 2D donnée la celle
std::vector<std::vector<double>> Flux::WENO(Grid &grid, int i, int j, int dir) {

    std::vector<double> ul(4), ur(4), ut(4), ub(4);
    std::vector<std::vector<double>> f(5, std::vector<double>(4));

    if(dir == 0){
        //UR
        for (int l = 0; l < 5; l++) {
            f[l] = grid.asCoords(i - l + 3,j);
        }

        ur = WENO5(f);

        //UL
        for (int l = 0; l < 5; l++) {
            f[l] = grid.asCoords(i + l - 2,j); //fi[0] =  i-2
        }

        ul = WENO5(f);
    } else if (dir==1) {
        //UT
        for (int l = 0; l < 5; l++) {
            f[l] = grid.asCoords(i,j - l + 3);
        }

        ut = WENO5(f);

        //UB
        for (int l = 0; l < 5; l++) {
            f[l] = grid.asCoords(i , j + l - 2);
        }

        ub = WENO5(f);
        return {ub, ut};
    } else {
        for (int l = 0; l < 5; l++) {
            f[l] = grid.asCoords(i - l + 3,j);
        }
        ur = WENO5(f);

        //UL
        for (int l = 0; l < 5; l++) {
            f[l] = grid.asCoords(i + l - 2,j); //fi[0] =  i-2
        }

        ul = WENO5(f);

        //UT
        for (int l = 0; l < 5; l++) {
            f[l] = grid.asCoords(i,j - l + 3);
        }

        ut = WENO5(f);

        //UB
        for (int l = 0; l < 5; l++) {
            f[l] = grid.asCoords(i, j + l - 2); //fi[0] =  i-2
        }

        ub = WENO5(f);
        return {ul,ur,ub,ut};
    }
    return {ul, ur};
}

void Flux::updateSmax(double s) {
    smax = s;
}

void Flux::setSxmax(double sxmax) {
    Flux::sxmax = sxmax;
}

void Flux::setSmax(double smax) {
    Flux::smax = smax;
}

void Flux::setSymax(double symax) {
    Flux::symax = symax;
}

