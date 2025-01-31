//
// Created by UTENTE on 19/12/2024.
//
#include <valarray>
#include <fstream>
#include "helpers.h"

std::vector<double> operator * (std::vector<double> vect, double mu){
    std::vector<double> out(vect.size());
    for(int i = 0; i<vect.size(); i++){
        out[i] = vect[i]*mu;
    }
    return out;
}

std::vector<double> operator * (double mu, std::vector<double> vect){
    std::vector<double> out(vect.size());
    for(int i = 0; i<vect.size(); i++){
        out[i] = vect[i]*mu;
    }
    return out;
}

std::vector<double> operator + (std::vector<double> vect1, std::vector<double> vect){
    if(vect1.size() == vect.size()){
        std::vector<double> out(vect.size());
        for(int i = 0; i<vect.size(); i++){
            out[i] = vect[i]+vect1[i];
        }
        return out;
    } else {
        std::vector<double> ciao;
        std::cerr << "Errore: qualcosa è andato storto!" << std::endl;
        return ciao;
    }

}

std::vector<double> operator - (std::vector<double> vect1, std::vector<double> vect){
    std::vector<double> out(vect.size());
    for(int i = 0; i<vect.size(); i++){
        out[i] = vect1[i]-vect[i];
    }
    return out;
}

void print(std::vector<double> vect){
    for(int i=0; i<vect.size(); i++){
        std::cout << vect[i] << "    ";
    };
    std::cout << std::endl;
}

double max(double a, double b) {
    if(a>b){
        return a;
    } else {
        return b;
    }
}

double min(double a, double b) {
    if(a<b){
        return a;
    } else {
        return b;
    }
}

std::vector<double> operator / (std::vector<double> vect, std::vector<double> vect2){
    std::vector<double> out (vect.size());
    for(int i=0; i<vect.size(); i++){
        out[i] = vect[i]/vect2[i];
    }
}



std::vector<double> operator ^ (std::vector<double> vect, int n){
    std::vector<double> out(vect.size());
    for(int i=0; i<vect.size(); i++){
        out[i] = pow(vect[i],n);
    }
    return out;
}

std::vector<double> operator * (std::vector<double> vect, std::vector<double> vect2){
    std::vector<double> out (vect.size());
    for(int i=0; i<vect.size(); i++){
        out[i] = vect[i]*vect2[i];
    }
}

double superBee(double r){
    return max(0,r);
}

std::vector<double> superBee(std::vector<double> & r){
    std::vector<double> out(r.size());
    for(int i=0; i<r.size(); i++){
        out[i] = superBee(r[i]);
    }
    return out;
}

double mediaIntegrale(double a, double b, double c, double d, double dt) {
    double integral = 0.0;
    integral = (b - a) * (d - c) - 0.2 / (M_PI * M_PI) * (
            (std::sin(M_PI * (b + d - 0.5 * dt)) - std::sin(M_PI * (a + d - 0.5 * dt))) -
            (std::sin(M_PI * (b + c - 0.5 * dt)) - std::sin(M_PI * (a + c - 0.5 * dt)))
    );
    /*
    // Somma sui rettangoli
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1000; ++j) {
            x = a + (i + 1) * dx; // Punto centrale del rettangolo in x
            y = c + (j + 1) * dy; // Punto centrale del rettangolo in y
            integral += (1 + 0.2 * sin(M_PI * ((x) + (y) - dt*(.5)))) * dx * dy;
        }
    }
    */
    // Calcolo della media integrale
    double area = (b - a) * (d - c);
    return integral / area;
}

double polynomial(double x, int mode){
    int n = 6;
    if(mode == 1){
        return pow(x,6)/6+pow(x,5)/5+pow(x,3)/3+pow(x,2)/2;
    } else {
        return pow(x,5)+pow(x,4)+pow(x,2)+pow(x,2)/2+x;
    }
}


using Real = double; // Tipo di dato per precisione

constexpr int R = 3;   // WENO r
constexpr int CN = 2;  // posizione centrale nello stencil
constexpr int DIMU = 10; // Dimensione del dominio (modificabile)

// Funzione per calcolare il quadrato
inline Real Square(const Real x) {
    return x * x;
}

void myWENO5(std::vector<std::vector<double>>& F, std::vector<double>& Fhat) {
    Real omega[R];  // pesi
    Real q[R];      // q vettori
    Real IS[R];     // misure di regolarità
    Real alpha[R];
    const Real C[R] = {1.0 / 10.0, 6.0 / 10.0, 3.0 / 10.0}; // Coefficienti WENO
    const Real epsilon = 1.0e-6; // Piccolo valore per evitare divisioni per zero

    for (int r = 0; r < DIMU; ++r) {
        // Calcolo delle misure di regolarità
        IS[0] = (13.0 / 12.0) * Square(F[CN-2][r] - 2.0 * F[CN-1][r] + F[CN][r]) +
                (1.0 / 4.0) * Square(F[CN-2][r] - 4.0 * F[CN-1][r] + 3.0 * F[CN][r]);
        IS[1] = (13.0 / 12.0) * Square(F[CN-1][r] - 2.0 * F[CN][r] + F[CN+1][r]) +
                (1.0 / 4.0) * Square(F[CN-1][r] - F[CN+1][r]);
        IS[2] = (13.0 / 12.0) * Square(F[CN][r] - 2.0 * F[CN+1][r] + F[CN+2][r]) +
                (1.0 / 4.0) * Square(3.0 * F[CN][r] - 4.0 * F[CN+1][r] + F[CN+2][r]);

        // Calcolo dei pesi non normalizzati
        alpha[0] = C[0] / Square(epsilon + IS[0]);
        alpha[1] = C[1] / Square(epsilon + IS[1]);
        alpha[2] = C[2] / Square(epsilon + IS[2]);

        // Normalizzazione dei pesi
        Real alpha_sum = alpha[0] + alpha[1] + alpha[2];
        omega[0] = alpha[0] / alpha_sum;
        omega[1] = alpha[1] / alpha_sum;
        omega[2] = alpha[2] / alpha_sum;

        // Calcolo di q
        q[0] = (1.0 / 6.0) * (2.0 * F[CN-2][r] - 7.0 * F[CN-1][r] + 11.0 * F[CN][r]);
        q[1] = (1.0 / 6.0) * (-F[CN-1][r] + 5.0 * F[CN][r] + 2.0 * F[CN+1][r]);
        q[2] = (1.0 / 6.0) * (2.0 * F[CN][r] + 5.0 * F[CN+1][r] - F[CN+2][r]);

        // Ricostruzione del valore finale
        Fhat[r] = omega[0] * q[0] + omega[1] * q[1] + omega[2] * q[2];
    }
}

double square(double x){
    return x*x;
}

double leggiNumeroDaMatrice(const std::string& nomeFile, int i, int j, int n) {
    std::ifstream file(nomeFile);
    if (!file.is_open()) {
        throw std::runtime_error("Errore: impossibile aprire il file " + nomeFile);
    }

    // Indice da raggiungere
    int indiceTarget = i * n + j;

    // Lettura sequenziale del file
    int valore;
    int indiceCorrente = 0;
    while (file >> valore) {
        if (indiceCorrente == indiceTarget) {
            return valore;
        }
        indiceCorrente++;
    }

    // Se esce dal ciclo, significa che non ha trovato l'indice
    throw std::out_of_range("Errore: indice fuori dai limiti della matrice");
}