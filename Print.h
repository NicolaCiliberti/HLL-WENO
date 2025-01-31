//
// Created by UTENTE on 20/12/2024.
//

#ifndef PROJET_PRINT_H
#define PROJET_PRINT_H

#endif //PROJET_PRINT_H

#include "Grid.h"
#include "Grid1D.h"

class Print {
public:
    void writeVariableToFile(Grid & grid, const std::string& filename, const std::string& variable);
    void writeVariableToFile1D(Grid1D & grid, const std::string& filename, const std::string& variable);
    void printGridToFiles(Grid & grid);
    void printGridToFiles1D(Grid1D & grid);
    void exportToVTK(Grid & grid, const std::string& filename);
    void printSmoothGrid(int n);
};
