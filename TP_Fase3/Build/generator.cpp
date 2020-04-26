#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include "generator.h"
//#define M_PI 3.14159265358979323846
std::ofstream outFile;



bool exists(const std::string fileName) {
    std::fstream checkFile(fileName);
    return checkFile.is_open();
}

bool open(const std::string strFile) {

    if (exists(strFile)) {
        std::cout << "There's already a file with this name. It will be replaced\n";
    }

    outFile.open(strFile);
    if (outFile.is_open()) {

        return true;
    }

    return false;
}

void close() {
    if (outFile.is_open()) {
        outFile.close();
    }
    else {
        std::cout << "Error 404.\n";
    }

}

void writeLine(const std::string line) {
    if (outFile.is_open()) {
        outFile << line << "\n";

    }
    else {
        std::cout << "File Already close.\n";
    }
}


void esfera(char* argv[]) {
    double pi = 3.14159265358979323846;
   // double x = M_PI;

    int radius = atoi(argv[2]);
    int slices = atoi(argv[3]);
    int stacks = atoi(argv[4]);
    std::string nomeFicheiro = argv[5];
    double angle1,angle3 = 0;
    double angle2 = (2 * pi) / slices;
    double angle4 = pi / stacks;
    double p1x, p1y, p1z;
    double p2x, p2y, p2z;
    double p3x, p3y, p3z;
    double p4x, p4y, p4z;
    std::cout << "File:" << nomeFicheiro << "\n";
    if (open(nomeFicheiro)) {
       
        for (int i = 0; i < stacks; i++) {
            angle1 = 0;
            for (int j = 0; j < slices; j++) {
                p1x = radius* sin(angle3)* sin(angle1);
                p1y = radius* cos(angle3);
                p1z = radius* sin(angle3)* cos(angle1);
                p2x = radius* sin(angle3 + angle4) * sin(angle1);
                p2y = radius* cos(angle3 + angle4);
                p2z = radius* sin(angle3 + angle4) * cos(angle1);
                p3x = radius* sin(angle3)* sin(angle1 + angle2);
                p3y = radius* cos(angle3);
                p3z = radius* sin(angle3)* cos(angle1 + angle2);
                p4x = radius* sin(angle3 + angle4) * sin(angle1 + angle2);
                p4y = radius* cos(angle3 + angle4);
                p4z = radius* sin(angle3 + angle4) * cos(angle1 + angle2);

                outFile << p1x << " " << p1y << " " << p1z << "\n";
                outFile << p2x << " " << p2y << " " << p2z << "\n";
                outFile << p4x << " " << p4y << " " << p4z << "\n";
                outFile << p1x << " " << p1y << " " << p1z << "\n";
                outFile << p4x << " " << p4y << " " << p4z << "\n";
                outFile << p3x << " " << p3y << " " << p3z << "\n";
                
                angle1 += angle2;
            }
            
            angle3 += angle4;
        }
    }

}

void plane(char* argv[]) {
    std::string slength = (argv[2]);
    double length = (double)atof(slength.c_str());
    double halfdouble = length / 2;
    std::string half = std::to_string(halfdouble);
    //Face voltada para cima 
    //Triangulo 1 => coordenadas do tri.esq.
    std::string triangulo1_1;
    triangulo1_1 = "-" + half + " 0 -" + half;

    std::string triangulo1_2;
    triangulo1_2 = "-" + half + " 0 " + half;

    std::string triangulo1_3;
    triangulo1_3 = half + " 0 " + half;


    // Triângulo 2 => Coordenadas do triângulo à direita.
    std::string triangulo2_1;//verificar
    triangulo2_1 = "-" + half + " 0 -" + half;

    std::string triangulo2_2;
    triangulo2_2 = half + " 0 " + half;

    std::string triangulo2_3;
    triangulo2_3 = half + " 0 -" + half;


    //Face voltada para baixo
    //Esq
    std::string triangulo1_3B;
    triangulo1_3B = half + " 0 " + half;

    std::string triangulo1_2B;
    triangulo1_2B = "-" + half + " 0 " + half;

    std::string triangulo1_1B;
    triangulo1_1B = "-" + half + " 0 -" + half;

    //Dir
    std::string triangulo2_1B;
    triangulo2_1B = "-" + half + " 0 -" + half;

    std::string triangulo2_3B;
    triangulo2_3B = half + " 0 -" + half;

    std::string triangulo2_2B;
    triangulo2_2B = half + " 0 " + half;


    std::string nomeFicheiro = argv[3];
    std::cout << "File:" << nomeFicheiro << "\n";
    if (open(nomeFicheiro)) {
       
        writeLine(triangulo1_1);
        writeLine(triangulo1_2);
        writeLine(triangulo1_3);
        writeLine(triangulo2_1);
        writeLine(triangulo2_2);
        writeLine(triangulo2_3);
        writeLine(triangulo1_1B);
        writeLine(triangulo1_3B);
        writeLine(triangulo1_2B);
        writeLine(triangulo2_1B);
        writeLine(triangulo2_3B);
        writeLine(triangulo2_2B);


        close();
    }
    else {
        std::cout << "Error opening file.\n";
    }

}

void box(int argc, char* argv[]) {
    std::string nomeFicheiro;
    if (argc == 6) { nomeFicheiro = argv[5]; }
    else { nomeFicheiro = argv[6]; }
    std::cout << "File:" << nomeFicheiro << "\n";
    if (open(nomeFicheiro)) {

        double length = (double)atof(argv[2]);
        double height = (double)atof(argv[3]);
        double width = (double)atof(argv[4]);
        int div = 1;
        if (argc == 7) div = (double)atof(argv[5]);

        double x, y, z;
        double divX = length / div;
        double divY = height / div;
        double divZ = width / div;

        //Face frontal
        x = length / 2;
        z = width / 2;
        std::string frontal;
        for (int i = 0; i < div; i++) {
            y = height;
            for (int j = 0; j < div; j++) {

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << (z - divZ) << "\n";

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << (z - divZ) << "\n";
                outFile << "" << x << " " << y << " " << (z - divZ) << "\n";

                y -= divY;
            }
            z -= divZ;
        }

        //Face direita
        x = -length / 2;
        z = -width / 2;
        std::string direita;
        for (int i = 0; i < div; i++) {
            y = height;
            for (int j = 0; j < div; j++) {

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << (x + divX) << " " << y << " " << z << "\n";
                outFile << "" << (x + divX) << " " << (y - divY) << " " << z << "\n";

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << (x + divX) << " " << (y - divY) << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << z << "\n";

                y -= divY;
            }
            x += divX;
        }

        //Face esquerda
        x = -length / 2;
        z = width / 2;
        std::string esquerda;
        for (int i = 0; i < div; i++) {
            y = height;
            for (int j = 0; j < div; j++) {

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << z << "\n";
                outFile << "" << (x + divX) << " " << y << " " << z << "\n";

                outFile << "" << x << " " << (y - divY) << " " << z << "\n";
                outFile << "" << (x + divX) << " " << (y - divY) << " " << z << "\n";
                outFile << "" << (x + divX) << " " << y << " " << z << "\n";

                y -= divY;
            }
            x += divX;
        }

        //Tampa superior
        y = height;
        z = -width / 2;
        std::string superior;
        for (int i = 0; i < div; i++) {
            x = -length / 2;
            for (int j = 0; j < div; j++) {

                outFile << "" << (x + divX) << " " << y << " " << z << "\n";
                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << x << " " << y << " " << (z + divZ) << "\n";

                outFile << "" << (x + divX) << " " << y << " " << z << "\n";
                outFile << "" << x << " " << y << " " << (z + divZ) << "\n";
                outFile << "" << (x + divX) << " " << y << " " << (z + divZ) << "\n";
                x += divX;
            }
            z += divZ;
        }


        //Base
        z = width / 2;
        std::string base;
        for (int i = 0; i < div; i++) {
            x = (-length / 2);
            for (int j = 0; j < div; j++) {

                outFile << "" << x << " 0 " << (z - divZ) << "\n";
                outFile << "" << (x + divX) << " 0 " << z << "\n";
                outFile << "" << x << " 0 " << z << "\n";

                outFile << "" << x << " 0 " << (z - divZ) << "\n";
                outFile << "" << (x + divX) << " 0 " << (z - divZ) << "\n";
                outFile << "" << (x + divX) << " 0 " << z << "\n";

                x += divX;
            }
            z -= divZ;
        }

        //Face traseira
        x = -length / 2;
        z = -width / 2;
        std::string traseira;
        for (int i = 0; i < div; i++) {
            y = height;
            for (int j = 0; j < div; j++) {

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << (z + divZ) << "\n";

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << (z + divZ) << "\n";
                outFile << "" << x << " " << y << " " << (z + divZ) << "\n";


                y -= divY;
            }
            z += divZ;
        }
    }
}
void cone(char* argv[])
{
    double pi = 3.14159265358979323846;
    double radius = atof(argv[2]);
    double height = atof(argv[3]);
    double slices = atof(argv[4]);
    double stacks = atof(argv[5]);
    double angle1 = 0;
    double angle2 = (2 * pi) / slices;
    std::string nomeFicheiro = argv[6];
    std::cout << "File:" << nomeFicheiro << "\n";
     if (open(nomeFicheiro)) {
       
        for (int i = 0; i < slices; i++) {
            outFile << "" << "0" << " 0 " << "0" << "\n";
            outFile << "" << (radius * sin(angle1 + angle2)) << " 0 " << (radius * cos(angle1 + angle2)) << "\n";
            outFile << "" << (radius * sin(angle1)) << " 0 " << (radius * cos(angle1)) << "\n";
            outFile << "" << "0 " << height << " 0" << "\n";
            outFile << "" << (radius * sin(angle1)) << " 0 " << (radius * cos(angle1)) << "\n";
            outFile << "" << (radius * sin(angle1 + angle2)) << " 0 " << (radius * cos(angle1 + angle2)) << "\n";
            angle1 += angle2;
        }
    }

}

int main(int argc, char* argv[]) {
    if (argc > 0) {
        std::string type = argv[1];
        if (type == "sphere") {
            if (argc == 6) {
                
                //open(file);
                esfera(argv);
                close();
                std::cout << "Success!\n";
            }
            else {
                std::cout << "Not enough arguments.\n";
            }

        }
        if (type == "cone") {
            if (argc == 7) {
                cone(argv);
                close();
                std::cout << "Success!\n";
                }
            else {
                    std::cout << "Not enough arguments.\n";
            }
            
        }
        if (type == "plane") {
            if (argc == 4) {
                plane(argv);
                std::cout << "Success!\n";
            }
            else {
                std::cout << "Not enough arguments.\n";
            }
        }
        if (type == "box") {
            if (argc == 6 || argc == 7) {
                box(argc, argv);
                close();
                std::cout << "Success!\n";
            }
            else {
                std::cout << "Not enough arguments.\n";
            }
        }

    }
    


    return 0;

}