#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <string>       // std::string
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <iostream>     // std::cout
#include <fstream>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include "generator.h"
//#define M_PI 3.14159265358979323846
std::ofstream outFile;

class Point {
public: float x, y, z;
};
int patches = 0;
int vertices = 0;
std::vector<std::vector<int> > indicesPatch;
std::vector<Point> points;

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
    double s_factor, t_factor;
    int radius = atoi(argv[2]);
    int slices = atoi(argv[3]);
    int stacks = atoi(argv[4]);
    std::string nomeFicheiro = argv[5];
    double angle1, angle3 = 0;
    double angle2 = (2 * pi) / slices;
    double angle4 = pi / stacks;
    double p1x, p1y, p1z;
    double p2x, p2y, p2z;
    double p3x, p3y, p3z;
    double p4x, p4y, p4z;
    std::cout << "File:" << nomeFicheiro << "\n";
    if (open(nomeFicheiro)) {
        for (double i = 0; i < stacks; i++) {
            t_factor = (stacks - i) * (1.d / stacks); 
            angle1 = 0;
            for (int j = 0; j < slices; j++) {
                s_factor = j * (1.d / slices);
                p1x = radius * sin(angle3) * sin(angle1);
                p1y = radius * cos(angle3);
                p1z = radius * sin(angle3) * cos(angle1);
                p2x = radius * sin(angle3 + angle4) * sin(angle1);
                p2y = radius * cos(angle3 + angle4);
                p2z = radius * sin(angle3 + angle4) * cos(angle1);
                p3x = radius * sin(angle3) * sin(angle1 + angle2);
                p3y = radius * cos(angle3);
                p3z = radius * sin(angle3) * cos(angle1 + angle2);
                p4x = radius * sin(angle3 + angle4) * sin(angle1 + angle2);
                p4y = radius * cos(angle3 + angle4);
                p4z = radius * sin(angle3 + angle4) * cos(angle1 + angle2);

                /////////////////////////////
                outFile << p1x << " " << p1y << " " << p1z << "\n";
                outFile << p2x << " " << p2y << " " << p2z << "\n";
                outFile << p4x << " " << p4y << " " << p4z << "\n";
                outFile << "" << radius * sin(angle3)*sin(angle1) << " " << radius * cos(angle3) << " " << radius * cos(angle1) * sin(angle3) << "\n";
                outFile << "" << radius * sin(angle3 + angle4) * sin(angle1) << " " << radius * cos(angle3 + angle4) << " " << radius * sin(angle3+angle4) * cos(angle1) << "\n";
                outFile << "" << radius * sin(angle3 + angle4) * sin(angle1+angle2) << " " << radius * cos(angle3 + angle4) << " " << radius * sin(angle3 + angle4) * cos(angle1+angle2) << "\n";
                outFile << s_factor << " " << t_factor << " 0\n";
                outFile << s_factor << " " << t_factor - (1.d/stacks) << " 0\n";
                outFile << s_factor + (1.d/slices) << " " << t_factor-(1.d/stacks) << " 0\n";


                outFile << p1x << " " << p1y << " " << p1z << "\n";
                outFile << p4x << " " << p4y << " " << p4z << "\n";
                outFile << p3x << " " << p3y << " " << p3z << "\n";
                outFile << "" << radius * sin(angle3) * sin(angle1) << " " << radius * cos(angle3) << " " << radius * cos(angle1) * sin(angle3) << "\n";
                outFile << "" << radius * sin(angle3 + angle4) * sin(angle1) << " " << radius * cos(angle3 + angle4) << " " << radius * sin(angle3 + angle4) * cos(angle1) << "\n";
                outFile << "" << radius * sin(angle3) * sin(angle1 + angle2) << " " << radius * cos(angle3) << " " << radius * sin(angle3) * cos(angle1 + angle2) << "\n";
                outFile << s_factor << " " << t_factor << " 0\n";
                outFile << s_factor + (1.d / slices) << " " << t_factor - (1.d / stacks) << " 0\n";
                outFile << s_factor + (1.d/slices) << " " << t_factor << " 0\n";




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



    std::string triangulo1_1Normal;
    triangulo1_1Normal = "0 1 0";


    std::string triangulo2_1Normal;
    triangulo2_1Normal = "0 -1 0";

    ////////////texturas

    std::string triangulo1_1Text;
    triangulo1_1Text = "0 1 0";

    std::string triangulo1_2Text;
    triangulo1_2Text = "1 1 0";

    std::string triangulo1_3Text;
    triangulo1_3Text = "0 1 0";

    /////////////////////////

    std::string triangulo2_1Text;
    triangulo2_1Text = "0 0 0";

    std::string triangulo2_2Text;
    triangulo2_2Text = "1 0 0";

    std::string triangulo2_3Text;
    triangulo2_3Text = "0 1 0";



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
        writeLine(triangulo1_1Normal);
        writeLine(triangulo1_1Normal);
        writeLine(triangulo1_1Normal);
        writeLine(triangulo1_1Text);
        writeLine(triangulo1_2Text);
        writeLine(triangulo1_3Text);


        writeLine(triangulo2_1);
        writeLine(triangulo2_2);
        writeLine(triangulo2_3);
        writeLine(triangulo1_1Normal);
        writeLine(triangulo1_1Normal);
        writeLine(triangulo1_1Normal);
        writeLine(triangulo2_1Text);
        writeLine(triangulo2_1Text);
        writeLine(triangulo2_1Text);


        writeLine(triangulo1_1B);
        writeLine(triangulo1_3B);
        writeLine(triangulo1_2B);
        writeLine(triangulo2_1Normal);
        writeLine(triangulo2_1Normal);
        writeLine(triangulo2_1Normal);
        writeLine(triangulo1_1Text);
        writeLine(triangulo1_2Text);
        writeLine(triangulo1_3Text);


        writeLine(triangulo2_1B);
        writeLine(triangulo2_3B);
        writeLine(triangulo2_2B);
        writeLine(triangulo2_1Normal);
        writeLine(triangulo2_1Normal);
        writeLine(triangulo2_1Normal);
        writeLine(triangulo2_1Text);
        writeLine(triangulo2_2Text);
        writeLine(triangulo2_3Text);


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



        double textura1 = 0.25 / div;
        double textura2 = 0.33 / div;
        //   1 __ __ __ __
        //    |  |tp|  |  |
        //0.66|__|__|__|__|
        //    |Lf|Fr|Ri|Bk|
        //0.33|__|__|__|__|
        //    |  |Ba|  |  |
        //   0|__|__|__|__|
        //   0 .25 .50 .75 1
        //
        //
        //

        double frenteX = 0.25;
        double frenteY = 0.33;
        double trasX = 0.75;
        double trasY = 0.33;
        float esquerdaX = 0;
        float esquerdaY = 0.33;
        float direitaX = 0.50;
        float direitaY = 0.33;
        float cimaX = 0.25;
        float cimaY = 0.66;
        float baixoX = 0.25;
        float baixoY = 0;
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
                outFile << "0 0 1\n";
                outFile << "0 0 1\n";
                outFile << "0 0 1\n";
                outFile << frenteX << " " << frenteY << " 0\n";
                outFile << frenteX + textura1 << " " << frenteY << " 0\n";
                outFile << frenteX + textura1 << " " << frenteY + textura2 << " 0\n";

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << (z - divZ) << "\n";
                outFile << "" << x << " " << y << " " << (z - divZ) << "\n";
                outFile << "0 0 1\n";
                outFile << "0 0 1\n";
                outFile << "0 0 1\n";
                outFile << frenteX + textura1 << " " << frenteY + textura2 << " 0\n";
                outFile << frenteX << " " << frenteY + textura2 << " 0\n";
                outFile << frenteX << " " << frenteY  << " 0\n";
                frenteX += textura1;


                y -= divY;
            }
            frenteY += textura2;
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
                outFile << "1 0 0\n";
                outFile << "1 0 0\n";
                outFile << "1 0 0\n";
                outFile << direitaX << " " << direitaY  << " 0\n";
                outFile << direitaX + textura1 << " " << direitaY << " 0\n";
                outFile << direitaX + textura1 << " " << direitaY + textura2 << " 0\n";

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << (x + divX) << " " << (y - divY) << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << z << "\n";
                outFile << "1 0 0\n";
                outFile << "1 0 0\n";
                outFile << "1 0 0\n";
                outFile << direitaX + textura1 << " " << direitaY + textura2 << " 0\n";
                outFile << direitaX << " " << direitaY + textura2<< " 0\n";
                outFile << direitaX << " " << direitaY  << " 0\n";

                direitaX += textura1;
                y -= divY;
            }
            direitaY += textura2;
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
                outFile << "-1 0 0\n";
                outFile << "-1 0 0\n";
                outFile << "-1 0 0\n";
                outFile << esquerdaX << " " << esquerdaY << " 0\n";
                outFile << esquerdaX + textura1 << " " << esquerdaY << " 0\n";
                outFile << esquerdaX + textura1 << " " << esquerdaY + textura2 << " 0\n";

                outFile << "" << x << " " << (y - divY) << " " << z << "\n";
                outFile << "" << (x + divX) << " " << (y - divY) << " " << z << "\n";
                outFile << "" << (x + divX) << " " << y << " " << z << "\n";
                outFile << "-1 0 0\n";
                outFile << "-1 0 0\n";
                outFile << "-1 0 0\n";
                outFile << esquerdaX + textura1 << " " << esquerdaY + textura2 << " 0\n";
                outFile << esquerdaX << " " << esquerdaY + textura2 << " 0\n";
                outFile << esquerdaX << " " << esquerdaY << " 0\n";


                esquerdaX += textura1;

                y -= divY;
            }
            esquerdaY += textura2;
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
                outFile << "0 1 0\n";
                outFile << "0 1 0\n";
                outFile << "0 1 0\n";
                outFile << cimaX << " " << cimaY << " 0\n";
                outFile << cimaX + textura1 << " " << cimaY << " 0\n";
                outFile << cimaX + textura1 << " " << cimaY + textura2 << " 0\n";


                outFile << "" << (x + divX) << " " << y << " " << z << "\n";
                outFile << "" << x << " " << y << " " << (z + divZ) << "\n";
                outFile << "" << (x + divX) << " " << y << " " << (z + divZ) << "\n";
                outFile << "0 1 0\n";
                outFile << "0 1 0\n";
                outFile << "0 1 0\n";
                outFile << cimaX + textura1 << " " << cimaY + textura2 << " 0\n";
                outFile << cimaX << " " << cimaY + textura2 << " 0\n";
                outFile << cimaX << " " << cimaY << " 0\n";



                cimaX += textura1;
                x += divX;
            }
            cimaY += textura2;
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
                outFile << "0 -1 0\n";
                outFile << "0 -1 0\n";
                outFile << "0 -1 0\n";
                outFile << baixoX << " " << baixoY << " 0\n";
                outFile << baixoX + textura1 << " " << baixoY << " 0\n";
                outFile << baixoX + textura1 << " " << baixoY + textura2 << " 0\n";

                outFile << "" << x << " 0 " << (z - divZ) << "\n";
                outFile << "" << (x + divX) << " 0 " << (z - divZ) << "\n";
                outFile << "" << (x + divX) << " 0 " << z << "\n";
                outFile << "0 -1 0\n";
                outFile << "0 -1 0\n";
                outFile << "0 -1 0\n";
                outFile << baixoX + textura1 << " " << baixoY + textura2 << " 0\n";
                outFile << baixoX << " " << baixoY + textura2 << " 0\n";
                outFile << baixoX << " " << baixoY << " 0\n";

                baixoX += textura1;
                x += divX;
            }
            baixoY += textura2;
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
                outFile << "0 0 -1\n";
                outFile << "0 0 -1\n";
                outFile << "0 0 -1\n";
                outFile << trasX << " " << trasY << " 0\n";
                outFile << trasX + textura1 << " " << trasY << " 0\n";
                outFile << trasX + textura1 << " " << trasY + textura2 << " 0\n";

                outFile << "" << x << " " << y << " " << z << "\n";
                outFile << "" << x << " " << (y - divY) << " " << (z + divZ) << "\n";
                outFile << "" << x << " " << y << " " << (z + divZ) << "\n";
                outFile << "0 0 -1\n";
                outFile << "0 0 -1\n";
                outFile << "0 0 -1\n";
                outFile << trasX + textura1 << " " << trasY + textura2 << " 0\n";
                outFile << trasX << " " << trasY + textura2 << " 0\n";
                outFile << trasX << " " << trasY << " 0\n";

                trasX += textura1;
                y -= divY;
            }
            trasY += textura2;
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
            outFile << "0 -1 0\n";
            outFile << "0 -1 0\n";
            outFile << "0 -1 0\n";
            ////texturas
            outFile << "0 -1 0\n";
            outFile << "0 -1 0\n";
            outFile << "0 -1 0\n";







            outFile << "" << "0 " << height << " 0" << "\n";
            outFile << "" << (radius * sin(angle1)) << " 0 " << (radius * cos(angle1)) << "\n";
            outFile << "" << (radius * sin(angle1 + angle2)) << " 0 " << (radius * cos(angle1 + angle2)) << "\n";
            outFile << "0 1 0\n";
            outFile << "" << radius * sin(i + 1) * angle1 << " 0 " << radius * cos(i + 1) * (angle1) << "\n";
            outFile << "" << radius * sin(i + 1) * (angle1+angle2) << " 0 " << radius * cos(i + 1) * (angle1 + angle2) << "\n";
            ////texturas
            outFile << "0 -1 0\n";
            outFile << "0 -1 0\n";
            outFile << "0 -1 0\n";



            angle1 += angle2;
        }
    }

}

void writeBezier(std::string newfile, int tesslevel) {
    ////////////////////////////////////////
   ////// Escrever Patches no ficheiro ////
   ////////////////////////////////////////
    float res[3];
    float t;
    int index, indices[4];
    float q[4][tesslevel][3];

    float tess = 1 / ((float)tesslevel - 1);
    float r[tesslevel][tesslevel][3];
    float pontos = patches * (tesslevel) * 2 * (tesslevel) * 3 * 3;


    if (open(newfile)) {
        /////////////////////////////
        ///escrever num de pontos////
        /////////////////////////////

        //outFile << pontos; outFile << '\n';
        for (int n = 0; n < patches; n++) {
            /////////////////////////////////////////////////////
            /////recolher os vértices de points para x y z //////
            /////////////////////////////////////////////////////
            float p[16][3];
            for (int m = 0; m < 16; m++) {
                p[m][0] = points[indicesPatch[n][m]].x;
                p[m][1] = points[indicesPatch[n][m]].y;
                p[m][2] = points[indicesPatch[n][m]].z;
            }
            int j = 0, k = 0;
            /////////////////////////////////
            ///desinahr 4 curvas de Bezier///
            /////////////////////////////////
            for (int i = 0; i < 15; i += 4) {
                indices[0] = i;
                indices[1] = i + 1;
                indices[2] = i + 2;
                indices[3] = i + 3;
                /////////////////////////////////////////////////////////////////////
                ///////////////calcular a curva//////////////////////////////////////
                ///  [-t^3 + 3t^2 - 3t,3t^3-6t^2+3t, -3t^3+3t^2,t^3][P0] ////////////
                /////////////////////////////////////////////////// [P1] ////////////
                /////////////////////////////////////////////////// [P2] ////////////
                //////////////////////////////////////////////////  [P3] ////////////
                /////////////////////////////////////////////////////////////////////
                for (int a = 0; a < tesslevel - 1; a++) {
                    t = a * tess;
                    index = floor(t);
                    t = t - index;
                    res[0] = (-p[indices[0]][0] + 3 * p[indices[1]][0] - 3 * p[indices[2]][0] + p[indices[3]][0]) * t * t * t + (3 * p[indices[0]][0] - 6 * p[indices[1]][0] + 3 * p[indices[2]][0]) * t * t + (-3 * p[indices[0]][0] + 3 * p[indices[1]][0]) * t + p[indices[0]][0];
                    res[1] = (-p[indices[0]][1] + 3 * p[indices[1]][1] - 3 * p[indices[2]][1] + p[indices[3]][1]) * t * t * t + (3 * p[indices[0]][1] - 6 * p[indices[1]][1] + 3 * p[indices[2]][1]) * t * t + (-3 * p[indices[0]][1] + 3 * p[indices[1]][1]) * t + p[indices[0]][1];
                    res[2] = (-p[indices[0]][2] + 3 * p[indices[1]][2] - 3 * p[indices[2]][2] + p[indices[3]][2]) * t * t * t + (3 * p[indices[0]][2] - 6 * p[indices[1]][2] + 3 * p[indices[2]][2]) * t * t + (-3 * p[indices[0]][2] + 3 * p[indices[1]][2]) * t + p[indices[0]][2];
                    q[j][k][0] = res[0];
                    q[j][k][1] = res[1];
                    q[j][k][2] = res[2];
                    k++;
                }
                /////////////////////////////////////////////////////////
                //////////////////// caso tesslevel seja um /////////////
                ///////////////////  slavaguarda ////////////////////////
                /////////////////////////////////////////////////////////
                t = 1;

                res[0] = (-p[indices[0]][0] + 3 * p[indices[1]][0] - 3 * p[indices[2]][0] + p[indices[3]][0]) * t * t * t + (3 * p[indices[0]][0] - 6 * p[indices[1]][0] + 3 * p[indices[2]][0]) * t * t + (-3 * p[indices[0]][0] + 3 * p[indices[1]][0]) * t + p[indices[0]][0];
                res[1] = (-p[indices[0]][1] + 3 * p[indices[1]][1] - 3 * p[indices[2]][1] + p[indices[3]][1]) * t * t * t + (3 * p[indices[0]][1] - 6 * p[indices[1]][1] + 3 * p[indices[2]][1]) * t * t + (-3 * p[indices[0]][1] + 3 * p[indices[1]][1]) * t + p[indices[0]][1];
                res[2] = (-p[indices[0]][2] + 3 * p[indices[1]][2] - 3 * p[indices[2]][2] + p[indices[3]][2]) * t * t * t + (3 * p[indices[0]][2] - 6 * p[indices[1]][2] + 3 * p[indices[2]][2]) * t * t + (-3 * p[indices[0]][2] + 3 * p[indices[1]][2]) * t + p[indices[0]][2];

                q[j][k][0] = res[0];
                q[j][k][1] = res[1];
                q[j][k][2] = res[2];
                j++;
                k = 0;
            }

            for (int j = 0; j < tesslevel; j++) {
                for (int a = 0; a < tesslevel - 1; a++) {
                    t = a * tess;;
                    index = floor(t);
                    t = t - index;

                    res[0] = (-q[0][j][0] + 3 * q[1][j][0] - 3 * q[2][j][0] + q[3][j][0]) * t * t * t + (3 * q[0][j][0] - 6 * q[1][j][0] + 3 * q[2][j][0]) * t * t + (-3 * q[0][j][0] + 3 * q[1][j][0]) * t + q[0][j][0];
                    res[1] = (-q[0][j][1] + 3 * q[1][j][1] - 3 * q[2][j][1] + q[3][j][1]) * t * t * t + (3 * q[0][j][1] - 6 * q[1][j][1] + 3 * q[2][j][1]) * t * t + (-3 * q[0][j][1] + 3 * q[1][j][1]) * t + q[0][j][1];
                    res[2] = (-q[0][j][2] + 3 * q[1][j][2] - 3 * q[2][j][2] + q[3][j][2]) * t * t * t + (3 * q[0][j][2] - 6 * q[1][j][2] + 3 * q[2][j][2]) * t * t + (-3 * q[0][j][2] + 3 * q[1][j][2]) * t + q[0][j][2];
                    r[j][k][0] = res[0];
                    r[j][k][1] = res[1];
                    r[j][k][2] = res[2];
                    k++;
                }

                t = 1;

                res[0] = (-q[0][j][0] + 3 * q[1][j][0] - 3 * q[2][j][0] + q[3][j][0]) * t * t * t + (3 * q[0][j][0] - 6 * q[1][j][0] + 3 * q[2][j][0]) * t * t + (-3 * q[0][j][0] + 3 * q[1][j][0]) * t + q[0][j][0];
                res[1] = (-q[0][j][1] + 3 * q[1][j][1] - 3 * q[2][j][1] + q[3][j][1]) * t * t * t + (3 * q[0][j][1] - 6 * q[1][j][1] + 3 * q[2][j][1]) * t * t + (-3 * q[0][j][1] + 3 * q[1][j][1]) * t + q[0][j][1];
                res[2] = (-q[0][j][2] + 3 * q[1][j][2] - 3 * q[2][j][2] + q[3][j][2]) * t * t * t + (3 * q[0][j][2] - 6 * q[1][j][2] + 3 * q[2][j][2]) * t * t + (-3 * q[0][j][2] + 3 * q[1][j][2]) * t + q[0][j][2];

                r[j][k][0] = res[0];
                r[j][k][1] = res[1];
                r[j][k][2] = res[2];
                k = 0;
            }
            /////////////////////////////
            ////escrever no newfile//////
            /////////////////////////////
            ////nao usamos a wrtiteline//
            ////pois tinha de se pssar///
            ////os floats para string////
            ////o que colocava casas/////
            ////decimais a mais /////////
            /////////////////////////////
            for (int i = 0; i < tesslevel - 1; i++)
                for (int j = 0; j < tesslevel - 1; j++) {
                    outFile << r[i][j][0];   outFile << ' '; outFile << r[i][j][1];   outFile << ' '; outFile << r[i][j][2];   outFile << '\n';
                    outFile << r[i + 1][j][0]; outFile << ' '; outFile << r[i + 1][j][1]; outFile << ' '; outFile << r[i + 1][j][2];  outFile << '\n';
                    outFile << r[i][j + 1][0]; outFile << ' '; outFile << r[i][j + 1][1]; outFile << ' '; outFile << r[i][j + 1][2]; outFile << '\n';

                    outFile << r[i + 1][j][0];   outFile << ' '; outFile << r[i + 1][j][1];   outFile << ' '; outFile << r[i + 1][j][2];   outFile << '\n';
                    outFile << r[i + 1][j + 1][0]; outFile << ' '; outFile << r[i + 1][j + 1][1]; outFile << ' '; outFile << r[i + 1][j + 1][2]; outFile << '\n';
                    outFile << r[i][j + 1][0];   outFile << ' '; outFile << r[i][j + 1][1];   outFile << ' '; outFile << r[i][j + 1][2];   outFile << '\n';
                }
        }
    }
    else { printf("Error: Cannot open file...\n"); exit(0); }
}

void readBezier(std::string patchfilename) {
    std::fstream f;
    f.open(patchfilename, std::ios::in);
    if (f.is_open()) {
        /////////////////// le numero de patches
        std::string line;
        if (getline(f, line)) sscanf(line.c_str(), "%d\n", &patches);
        for (int i = 0; i < patches; i++) {
            std::vector<int> iPatch;
            std::getline(f, line);
            std::istringstream str1(line);
            std::string str2;
            /////////////////// patches
            while (std::getline(str1, str2, ',')) {
                /////////// adiciona indices dos patchs
                iPatch.push_back(atoi(str2.c_str()));

            }

            indicesPatch.push_back(iPatch);
        }
        //////////////////// cria os pontos e coloca no vetor points
        if (getline(f, line)) sscanf(line.c_str(), "%d\n", &vertices);
        for (int i = 0; i < vertices; i++) {
            float x = 0, y = 0, z = 0;
            if (getline(f, line)) sscanf(line.c_str(), "%f, %f, %f\n", &x, &y, &z);

            Point p;
            p.x = x; p.y = y; p.z = z;

            points.push_back(p);
        }

        f.close();
    }
    else { printf("Error: Not possible acess patch file..."); exit(0); }

}

void besier(const char* patchfilename, int tesslevel, const char* newfile) {

    //////////////////////////
    /////ler patch bezier/////
    //////////////////////////
    readBezier(patchfilename);
    //////////////////////////
    /////Desinhar bezier//////
    //////////////////////////
    writeBezier(newfile, tesslevel);

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
        if (type == "Bezier") {

            if (argc == 5) {
                std::ofstream file;
                std::string nomeFicheiro;
                nomeFicheiro = argv[4];
                int tess = atoi(argv[3]);
                besier(argv[2], tess, argv[4]);
                //writeLine( readbezier(argv));
                close();
                std::cout << "Success!\n";
            }
        }
    }



    return 0;

}