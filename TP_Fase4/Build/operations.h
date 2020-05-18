#ifndef CG_TP_OPERATIONS_H
#define CG_TP_OPERATIONS_H

#include <iostream>
#include <vector>
#include <string.h>

using namespace std;

/*
 Coordenadas de todos os triângulos.
*/
struct Ponto {
    double x;
    double y;
    double z;
};

/*
 Estrutura com informação relativa à iluminação.
 Tipo (POINT, DIRECTIONAL ou SPOT) e argumentos.
*/
typedef struct light {
    char* type;
    double x, y, z;
} Light;

/*
 Estrutura com informação relativa às cores.
 Component (diffuse, specular, emissive ou ambient) e RGB (diff,spec,emis,ambi).
*/
typedef struct color {
    char* component;
    double r, g, b;
} Color;

/*
 Estrutura com informação relativa a uma transformação dependente de um tempo.
 Tempo e vetor de pontos associados à transformação.
*/
typedef struct transform {
    double time;
    vector<Ponto> points;
} Transform;

/*
 Estrutura com informação realtiva a uma operação.
 Nome da operação e suas variáveis.
*/
typedef struct oper {
    char* operation;
    double x, y, z, angle;
    Transform* transform;
} Oper;


/*
 Estrutura que atribui a cada ficheiro .3d as respetivas
operações a realizar.
 */
typedef struct operFile {
    char* fileName;
    vector<Oper*> operations;
    char* texture;
    Color* color;
} OperFile;


/*
 Lista ligada auxiliar com informação realtiva a uma operação.
 Apoio à leitura de grupos secundários.
 */
typedef struct operationAux {
    char* operation;
    double x, y, z, angle;
    Transform* transform;
    struct operationAux* next;
} OperAux;


/*
 Array de hierarquias, que contém OperAux's, ou seja,
listas ligadas de operações e suas variáveis.
 */
typedef OperAux* Group[40];



void initGroup(Group g);

void addPoint(double x, double y, double z, vector<Ponto>& points);

void addOpGroup(Group g, int i, char opName[], double x, double y, double z, double angle, Transform* transform);

void addOperFile(char* fileName, Group group, vector<OperFile*>& files, char* texture, Color* color);

#endif //CG_TP_OPERATIONS_H