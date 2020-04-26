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
typedef OperAux* Group[10];



void initGroup(Group g);

void addPoint(double x, double y, double z, vector<Ponto>& points);

void addOpGroup(Group g, int i, char opName[], double x, double y, double z, double angle, Transform* transform);

void addOperFile(char* fileName, Group g, vector<OperFile*>& files);

#endif //CG_TP_OPERATIONS_H