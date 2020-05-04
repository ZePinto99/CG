#include <stdlib.h>
#include <stdio.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <string>
#include "parser.h"


std::vector<std::string> Trabalho;
std::vector<Ponto> triangles;

vector<OperFile*> files; // Vector de OperFiles (que relacionam os ficheiros 
					     // com as suas respetivas transformações).

int* fiVertexCount;
GLuint vertexCount;
GLuint buffers[20];
double** vertexB;

////////////////////////////////Curvas///////////////////////////////////////
#define POINT_COUNT 5
// Points that make up the loop for catmull-rom interpolation
float p[POINT_COUNT][3] = { {-1,-1,0},{-1,1,0},{1,1,0},{0,0,0},{1,-1,0} };

void buildRotMatrix(float* x, float* y, float* z, float* m) {

	m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
	m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
	m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
	m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}


void cross(float* a, float* b, float* res) {

	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}


void normalize(float* a) {

	float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0] / l;
	a[1] = a[1] / l;
	a[2] = a[2] / l;
}


float length(float* v) {

	float res = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	return res;

}

void multMatrixVector(float* m, float* v, float* res) {

	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}

}

void multMatrix14(double uW[1][4], double m[4][4], double res[1][4])
{
	res[0][0] = uW[0][0] * m[0][0] + uW[0][1] * m[1][0] + uW[0][2] * m[2][0] + uW[0][3] * m[3][0];
	res[0][1] = uW[0][0] * m[0][1] + uW[0][1] * m[1][1] + uW[0][2] * m[2][1] + uW[0][3] * m[3][1];
	res[0][2] = uW[0][0] * m[0][2] + uW[0][1] * m[1][2] + uW[0][2] * m[2][2] + uW[0][3] * m[3][2];
	res[0][3] = uW[0][0] * m[0][3] + uW[0][1] * m[1][3] + uW[0][2] * m[2][3] + uW[0][3] * m[3][3];
}

void multMatrix1441(double trans[1][4], double vH[4][1], double* res)
{
	*res = trans[0][0] * vH[0][0] + trans[0][1] * vH[1][0] + trans[0][2] * vH[2][0] + trans[0][3] * vH[3][0];
}


void getCatmullRomPoint(double t, double* p0, double* p1, double* p2, double* p3, double* pos, float* deriv) {

	double x, y, z;
	double res[1][4];
	double time[1][4] = { { powf(t,3), powf(t,2), t, 1 } };

	double m[4][4] = { {-0.5f,  1.5f, -1.5f,  0.5f},
					   { 1.0f, -2.5f,  2.0f, -0.5f},
					   {-0.5f,  0.0f,  0.5f,  0.0f},
					   { 0.0f,  1.0f,  0.0f,  0.0f} };

	double px[4][1] = { { p0[0] }, { p1[0] }, { p2[0] }, { p3[0] } };
	double py[4][1] = { { p0[1] }, { p1[1] }, { p2[1] }, { p3[1] } };
	double pz[4][1] = { { p0[2] }, { p1[2] }, { p2[2] }, { p3[2] } };

	multMatrix14(time, m, res);
	multMatrix1441(res, px, &x);
	multMatrix1441(res, py, &y);
	multMatrix1441(res, pz, &z);

	pos[0] = x;
	pos[1] = y;
	pos[2] = z;

	float a[4][4] = { 0 };
	float pxx[4];
	float pyy[4];
	float pzz[4];
	float dT[4] = { 3 * powf(t,2), static_cast<float>(2 * t), 1, 0 };

	pxx[0] = px[0][0]; pxx[1] = px[1][0]; pxx[2] = px[2][0]; pxx[3] = px[3][0];
	pyy[0] = py[0][0]; pyy[1] = py[1][0]; pyy[2] = py[2][0]; pyy[3] = py[3][0];
	pzz[0] = pz[0][0]; pzz[1] = pz[1][0]; pzz[2] = pz[2][0]; pzz[3] = pz[3][0];

	multMatrixVector(reinterpret_cast<float*>(*m), pxx, a[0]);
	multMatrixVector(reinterpret_cast<float*>(*m), pyy, a[1]);
	multMatrixVector(reinterpret_cast<float*>(*m), pzz, a[2]);
	multMatrixVector(*a, dT, deriv);
}


// given  global t, returns the point in the curve
void getGlobalCatmullRomPoint(double gt, double* pos, float* deriv, double** curve, int size) {

	double t = gt * size; 
	int index = floor(t);
	t = t - index;

	int indices[4];
	indices[0] = (index + size - 1) % size;
	indices[1] = (indices[0] + 1) % size;
	indices[2] = (indices[1] + 1) % size;
	indices[3] = (indices[2] + 1) % size;

	getCatmullRomPoint(t, curve[indices[0]], curve[indices[1]], curve[indices[2]], curve[indices[3]], pos, deriv);
}

////////////////////////////////////////////////////////////////

void createVBO(int i) {
	int p = 0;  int vertex = 0;

	vector<Ponto>::iterator tri = triangles.begin();
	vector<OperFile*>::iterator itFile;

	vertexB[i] = (double*)malloc(sizeof(double) * triangles.size() * 3);

	while (tri != triangles.end()) {

		Ponto aux_1 = *tri; tri++;
		Ponto aux_2 = *tri; tri++;
		Ponto aux_3 = *tri; tri++;

		vertexB[i][p] = aux_1.x;
		vertexB[i][p + 1] = aux_1.y;
		vertexB[i][p + 2] = aux_1.z;
		vertex++;

		vertexB[i][p + 3] = aux_2.x;
		vertexB[i][p + 4] = aux_2.y;
		vertexB[i][p + 5] = aux_2.z;
		vertex++;

		vertexB[i][p + 6] = aux_3.x;
		vertexB[i][p + 7] = aux_3.y;
		vertexB[i][p + 8] = aux_3.z;
		vertex++;

		p += 9;
	}
	vertexCount = vertex;
}


void lerficheiro(std::string nomeficheiro)
{
	int c = 0;

	double storedouble[3];
	for (int i = 0; i < 3; i++) storedouble[i] = 0;

	std::ifstream trigsFile;
	trigsFile.open(nomeficheiro);

	std::string linha;
	if (trigsFile.is_open()) {
		while (getline(trigsFile, linha)) {
			std::string sTmp;
			for (int i = 0; i <= linha.length(); i++) {
				if (linha[i] == ' ' || linha[i] == '\0') {
					storedouble[c] = stod(sTmp);
					sTmp.clear();
					c++;
				}
				else {
					sTmp.push_back(linha[i]);
				}
			}
			Ponto aux;
			aux.x = storedouble[0];
			aux.y = storedouble[1];
			aux.z = storedouble[2];
			triangles.push_back(aux);
			c = 0;
		}
	}
}

void lertudoemaisalgumacoisa() {
	vertexB = (double**)malloc(sizeof(double*) * files.size());
	std::vector<OperFile*>::iterator itout;
	itout = files.begin();
	int i=0;
	while (itout != files.end()) {
		triangles.clear();
		OperFile* op = *itout;
		OperFile* opcopy = op;
		char* stds = op->fileName;
		std::string std(stds);
		lerficheiro(std);
		createVBO(i);
		itout++;
		i++;
	}
}


void dynamicTranslate(Oper* oper, int i)
{
	Transform* t = oper->transform;
	static double intervalos[20];
	static double checkpoints[20];
	static double time = 0;

	//criar matriz aqui
	int h = (t->points).size();
	int cols = 3;
	double** curve = new double* [h];
	
	for (int i = 0; i < h; ++i)
		curve[i] = new double[cols];


	int size = 0;
	vector<Ponto>::iterator it;
	for (it = t->points.begin(); it != t->points.end(); it++) {
		Ponto aux = *it;

		curve[size][0] = aux.x;
		curve[size][1] = aux.y;
		curve[size][2] = aux.z;

		size++;
	}

	double pos[4];
	float deriv[4];

	glBegin(GL_LINE_LOOP);
		int n = 100;
		for (int j = 1; j < n; j++) {
			getGlobalCatmullRomPoint((double) j / n, pos, deriv, curve, size);
			glVertex3d(pos[0], pos[1], pos[2]);
		}
	glEnd();

	getGlobalCatmullRomPoint(intervalos[i], pos, deriv, curve, size);
	
	if (time < files.size()) {
		double ttt = glutGet(GLUT_ELAPSED_TIME);
		checkpoints[i] = 1 / (t->time * 1000);
		time++;
	}
	intervalos[i] += checkpoints[i];

	//apagar matriz aqui
	for (int i = 0; i < h; ++i)
		delete[] curve[i];
	delete[] curve;
	

	glTranslated(pos[0], pos[1], pos[2]);
}


void dynamicRotate(Oper* oper, int i)
{
	Transform* t = oper->transform;

	double angles[20];
	double time = t->time;

	angles[i] += 360 / (time * 1000);
	glRotated(angles[i], oper->x, oper->y, oper->y);
}

void desenhar(void)
{
	std::vector<OperFile*>::iterator itout;
	itout = files.begin();
	int i = 0;
	
	while (itout != files.end()) {
		triangles.clear();
		OperFile* op = *itout;
		itout++;
		
		char* stds = op->fileName;
		glPushMatrix();

			std::vector<Oper*>::iterator it2;
			it2 = op->operations.begin();
			while(it2 != op->operations.end())
			{	
				Oper* oper = *it2;
				it2++;
				if (strcmp(oper->operation, "translate")==0)
				{
					if (oper->transform == nullptr) {
						glTranslatef(oper->x, oper->y, oper->z);
					}
					else {
						dynamicTranslate(oper, i);
					}
				}
				if (strcmp(oper->operation, "rotate") == 0)
				{
					if (oper->transform == nullptr) {
						glRotatef(oper->angle, oper->x, oper->y, oper->z);
					}
					else {
						dynamicRotate(oper, i);
					}
				}
				if (strcmp(oper->operation, "scale") == 0) {
					glScalef(oper->x, oper->y, oper->z) ;
				}
			}
			glColor3f(1, 1, 1);
			glBindBuffer(GL_ARRAY_BUFFER, buffers[i]);
			glBufferData(GL_ARRAY_BUFFER, vertexCount * sizeof(double) * 3, vertexB[i], GL_STATIC_DRAW);
			glVertexPointer(3, GL_DOUBLE, 0, 0);
			glDrawArrays(GL_TRIANGLES, 0, vertexCount);
			glPopMatrix();
			i++;
	}

}

void changeSize(int w, int h)
{
	// Prevent a divide by zero, when window is too short
	// (you can’t make a window with zero width).
	if (h == 0)
		h = 1;
	// compute window's aspect ratio
	float ratio = w * 1.0f / h;
	// Set the projection matrix as current
	glMatrixMode(GL_PROJECTION);
	// Load the identity matrix
	glLoadIdentity();
	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);
	// Set the perspective
	gluPerspective(45.0f, ratio, 10.0f, 10000.0f);
	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}

float lx = 30.0, ly = 30.0, lz = 30.0;

void renderScene(void)
{
	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set the camera
	glLoadIdentity();


	gluLookAt(lx, ly, lz,
		20.0, 5.0, 0.0,
		0.0f, 1.0f, 0.0f);

	// put drawing instructions here
	glEnableClientState(GL_VERTEX_ARRAY);
	int size = files.size();
	glGenBuffers(size, buffers);
	desenhar();
	

	// End of frame
	glutSwapBuffers();
}

void processSpecialKeys(int key, int xx, int yy)
{
	switch (key) {
	case GLUT_KEY_RIGHT:
		glTranslatef(lx += 1, ly += 0, lz += 0);
		glutPostRedisplay();
		break;
	case GLUT_KEY_LEFT:
		glTranslatef(lx -= 1, ly -= 0, lz -= 0);
		glutPostRedisplay();
		break;
	case GLUT_KEY_UP:
		glTranslatef(lx -= 1, ly -= 1, lz -= 1);
		glutPostRedisplay();
		break;
	case GLUT_KEY_DOWN:
		glTranslatef(lx += 1, ly += 1, lz += 1);
		glutPostRedisplay();
		break;

	default:
		std::cout << "Não conheço esse comando!" << "\n";
	}
}



int main(int argc, char** argv)
{
	// put GLUT’s init here
	xmlParser("sistemaSolarDinamico.xml", files);

	lertudoemaisalgumacoisa();

	//if (!files.empty()) std::cout << "ola";
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(1000, 1000);
	glutCreateWindow("CG_Trabalho_prático");
	// put callback registry here
	glutReshapeFunc(changeSize);
	glutSpecialFunc(processSpecialKeys);
	glutIdleFunc(renderScene);
	glutDisplayFunc(renderScene);
	// some OpenGL settings

#ifndef __APPLE__
	glewInit();
#endif

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	// enter GLUT’s main cycle
	glutMainLoop();
	return 1;
}