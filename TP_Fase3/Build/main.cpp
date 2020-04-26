#include <stdlib.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
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

GLuint vertexCount;
GLuint buffers;
float* vertexB;

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


void getCatmullRomPoint(float t, float* p0, float* p1, float* p2, float* p3, float* pos, float* deriv) {

	// catmull-rom matrix
	float m[4][4] = { {-0.5f,  1.5f, -1.5f,  0.5f},
						{ 1.0f, -2.5f,  2.0f, -0.5f},
						{-0.5f,  0.0f,  0.5f,  0.0f},
						{ 0.0f,  1.0f,  0.0f,  0.0f} };
	float p[4];
	float a[4];
	for (int i = 0; i < 4; i++) {
		// Compute A = M * P
		p[0] = p0[i];
		p[1] = p1[i];
		p[2] = p2[i];
		p[3] = p3[i];
		multMatrixVector(*m, p, a);

		// Compute pos[i] = T * A
		pos[i] = pow(t, 3) * a[0] + pow(t, 2) * a[1] + t * a[2] + a[3];

		// Compute pos[i] = T * A
		deriv[i] = 3 * pow(t, 2) * a[0] + 2 * t * a[1] + a[2];
	}
}



// given  global t, returns the point in the curve
void getGlobalCatmullRomPoint(float gt, float* pos, float* deriv) {

	float t = gt * POINT_COUNT; // this is the real global t
	int index = floor(t);  // which segment
	t = t - index; // where within  the segment

	// indices store the points
	int indices[4];
	indices[0] = (index + POINT_COUNT - 1) % POINT_COUNT;
	indices[1] = (indices[0] + 1) % POINT_COUNT;
	indices[2] = (indices[1] + 1) % POINT_COUNT;
	indices[3] = (indices[2] + 1) % POINT_COUNT;

	getCatmullRomPoint(t, p[indices[0]], p[indices[1]], p[indices[2]], p[indices[3]], pos, deriv);
}

void renderCatmullRomCurve() {

	// draw curve using line segments with GL_LINE_LOOP
	float p[4][3] = { { 0.0f, 0.0f, 0.0f}, // A
					  { 1.0f, 0.0f, 1.0f}, // B
					  {-1.0f, 0.0f, 1.0f}, // C
					  {-1.0f, 0.0f, 0.0f} }; // D

	float pos[4];
	float deriv[4];
	glBegin(GL_LINE_LOOP);
	for (float gt = 0; gt < 1; gt += 0.01) {
		getGlobalCatmullRomPoint(gt, pos, deriv);
		glVertex3f(pos[0], pos[1], pos[2]);
	}
	glEnd();
}

////////////////////////////////////////////////////////////////


void lerficheiro(std::string nomeficheiro)
{
	int c = 0;
	float storefloat[3];

	std::ifstream trigsFile;
	trigsFile.open(nomeficheiro);

	std::string linha;
	if (trigsFile.is_open()) {
		while (getline(trigsFile, linha)) {
			std::string sTmp;
			for (int i = 0; i <= linha.length(); i++) {
				if (linha[i] == ' ' || linha[i] == '\0') {
					storefloat[c] = stod(sTmp);
					sTmp.clear();
					c++;
				}
				else {
					sTmp.push_back(linha[i]);
				}
			}
			Ponto aux;
			aux.x = storefloat[0];
			aux.y = storefloat[1];
			aux.z = storefloat[2];
			triangles.push_back(aux);
			c = 0;
		}
		/*
		glBindBuffer(GL_ARRAY_BUFFER, buffers);
		glBufferData(GL_ARRAY_BUFFER, vertexCount, vertexB, GL_STATIC_DRAW);
		glVertexPointer(3, GL_DOUBLE, 0, 0);
		glDrawArrays(GL_TRIANGLES, 0, vertexCount);

		glPopMatrix();
		*/
	}
}

void createVBO() {
	vertexB = (float*)malloc(sizeof(float) * triangles.size() * 3);
	int i = 0;  int vertex = 0;

	vector<Ponto>::iterator tri = triangles.begin();

	for (int j = 0; tri != triangles.end(); j++, vertex++) {
		vertexB[i] = tri[j].x;
		i++;
		vertexB[i] = tri[j].y;
		i++;
		vertexB[i] = tri[j].z;

		vertex++;
	}

	vertexCount = vertex;
}


void desenhar(void)
{
	
	std::vector<OperFile*>::iterator itout;
	itout = files.begin();

	
	while (itout != files.end()) {
		triangles.clear();
		OperFile* op = *itout;
		itout++;
		
		char* stds = op->fileName;
		std::string std(stds);
		lerficheiro(std);
		std::vector<Ponto>::iterator it;
		it = triangles.begin();

		while (it != triangles.end()) {

			glPushMatrix();
			std::vector<Oper*>::iterator it2;
			it2 = op->operations.begin();
			while(it2 != op->operations.end())
			{	
				Oper* oper = *it2;
				it2++;
				if (strcmp(oper->operation, "translate")==0)
				{
					glTranslatef(oper->x, oper->y, oper->z);
				}
				if (strcmp(oper->operation, "rotate") == 0)
				{
					glRotatef(oper->angle, oper->x, oper->y, oper->z);
				}
				if (strcmp(oper->operation, "scale") == 0) {
					glScalef(oper->x, oper->y, oper->z) ;
				}
			}
			Ponto aux_1 = *it; it++;
			Ponto aux_2 = *it; it++;
			Ponto aux_3 = *it; it++;

			glBegin(GL_TRIANGLES);
			//glColor3f(0.3, 0.1, 0);
			glColor3f(1, 1, 1);
			glVertex3d(aux_1.x, aux_1.y, aux_1.z);
			glVertex3d(aux_2.x, aux_2.y, aux_2.z);
			glVertex3d(aux_3.x, aux_3.y, aux_3.z);
			glEnd();

			glPopMatrix();
		}
		
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
	xmlParser("config.xml", files);
	createVBO();
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
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	// enter GLUT’s main cycle
	glutMainLoop();
	return 1;
}