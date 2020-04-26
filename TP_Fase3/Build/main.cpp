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
		//	}
	}
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