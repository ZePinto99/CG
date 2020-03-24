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
#include "tinyxml.h"

float alpha = 1;
float beta = 1;
float r = 20;
double colorR = 1, colorG = 0, colorB = 1;

// Coordenadas de todos os triângulos.

struct Ponto {
	float x;
	float y;
	float z;
};

std::vector<std::string> Trabalho;
std::vector<Ponto> triangles;


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





/*
void extractPonto(std::vector<Ponto>& triangles)
{
	int cc = 0;
	float tmpPonto[3];

	//for (int i = 0; files[i] != nullptr; i++) {
		std::ifstream trigsFile;
		trigsFile.open("sphere.txt");

		std::string linha;
		if (trigsFile.is_open()) {
			while (getline(trigsFile, linha)) {
				std::string sTmp;
				for (int i = 0; i <= linha.length(); i++) {
					if (line[i] == ' ' || linha[i] == '\0') {
						tmpPonto[cc] = stod(sTmp);
						sTmp.clear();
						cc++;
					}
					else {
						sTmp.push_back(line[i]);
					}
				}
				Ponto aux = std::make_tuple(tmpPonto[0], tmpPonto[1], tmpPonto[2]);
				triangles.push_back(aux);
				cc = 0;
			}
	//	}
	}
}*/




void desenhar(void)
{
	//Trabalho.push_back("cone.txt");
	//Trabalho.push_back("plane.txt");
	std::vector<std::string>::iterator itout;
	itout = Trabalho.begin();
	while (itout != Trabalho.end()) {
		triangles.clear();
		std::string std = *itout;
		*itout++;
		lerficheiro(std);
		std::vector<Ponto>::iterator it;
		it = triangles.begin();

		while (it != triangles.end()) {


			Ponto aux_1 = *it; it++;
			Ponto aux_2 = *it; it++;
			Ponto aux_3 = *it; it++;

			glBegin(GL_TRIANGLES);
			glColor3f(1, 0, 0);
			//glColor3f(1, 1, 1);
			glVertex3d(aux_1.x, aux_1.y, aux_1.z);
			glVertex3d(aux_2.x, aux_2.y, aux_2.z);
			glVertex3d(aux_3.x, aux_3.y, aux_3.z);
			colorR -= 0.1; colorG += 0.1; colorB -= 0.1;
			glEnd();
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
	gluPerspective(45.0f, ratio, 1.0f, 1000.0f);
	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}

void renderScene(void)
{
	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set the camera
	glLoadIdentity();


	gluLookAt(r * cos(beta) * sin(alpha), r * sin(beta), r * cos(beta) * cos(alpha),
		0.0, 0.0, 0.0,
		0.0f, 1.0f, 0.0f);


	// put the geometric transformations here
	glBegin(GL_LINES);
	// X axis in red
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(-100.0f, 0.0f, 0.0f);
	glVertex3f(100.0f, 0.0f, 0.0f);
	// Y Axis in Green
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, -100.0f, 0.0f);
	glVertex3f(0.0f, 100.0f, 0.0f);
	// Z Axis in Blue
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, -100.0f);
	glVertex3f(0.0f, 0.0f, 100.0f);
	glEnd();
	// put drawing instructions here
	desenhar();
	

	// End of frame
	glutSwapBuffers();
}

void processSpecialKeys(int key, int xx, int yy)
{
	switch (key) {
	case GLUT_KEY_RIGHT:
		alpha += 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_LEFT:
		alpha -= 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_UP:
		beta += 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_DOWN:
		beta -= 0.1;
		glutPostRedisplay();
		break;
	default:
		std::cout << "Não conheço esse comando!" << "\n";
	}
}



void xmlParser(const char* configFile)
{
	TiXmlDocument doc("config.xml");
	Trabalho.clear();

	if (doc.LoadFile())
	{
		TiXmlHandle hDoc(&doc);
		TiXmlElement* pRoot, * pParm;
		pRoot = doc.FirstChildElement("scene");
		if (pRoot)
		{
			pParm = pRoot->FirstChildElement("model");
			while (pParm)
			{
				std::string filename = pParm->Attribute("file");
				std::cout << "nomes:" << "\n";
				std::cout << filename << "\n";
				Trabalho.push_back(filename);
				pParm = pParm->NextSiblingElement("model");
			}
		}
	}
	else
	{
		std::cout << "Could not load XML File. \n";
		return;
	}
}



int main(int argc, char** argv)
{
	///////////////////////////////////////////////////////////////////////////////////////

//	Ler codigo xml aqui

	xmlParser("config.xml");

	/*std::cout << "myvector contains:";
	for (std::vector<std::string>::iterator it = Trabalho.begin(); it != Trabalho.end(); ++it)
		std::cout << ' ' << *it;
	std::cout << '\n';*/


	///////////////////////////////////////////////////////////////////////////////////////
	//extractPonto(triangles);
	// put GLUT’s init here
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(800, 800);
	glutCreateWindow("CG_Trabalho_prático");
	// put callback registry here
	glutReshapeFunc(changeSize);
	glutSpecialFunc(processSpecialKeys);
	glutIdleFunc(renderScene);
	glutDisplayFunc(renderScene);
	// some OpenGL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.184314f, 0.309804f, 0.309804f, 0.0f);
	// enter GLUT’s main cycle
	glutMainLoop();
	return 1;
}