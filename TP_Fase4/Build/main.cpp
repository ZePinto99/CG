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
std::vector<Ponto> normal;
std::vector<Ponto> texture;
vector<Light*> lightVector;

vector<OperFile*> files; // Vector de OperFiles (que relacionam os ficheiros
						 // com as suas respetivas transforma��es).

int* fiVertexCount;
GLuint vertexCount;
GLuint buffers[20];
GLuint normals[20];
double** vertexB;
double** normais;
double** textures;

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
	int p = 0, n = 0, t = 0;  int vertex = 0;
	
	vector<Ponto>::iterator tri = triangles.begin();
	vector<Ponto>::iterator norm = normal.begin();
	vector<Ponto>::iterator text = texture.begin();
	//cout << texture.size() << "!!";

	vertexB[i]    = (double*)malloc(sizeof(double) * triangles.size() * 3);
    normais[i]    = (double*)malloc(sizeof(double) * normal.size() * 3);
  	textures[i]   = (double*)malloc(sizeof(double) * texture.size() * 2);
	
	printf("%d\n", normal.size());

	while (tri != triangles.end()) {

		Ponto aux_1 = *tri; tri++;
		Ponto aux_2 = *tri; tri++;
		Ponto aux_3 = *tri; tri++;

		Ponto norm_1 = *norm; norm++;
		Ponto norm_2 = *norm; norm++;
		Ponto norm_3 = *norm; norm++;
		Ponto text_1 = *text; text++;
		Ponto text_2 = *text; text++;
     	Ponto text_3 = *text; text++; 

		vertexB[i][p] = aux_1.x;
		vertexB[i][p + 1] = aux_1.y;
		vertexB[i][p + 2] = aux_1.z;

		normais[i][n] = norm_1.x;
		normais[i][n + 1] = norm_1.y;
		normais[i][n + 2] = norm_1.z;
		
		textures[i][t] = 1.f;
		textures[i][t + 1] = text_1.y;
		textures[i][t + 1] = text_1.z;
		vertex++;

		vertexB[i][p + 3] = aux_2.x;
		vertexB[i][p + 4] = aux_2.y;
		vertexB[i][p + 5] = aux_2.z;

		
		normais[i][n + 3] = norm_1.x;
		normais[i][n + 4] = norm_1.y;
		normais[i][n + 5] = norm_1.z;
		
		textures[i][t + 2] = text_2.x;
		textures[i][t + 3] = text_2.y;
		vertex++;

		vertexB[i][p + 6] = aux_3.x;
		vertexB[i][p + 7] = aux_3.y;
		vertexB[i][p + 8] = aux_3.z;

		
		normais[i][n + 6] = norm_1.x;
		normais[i][n + 7] = norm_1.y;
		normais[i][n + 8] = norm_1.z;
		
		textures[i][t + 4] = text_3.x;
		textures[i][t + 5] = text_3.y;
		vertex++;

		p += 9; n += 9; t += 6;
	}
	vertexCount = vertex;
	normal.end();
	texture.end();
	triangles.end();
}

void lerficheiro(string nomeficheiro)
{
	int c = 0;
	cout << nomeficheiro << "\n";

	double storedouble[3];
	for (int i = 0; i < 3; i++) storedouble[i] = 0;

	ifstream trigsFile;
	trigsFile.open(nomeficheiro);
	int type = 0;
	string linha;
	if (trigsFile.is_open()) {
		int numlinhas = 0;
		while (getline(trigsFile, linha)) {
			type++;
			string sTmp;
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

			if (type == 1 || type == 2 || type == 3) {
				////cout << "1 \n";
				triangles.push_back(aux);
				////cout << "1.1 \n";
			}
			if (type == 4 || type == 5 || type == 6)
			{
				////cout << "2 \n";
				normal.push_back(aux);
				////cout << "2.1 \n";
			}
			if (type == 7 || type == 8 || type == 9)
			{
				////cout << "3 \n";
				texture.push_back(aux);
				////cout << "3.1 \n";
			}
			if (type == 9)
			{
				////cout << "4 \n";
				type = 0;
				////cout << "4.1 \n";
			}
			c = 0;
			numlinhas++;
		}
		cout << numlinhas << "\n";
	}
}


void lertudoemaisalgumacoisa() {
	vertexB = (double**)malloc(sizeof(double*) * files.size());
	normais = (double**)malloc(sizeof(double*) * files.size());
	textures = (double**)malloc(sizeof(double*) * files.size());
	std::vector<OperFile*>::iterator itout;
	itout = files.begin();
	int i = 0;
	while (itout != files.end()) {
		triangles.clear();
		normal.clear();
		texture.clear();
		//cout << texture.size() << "........";
		OperFile* op = *itout;
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
		getGlobalCatmullRomPoint((double)j / n, pos, deriv, curve, size);
		glVertex3d(pos[0], pos[1], pos[2]);
	}
	glEnd();

	getGlobalCatmullRomPoint(intervalos[i], pos, deriv, curve, size);

	//if (time < files.size()) {
	double ttt = glutGet(GLUT_ELAPSED_TIME);
	checkpoints[i] = 1 / (t->time * 1000);
	time++;
	//}
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

void processLight() {
	float pos[4] = { 1.0, 1.0, 1.0, 0.0 };
	GLfloat ambi[4] = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat diff[4] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat direction[3] = { 0.0, 0.0, -1.0 };
	vector<Light*>::iterator it;
	for (it = lightVector.begin(); it != lightVector.end(); it++) {
		Light* light = *it;
		if (strcmp(light->type, "POINT") == 0) {
			pos[0] = static_cast<GLfloat>(light->x);
			pos[1] = static_cast<GLfloat>(light->y);
			pos[2] = static_cast<GLfloat>(light->z);
			pos[3] = 1;
			glEnable(GL_LIGHT0);
			glLightfv(GL_LIGHT0, GL_POSITION, pos);
			glLightfv(GL_LIGHT0, GL_AMBIENT, ambi);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diff);
		}
		else if (strcmp(light->type, "DIRECTIONAL") == 0) {
			pos[0] = static_cast<GLfloat>(light->x);
			pos[1] = static_cast<GLfloat>(light->y);
			pos[2] = static_cast<GLfloat>(light->z);
			pos[3] = 0;
			glEnable(GL_LIGHT0);
			glLightfv(GL_LIGHT0, GL_POSITION, pos);
			glLightfv(GL_LIGHT0, GL_AMBIENT, ambi);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diff);
		}
		else if (strcmp(light->type, "SPOT") == 0) {
			pos[0] = static_cast<GLfloat>(light->x);
			pos[1] = static_cast<GLfloat>(light->y);
			pos[2] = static_cast<GLfloat>(light->z);
			pos[3] = 1;
			direction[0] = static_cast<GLfloat>(light->x);
			direction[1] = static_cast<GLfloat>(light->y);
			direction[2] = static_cast<GLfloat>(light->z);
			glEnable(GL_LIGHT0);
			glLightfv(GL_LIGHT0, GL_POSITION, pos);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diff);
			glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, direction);
			glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 45.0);
			glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 0.0);
		}
	}
	cout << "processLight check \n";
}

void processColor(Color* color)
{
	GLfloat diff[4] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat spec[4] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat emis[4] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat ambi[4] = { 0.2, 0.2, 0.2, 1.0 };
	if (strcmp(color->component, "diffuse") == 0) {
		diff[0] = static_cast<GLfloat>(color->r);
		diff[1] = static_cast<GLfloat>(color->g);
		diff[2] = static_cast<GLfloat>(color->b);
		diff[3] = 1;
		glMaterialfv(GL_FRONT, GL_DIFFUSE, diff);
	}
	if (strcmp(color->component, "specular") == 0) {
		spec[0] = static_cast<GLfloat>(color->r);
		spec[1] = static_cast<GLfloat>(color->g);
		spec[2] = static_cast<GLfloat>(color->b);
		spec[3] = 1;
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
	}
	if (strcmp(color->component, "emissive") == 0) {
		emis[0] = static_cast<GLfloat>(color->r);
		emis[1] = static_cast<GLfloat>(color->g);
		emis[2] = static_cast<GLfloat>(color->b);
		emis[3] = 1;
		glMaterialfv(GL_FRONT, GL_EMISSION, emis);
	}
	if (strcmp(color->component, "ambient") == 0) {
		ambi[0] = static_cast<GLfloat>(color->r);
		ambi[1] = static_cast<GLfloat>(color->g);
		ambi[2] = static_cast<GLfloat>(color->b);
		ambi[3] = 1;
		glMaterialfv(GL_FRONT, GL_AMBIENT, ambi);
	}
	cout << "processColor check \n";
}

/*
int loadTexture(string texture) {
	unsigned int t, tw, th;
	unsigned char* texData;
	unsigned int id;
	ilInit();
	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_LOWER_LEFT);
	ilGenImages(1, &t);
	ilBindImage(t);
	ilLoadImage((ILstring)texture.c_str());
	tw = ilGetInteger(IL_IMAGE_WIDTH);
	th = ilGetInteger(IL_IMAGE_HEIGHT);
	ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE);
	texData = ilGetData();
	glGenTextures(1, &id);
	glBindTexture(GL_TEXTURE_2D, id);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tw, th, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
	glGenerateMipmap(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	return id;
	cout << "loadTexture check \n";
}

void initTexturesByID() {
	int i = 0;
	texturesByID = (GLuint*)malloc(sizeof(GLuint) * files.size());
	for (int i = 0; i < files.size(); i++) texturesByID[i] = NULL;
	vector<OperFile*>::iterator it;
	for (it = files.begin(); it != files.end(); i++, it++) {
		OperFile* oper = *it;
		if (oper->texture != NULL) {
			char dir[60] = "../texturas/";
			char* textureFileName = strcat(dir, oper->texture);
			string textureString(textureFileName);
			texturesByID[i] = loadTexture(textureString);
		}
	}
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_TEXTURE_2D);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	cout << "initTexturesByID check \n";
}
*/
void desenhar(void)
{
	vector<OperFile*>::iterator itout;
	itout = files.begin();
	int i = 0;

	while (itout != files.end()) {
		triangles.clear();
		OperFile* op = *itout;
		itout++;

		char* stds = op->fileName;
		glPushMatrix();

		vector<Oper*>::iterator it2;
		it2 = op->operations.begin();
		while (it2 != op->operations.end())
		{
			Oper* oper = *it2;
			it2++;
			if (strcmp(oper->operation, "translate") == 0)
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
				glScalef(oper->x, oper->y, oper->z);
			}
		}
		glColor3f(1, 1, 1);

//		glBindTexture(GL_TEXTURE_2D, textures[i][i]);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[i]);
		glBufferData(GL_ARRAY_BUFFER, vertexCount * sizeof(double) * 3, vertexB[i], GL_STATIC_DRAW);
		glVertexPointer(3, GL_DOUBLE, 0, 0);
		glDrawArrays(GL_TRIANGLES, 0, vertexCount);

		glBindBuffer(GL_ARRAY_BUFFER, normals[0]);
		glBufferData(GL_ARRAY_BUFFER, vertexCount * 8 * 3, normais[i], GL_STATIC_DRAW);
		glPopMatrix();
		i++;
	}

}

void changeSize(int w, int h)
{
	// Prevent a divide by zero, when window is too short
	// (you can�t make a window with zero width).
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
	glGenBuffers(size, normals);
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
		std::cout << "N�o conhe�o esse comando!" << "\n";
	}
}

void init() {
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
}



int main(int argc, char** argv)
{
	// put GLUT�s init here
	xmlParser("sistemaSolarDinamico.xml", files, lightVector);

	lertudoemaisalgumacoisa();
	init();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(1920, 1080);
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
	// enter GLUT�s main cycle
	glutMainLoop();
	return 1;
}
