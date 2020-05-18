#ifndef CG_TP_PARSER_H
#define CG_TP_PARSER_H

#include <iostream>
#include <string.h>
#include "tinyxml2.h"
#include "operations.h"

using namespace std;
using namespace tinyxml2;

void parseGroup(const XMLElement* child, Group g, int i, vector<OperFile*>& files);

void xmlParser(const char* config, vector<OperFile*>& files);

#endif //CG_TP_PARSER_H