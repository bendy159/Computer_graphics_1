#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

struct Face {               //Такая же ситуация, что с вершинами
    int v1, v2, v3;
};


vector<Face> read_obj_faces() {     //Считываем полигоны
    vector<Face> faces;
    ifstream file("objects/model_1.obj");
    
    if (!file.is_open()) {
        cout << "Error with opening file\n";
        return faces;
    }
    
    string line;
    
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'v') continue;     //Ультра оптимизация
        
        stringstream ss(line);      //Так же префикс смотрю
        string prefix;
        ss >> prefix;
        
        if (prefix == "f") {
            Face face;
            string v1_str, v2_str, v3_str;
            
            ss >> v1_str >> v2_str >> v3_str;       //Три тройки записываю
            
            stringstream ss1(v1_str);
            stringstream ss2(v2_str);
            stringstream ss3(v3_str);
            
            string token;
            
            getline(ss1, token, '/');       //И затем считываю строку до первого /
            face.v1 = stoi(token);
            
            getline(ss2, token, '/');
            face.v2 = stoi(token);
            
            getline(ss3, token, '/');
            face.v3 = stoi(token);
            
            faces.push_back(face);
        }
    }
    
    file.close();
    return faces;
}

void task5() {
    vector<Face> faces = read_obj_faces();
    
    if (faces.empty()) {
        cout << "Error with open file" << endl;
        return;
    }
    
    for (int i = 0; i < 10; i++) {      //Проверка
        cout << "f " << faces[i].v1 << " " << faces[i].v2 << " " << faces[i].v3 << endl;
    }
    
}

int main() {
    task5();
    return 0;
}