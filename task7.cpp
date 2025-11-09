#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "lodepng.h"

using namespace std;

struct Vertex {
    float x, y, z;
};

struct Face {
    int v1, v2, v3;
};


vector<Vertex> read_obj_vertices() {        //Взял из 3 задачи
    vector<Vertex> vertices;
    ifstream file("objects/model_1.obj");
    
    if (!file.is_open()) {
        cout << "Error with open file" << endl;
        return vertices;
    }
    
    string line;
    
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'f') continue;
        
        stringstream ss(line);
        string prefix;
        ss >> prefix;
        
        if (prefix == "v") {
            Vertex v;
            if (ss >> v.x >> v.y >> v.z) {
                vertices.push_back(v);
            }
        }
    }
    
    file.close();
    return vertices;
}

vector<Face> read_obj_faces() {     //Взял из 5 задачи
    vector<Face> faces;
    ifstream file("objects/model_1.obj");
    
    if (!file.is_open()) {
        cout << "Ошибка открытия файла" << endl;
        return faces;
    }
    
    string line;
    
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'v') continue;
        
        stringstream ss(line);
        string prefix;
        ss >> prefix;
        
        if (prefix == "f") {
            Face face;
            string v1_str, v2_str, v3_str;
            
            ss >> v1_str >> v2_str >> v3_str;
            
            stringstream ss1(v1_str);
            stringstream ss2(v2_str);
            stringstream ss3(v3_str);
            
            string token;
            
            getline(ss1, token, '/');
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

vector<Vertex> bar_cord(double x, double y) {

    vector<Vertex> bar_cords;

    vector<Vertex> vertices = read_obj_vertices();
    vector<Face> faces = read_obj_faces();
    
    vector<Vertex> bar_cord;

    for(int i = 0; i < faces.size(); i++) {
        
        Vertex bar_cord;

        Vertex& v0 = vertices[faces[i].v1 - 1];
        Vertex& v1 = vertices[faces[i].v2 - 1];         //Индексация с 1, но нам надо с 0
        Vertex& v2 = vertices[faces[i].v3 - 1];
        
        double lambda0 = ((x - v2.x) * (v1.y - v2.y) - (v1.x - v2.x) * (y - v2.y)) / 
                        ((v0.x - v2.x) * (v1.y - v2.y) - (v1.x - v2.x) * (v0.y - v2.y));
        
        double lambda1 = ((v0.x - v2.x) * (y - v2.y) - (x - v2.x) * (v0.y - v2.y)) / 
                        ((v0.y - v2.y) * (v1.x - v2.x) - (v0.x - v2.x) * (v1.y - v2.y));
        
        double lambda2 = 1 - lambda0 - lambda1;

        bar_cord.x = lambda0;
        bar_cord.y = lambda1;
        bar_cord.z = lambda2;

        bar_cords.push_back(bar_cord);
    }

    return bar_cords;
}

int main() {
    vector<Vertex> vertices = read_obj_vertices();
    vector<Face> faces = read_obj_faces();


    for(int i = 0; i < faces.size(); ++i) {
        Vertex& v1 = vertices[faces[i].v1 - 1];
        Vertex& v2 = vertices[faces[i].v2 - 1];
        Vertex& v3 = vertices[faces[i].v3 - 1];
    }
    
    
    return 0;
}