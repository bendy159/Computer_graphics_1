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
    ifstream file("objects/Formula_ferrariloshary.obj");
    
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
    ifstream file("objects/Formula_ferrariloshary.obj");
    
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

void alg_bresenhema(vector<unsigned char>& image, unsigned int W, unsigned int H,       //Взял из второй задачи
                   int x0, int y0, int x1, int y1,
                   unsigned char r, unsigned char g, unsigned char b) {
    
    if(x0 > x1) {
        swap(x0, x1);
        swap(y0, y1);
    }

    bool xChange = false;
    
    if(abs((int)(x0 - x1)) < abs((int)(y0 - y1))) {
        swap(x0, y0);
        swap(x1, y1);
        xChange = true;

        if(x0 > x1) {
            swap(x0, x1);
            swap(y0, y1);
        }
    }

    unsigned int y = y0;
    double dy = 2 * abs((int)(y1 - y0));
    double derror = 0.0;
    int y_update = (y1 > y0) ? 1 : -1;

    for (int x = x0; x <= x1; ++x) {
        if (xChange) {
            if (y < W && x < (int)H) {
                image[3 * (x * W + y)] = r;
                image[3 * (x * W + y) + 1] = g;
                image[3 * (x * W + y) + 2] = b;
            }
        }
        else {
            if (x < (int)W && y < H) {
                image[3 * (y * W + x)] = r;
                image[3 * (y * W + x) + 1] = g;
                image[3 * (y * W + x) + 2] = b;
            }
        }

        derror += dy;
        if (derror > (x1 - x0)) {
            derror -= 2 * (x1 - x0);
            y += y_update;
        }
    }
}

void task6(const int H, const int W) {

    vector<Vertex> vertices = read_obj_vertices();
    vector<Face> faces = read_obj_faces();


    float min_x = vertices[0].x, max_x = vertices[0].x;
    float min_y = vertices[0].y, max_y = vertices[0].y;
    
    for(int i = 0; i < (int)vertices.size(); i++) {                 //Взял из 4 задачи
        if (vertices[i].x < min_x) min_x = vertices[i].x;
        if (vertices[i].x > max_x) max_x = vertices[i].x;
        if (vertices[i].y < min_y) min_y = vertices[i].y;
        if (vertices[i].y > max_y) max_y = vertices[i].y;
    }
    

    float scale_x = (float)(0.8 * W) / (max_x - min_x);
    float scale_y = (float)(0.8 * H) / (max_y - min_y);
    float scale = min(scale_x, scale_y);
    
    float offset_x = -min_x * scale + 100;
    float offset_y = -min_y * scale + 100;
    

    vector<unsigned char> image(H * W * 3, 0);
    

    for(int i = 0; i < (int)faces.size(); i++) {

        Vertex& v1 = vertices[faces[i].v1 - 1];
        Vertex& v2 = vertices[faces[i].v2 - 1];         //Индексация с 1, но нам надо с 0
        Vertex& v3 = vertices[faces[i].v3 - 1];
        

        int x1 = (int)round(v1.x * scale + offset_x);           //Масштабирование
        int y1 = (int)round(H - (v1.y * scale + offset_y));
        int x2 = (int)round(v2.x * scale + offset_x);
        int y2 = (int)round(H - (v2.y * scale + offset_y));
        int x3 = (int)round(v3.x * scale + offset_x);
        int y3 = (int)round(H - (v3.y * scale + offset_y));
        
        alg_bresenhema(image, W, H, x1, y1, x2, y2, 0, 0, 255);         //Строим треугольник
        alg_bresenhema(image, W, H, x2, y2, x3, y3, 0, 0, 255);  
        alg_bresenhema(image, W, H, x3, y3, x1, y1, 0, 0, 255);  
    }
    
    lodepng::encode("rabbit2.png", image, W, H, LCT_RGB, 8);
}

int main() {
    const unsigned int W = 10000;
    const unsigned int H = 10000;
    
    task6(H, W);
    
    return 0;
}