#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "lodepng.h"
#include <ctime>

using namespace std;

struct Vertex {
    double x, y, z;
};

struct Face {
    int v1, v2, v3;
};

struct Vec {
    double x, y, z;

    Vec() {

    }

    Vec(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

Vec LIGHT_DIR(0, 0, 1);


Vec normal(const Vertex& v0, const Vertex& v1, const Vertex& v2)  {
    Vec result;
    result.x = (v1.y - v2.y) * (v1.z - v0.z) - (v1.z - v2.z) * (v1.y - v0.y);
    result.y = -((v1.x - v2.x) * (v1.z -  v0.z) - (v1.z - v2.z) * (v1.x - v0.x));
    result.z = (v1.x - v2.x) * (v1.y - v0.y) - (v1.y - v2.y) * (v1.x - v0.x);

    return result;
}

double vec_mod(Vec& v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double scal_vec(Vec& v1, Vec& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double back_cos(Vec& norm, Vec& light) {
    double n_mod = vec_mod(norm);
    double l_mod = vec_mod(light);

    return scal_vec(norm, light) / (n_mod * l_mod);
}

vector<Vertex> read_obj_vertices() {        //Взял из 3 задачи
    vector<Vertex> vertices;
    ifstream file("objects/Formula_1.obj");
    
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
    ifstream file("objects/Formula_1.obj");
    
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

Vertex bar_cord(double x, double y, Vertex& v0, Vertex& v1, Vertex& v2) {
        
    Vertex bar_cord;
    
    double denom = (v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y);

    double lambda0 = ((x - v2.x) * (v1.y - v2.y) - (v1.x - v2.x) * (y - v2.y)) / denom;   
    double lambda1 = ((v0.x - v2.x) * (y - v2.y) - (x - v2.x) * (v0.y - v2.y)) / denom; 
    double lambda2 = 1 - lambda0 - lambda1;

    bar_cord.x = lambda0;
    bar_cord.y = lambda1;
    bar_cord.z = lambda2;

    return bar_cord;
}

void draw_triangle(vector<unsigned char>& image, double W, double H, vector<Vertex>& vertices, vector<Face>& faces, double scale, double offset_x, double offset_y) {
    
    for(int i = 0; i < faces.size(); ++i) {
        Vertex v0 = vertices[faces[i].v1 - 1];
        Vertex v1 = vertices[faces[i].v2 - 1];
        Vertex v2 = vertices[faces[i].v3 - 1];

        Vec norm = normal(v0, v1, v2);

        double dive = back_cos(norm, LIGHT_DIR);

        if(dive <= 0) {
            continue;
        }

        Vertex resize_v0, resize_v1, resize_v2;
        resize_v0.x = v0.x * scale + offset_x;
        resize_v0.y = v0.y * scale + offset_y;
        resize_v1.x = v1.x * scale + offset_x;
        resize_v1.y = v1.y * scale + offset_y;
        resize_v2.x = v2.x * scale + offset_x;
        resize_v2.y = v2.y * scale + offset_y;

        int xmin = max(0, (int)(min(min(resize_v0.x, resize_v1.x), resize_v2.x)));
        int xmax = min((int)W-1, (int)(max(max(resize_v0.x, resize_v1.x), resize_v2.x)));
        int ymin = max(0, (int)(min(min(resize_v0.y, resize_v1.y), resize_v2.y)));
        int ymax = min((int)H-1, (int)(max(max(resize_v0.y, resize_v1.y), resize_v2.y)));


        char red = 0 + rand() % 256;
        char green = 0 + rand() % 256;
        char blue = 0 + rand() % 256;
        for(int x = xmin; x <= xmax; ++x) {
            for(int y = ymin; y <= ymax; ++y) {

                Vertex bar_cords = bar_cord(x, y, resize_v0, resize_v1, resize_v2);
                if(bar_cords.x >= 0 && bar_cords.y >= 0 && bar_cords.z >= 0) {
                    int index = 3 * ((H - y) * W + x);
                    image[index] = red;     // R
                    image[index + 1] = green; // G  
                    image[index + 2] = blue; // B
                }

            }
        }
        


        
    }
}

int main() {

    srand(0);

    vector<Vertex> vertices = read_obj_vertices();
    vector<Face> faces = read_obj_faces();
    
    const unsigned int W = 10000;
    const unsigned int H = 10000;

    float min_x = vertices[0].x, max_x = vertices[0].x;     //Границы картинки
    float min_y = vertices[0].y, max_y = vertices[0].y;

    for(int i = 0; i < (int)vertices.size(); i++) {         //Ну и теперь ищем минимальную и максимальную границу. Долго но зато красиво
        if (vertices[i].x < min_x) min_x = vertices[i].x;
        if (vertices[i].x > max_x) max_x = vertices[i].x;
        if (vertices[i].y < min_y) min_y = vertices[i].y;
        if (vertices[i].y > max_y) max_y = vertices[i].y;
    }

    float scale_x = (float)(0.8 * W) / (max_x - min_x);  //Вычисляем масштаб так, чтобы отступ по 10% с каждой стороны если в центре
    float scale_y = (float)(0.8 * H) / (max_y - min_y);
    float scale = min(scale_x, scale_y);  // Сохраняем пропорции

    float offset_x = -min_x * scale + 0.1 * W;  // Смещаем в центр + отступ
    float offset_y = -min_y * scale + 0.1 * H;

    vector<unsigned char> image(H * W * 3, 0);  //Хачу черный фон
    
    draw_triangle(image, W, H, vertices, faces, scale, offset_x, offset_y);

    lodepng::encode("poly_2.png", image, W, H, LCT_RGB, 8);

    return 0;
}