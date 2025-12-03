#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "lodepng.h"
#include <ctime>
#include <limits>

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
double Z_BUF = numeric_limits<double>::max();

vector<vector<double>> matrix_mult(vector<vector<double>> first, vector<vector<double>> second) {
    vector<vector<double>> result(3, vector<double>(3, 0));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                result[i][j] += first[i][k] * second[k][j];
            }
        }
    }
    return result;
}

Vertex matrix_vector_mult(vector<vector<double>> matrix, Vertex vertex) {
    Vertex result;
    result.x = matrix[0][0] * vertex.x + matrix[0][1] * vertex.y + matrix[0][2] * vertex.z;
    result.y = matrix[1][0] * vertex.x + matrix[1][1] * vertex.y + matrix[1][2] * vertex.z;
    result.z = matrix[2][0] * vertex.x + matrix[2][1] * vertex.y + matrix[2][2] * vertex.z;
    return result;
}

vector<vector<double>> rotation_matrix_x(double alpha) {
    vector<vector<double>> rx = {
        {1, 0, 0},
        {0, cos(alpha), sin(alpha)},
        {0, -sin(alpha), cos(alpha)}
    };
    return rx;
}

vector<vector<double>> rotation_matrix_y(double beta) {
    vector<vector<double>> ry = {
        {cos(beta), 0, sin(beta)},
        {0, 1, 0},
        {-sin(beta), 0, cos(beta)}
    };
    return ry;
}

vector<vector<double>> rotation_matrix_z(double gamma) {
    vector<vector<double>> rz = {
        {cos(gamma), sin(gamma), 0},
        {-sin(gamma), cos(gamma), 0},
        {0, 0, 1}
    };
    return rz;
}

void center_model(vector<Vertex>& vertices) {
    double center_x = 0;
    double center_y = 0;
    double center_z = 0;

    for(Vertex& v : vertices) {
        center_x += v.x;
        center_y += v.y;
        center_z += v.z;
    }

    center_x /= vertices.size();
    center_y /= vertices.size();
    center_z /= vertices.size();

    for(Vertex& v : vertices) {
        v.x = v.x - center_x;
        v.y = v.y - center_y;
        v.z = v.z - center_z;
    }
}

vector<Vertex> transform_vertices(vector<Vertex> vertices, double alpha, double beta, double gamma, double tx, double ty, double tz) {

    vector<vector<double>> rx = rotation_matrix_x(alpha);
    vector<vector<double>> ry = rotation_matrix_y(beta);
    vector<vector<double>> rz = rotation_matrix_z(gamma);
    
    vector<vector<double>> r_xy = matrix_mult(rx, ry);
    vector<vector<double>> r = matrix_mult(r_xy, rz);
    
    vector<Vertex> transformed_vertices;
    for (Vertex v : vertices) {
        Vertex rotated = matrix_vector_mult(r, v);
        Vertex transformed;
        transformed.x = rotated.x + tx;
        transformed.y = rotated.y + ty;
        transformed.z = rotated.z + tz;
        transformed_vertices.push_back(transformed);
    }
    
    return transformed_vertices;
}


Vec normal(const Vertex& v0, const Vertex& v1, const Vertex& v2)  {
    Vec result;
    result.x = (v1.y - v2.y) * (v1.z - v0.z) - (v1.z - v2.z) * (v1.y - v0.y);
    result.y = -((v1.x - v2.x) * (v1.z -  v0.z) - (v1.z - v2.z) * (v1.x - v0.x));
    result.z = (v1.x - v2.x) * (v1.y - v0.y) - (v1.y - v2.y) * (v1.x - v0.x);

    return result;
}

double compute_z_coordinate(const Vertex& bar_cords, const Vertex& v0, const Vertex& v1, const Vertex& v2) {        //Глубина точки
    return bar_cords.x * v0.z + bar_cords.y * v1.z + bar_cords.z * v2.z;
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

vector<Vertex> read_obj_vertices(string filename) {        //Взял из 3 задачи
    vector<Vertex> vertices;
    ifstream file(filename);
    
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

vector<Face> read_obj_faces(string filename) {     //Взял из 5 задачи
    vector<Face> faces;
    ifstream file(filename);
    
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

vector<Vec> compute_v_norm(vector<Vertex>& veritices, vector<Face> faces) {
    vector<Vec> v_n(veritices.size(), Vec(0, 0, 0));
    for(Face face : faces) {
        Vertex v0 = veritices[face.v1 - 1];
        Vertex v1 = veritices[face.v2 - 1];
        Vertex v2 = veritices[face.v3 - 1];

        Vec norm = normal(v0, v1, v2);

        v_n[face.v1 - 1].x += norm.x;
        v_n[face.v1 - 1].y += norm.y;
        v_n[face.v1 - 1].z += norm.z;

        v_n[face.v2 - 1].x += norm.x;
        v_n[face.v2 - 1].y += norm.y;
        v_n[face.v2 - 1].z += norm.z;
        
        v_n[face.v3 - 1].x += norm.x;
        v_n[face.v3 - 1].y += norm.y;
        v_n[face.v3 - 1].z += norm.z;

    }
    for(Vec& vert : v_n) {
        double length = sqrt(vert.x * vert.x + vert.y * vert.y + vert.z * vert.z);
        if(length > 0) {
            vert.x /= length;
            vert.y /= length;
            vert.z /= length;
        }
    }
    return v_n;
}

Vec light(vector<Vec>& v_n, Face& face) {
    Vec n0 = v_n[face.v1 - 1];
    Vec n1 = v_n[face.v2 - 1];
    Vec n2 = v_n[face.v3 - 1];

    double i0 = scal_vec(n0, LIGHT_DIR);
    double i1 = scal_vec(n1, LIGHT_DIR);
    double i2 = scal_vec(n2, LIGHT_DIR);

    Vec res;
    res.x = i0;
    res.y = i1;
    res.z = i2;

    return res;
}

void draw_triangle(vector<unsigned char>& image, vector<double>& z_buffer, double W, double H,
    vector<Vertex>& vertices, vector<Face>& faces, double a_X, double a_Y, double u_0, double v_0) {
    vector<Vec> v_n = compute_v_norm(vertices, faces);
    for(int i = 0; i < faces.size(); ++i) {
        Vertex v0 = vertices[faces[i].v1 - 1];
        Vertex v1 = vertices[faces[i].v2 - 1];
        Vertex v2 = vertices[faces[i].v3 - 1];

        Vec norm = normal(v0, v1, v2);

        double dive = back_cos(norm, LIGHT_DIR);

        if(dive >= 0) {
            continue;
        }

        Vertex resize_v0, resize_v1, resize_v2;
        
        resize_v0.x = (a_X * v0.x / v0.z) + u_0;
        resize_v0.y = (a_Y * v0.y / v0.z) + v_0;
        resize_v0.z = v0.z;

        resize_v1.x = (a_X * v1.x / v1.z) + u_0;
        resize_v1.y = (a_Y * v1.y / v1.z) + v_0;
        resize_v1.z = v1.z;

        resize_v2.x = (a_X * v2.x / v2.z) + u_0;
        resize_v2.y = (a_Y * v2.y / v2.z) + v_0;
        resize_v2.z = v2.z;

        int xmin = max(0, (int)(min(min(resize_v0.x, resize_v1.x), resize_v2.x)));
        int xmax = min((int)W-1, (int)(max(max(resize_v0.x, resize_v1.x), resize_v2.x)));
        int ymin = max(0, (int)(min(min(resize_v0.y, resize_v1.y), resize_v2.y)));
        int ymax = min((int)H-1, (int)(max(max(resize_v0.y, resize_v1.y), resize_v2.y)));

        Vec l = light(v_n, faces[i]);
        cout << l.x << ' ' << l.y << l.z << '\n';
        for(int x = xmin; x <= xmax; ++x) {
            for(int y = ymin; y <= ymax; ++y) {

                Vertex bar_cords = bar_cord(x, y, resize_v0, resize_v1, resize_v2);
                if(bar_cords.x >= 0 && bar_cords.y >= 0 && bar_cords.z >= 0) {
                    double z_value = compute_z_coordinate(bar_cords, v0, v1, v2);
                    
                    int index = 3 * ((H - y) * W + x);
                    char red = -255.0 * (l.x * bar_cords.x + l.y * bar_cords.y + l.z * bar_cords.z);
                    char green = 0;
                    char blue = 0;

                    if(z_value < z_buffer[index / 3]) {
                        z_buffer[index / 3] = z_value;
                        image[index] = red;     // R
                        image[index + 1] = green; // G  
                        image[index + 2] = blue; // B

                    }
                }
            }
        }
    }
}

int main() {

    string filename = "objects/couch.obj";

    vector<Vertex> vertices = read_obj_vertices(filename);
    vector<Face> faces = read_obj_faces(filename);
    
    const unsigned int W = 10000;
    const unsigned int H = 10000;

    center_model(vertices);

    double a_X = 1000;
    double a_Y = 1000;
    double u_0 = (double)W / 2.0;
    double v_0 = (double)H / 2.0;

    double alpha = 0;  // Поворот вокруг оси X
    double beta = 0.5 + 3.14;   // Поворот вокруг оси Y
    double gamma = 0.1;  // Поворот вокруг оси Z

    vector<Vertex> transformed_vertices = transform_vertices(vertices, alpha, beta, gamma, 0, 0, 10000);

    vector<unsigned char> image(H * W * 3, 0);  //Хачу черный фон
    
    vector<double> z_buffer(H * W, Z_BUF);

    draw_triangle(image, z_buffer, W, H, transformed_vertices, faces, a_X, a_Y, u_0, v_0);

    lodepng::encode("poly_12.png", image, W, H, LCT_RGB, 8);

    return 0;
}