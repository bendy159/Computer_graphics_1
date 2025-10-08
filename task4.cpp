#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include "lodepng.h"

using namespace std;

struct Vertex {
    float x, y, z;
};

vector<Vertex> read_obj_vertices() {        //Читаем Вершины
    vector<Vertex> vertices;
    ifstream file("objects/IronMan.obj");
    
    if (!file.is_open()) {
        cout << "Error with opening file" << endl;
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

void task4(const int H, const int W) {              //Не лень было реализовать динамическое масштабирование
    vector<Vertex> vertices = read_obj_vertices();

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
    
    
    for(int i = 0; i < (int)vertices.size(); i++) { // Отрисовка вершин вершины
        int x = (int)(vertices[i].x * scale + offset_x);
        int y = (int)(H - (vertices[i].y * scale + offset_y));
        
        
        if (x >= 0 && x < W && y >= 0 && y < H) {   //На всякий пожарный проверяем границы и рисуем точку
            int index = 3 * (y * W + x);
            image[index] = 0;      // Red
            image[index + 1] = 0;  // Green  
            image[index + 2] = 255; // Blue
        }
    }
    
    lodepng::encode("rabbit1.png", image, W, H, LCT_RGB, 8);
}

int main() {
    const unsigned int W = 1000;
    const unsigned int H = 1000;
    
    task4(H, W);
    
    return 0;
}