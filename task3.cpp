#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

struct Vertex {             //Уж лучше структуру точек создам, чем буду вектор в векторе писать
    double x, y, z;
};


vector<Vertex> read_obj_vertices(const string& filename) {      
    vector<Vertex> vertices;        //Как раз вектор вершин сразу
    ifstream file(filename);
    
    if (!file.is_open()) {
        cout << "Error with opening file";
        return vertices;
    }
    
    string line;        //Наша линия, в которой будут все вершина
    
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'f') continue;     //Просто немного ускоряем перебор
        
        stringstream ss(line);
        string prefix;
        ss >> prefix;       //Смотрю первый элемент. Но не через line[0], так как может потом понадобиться
        
        if (prefix == "v") {
            Vertex v;
            if (ss >> v.x >> v.y >> v.z) {      //Вот записываю три координаты
                vertices.push_back(v);
            }
        }
    }
    
    file.close();
    return vertices;
}

void task3() {
    string filename = "objects/model_1.obj";
    vector<Vertex> vertices = read_obj_vertices(filename);

    for (int i = 0; i < 5; ++i) {           //Работает ли лаба вообще
        cout << "v " << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << endl;
    }
}

int main() {
    const unsigned int W = 200;
    const unsigned int H = 200;
    
    task3();
    
    return 0;
}