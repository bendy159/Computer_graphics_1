#define _USE_MATH_DEFINES

#include "lodepng.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void dotted_line(vector<unsigned char>& image, unsigned int W, unsigned int H, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, unsigned int count, unsigned int r, unsigned int g, unsigned int b) {
    double step = 1.0 / count;          //самая простая линия, в которой интерполируется и вычисляется при помощи константы

    for (double t = 0; t < 1.0; t += step) {
        unsigned int x = round((1.0 - t) * x0 + t * x1);        //Плавно перемещается от первой точки ко второй
        unsigned int y = round((1.0 - t) * y0 + t * y1);

        if (x < W && y < H) {
            image[3 * (y * W + x)] = r;
            image[3 * (y * W + x) + 1] = g;
            image[3 * (y * W + x) + 2] = b;
        }

    }
}

void dotted_line2(vector<unsigned char>& image, unsigned int W, unsigned int H, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, unsigned int r, unsigned int g, unsigned int b) {
    double count = sqrt(pow(abs((int)(x0 - x1)), 2) + pow(abs((int)(y0 - y1)), 2));
    double step = 1.0 / count;      //То же самое что и первое, но константа определяется количеством некоторых отрезков между двумя точками

    for (double t = 0; t < 1.0; t += step) {           
        unsigned int x = round((1.0 - t) * x0 + t * x1);
        unsigned int y = round((1.0 - t) * y0 + t * y1);

        if (x < W && y < H) {
            image[3 * (y * W + x)] = r;
            image[3 * (y * W + x) + 1] = g;
            image[3 * (y * W + x) + 2] = b;
        }

    }
}

void loop_line(vector<unsigned char>& image, unsigned int W, unsigned int H, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, unsigned int r, unsigned int g, unsigned int b) {
                                                    //Уже другой принцип. идём по иксу и не используем тогда константу
    for (unsigned int x = x0; x <= x1; ++x) {
        double t = (double)(x - x0) / (double)(x1 - x0);
        unsigned int y = round((1.0 - t) * y0 + t*y1);


        if (x < W && y < H) {
            image[3 * (y * W + x)] = r;
            image[3 * (y * W + x) + 1] = g;
            image[3 * (y * W + x) + 2] = b;
        }
    }
    
}

void loop_line2(vector<unsigned char>& image, unsigned int W, unsigned int H, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, unsigned int r, unsigned int g, unsigned int b) {
                    //Устраняем проблему, когда x0 < x1
   if(x0 > x1) {
        swap(x0, x1);
        swap(y0, y1);
    }

    for (unsigned int x = x0; x <= x1; ++x) {           
        double t = (double)(x - x0) / (double)(x1 - x0);
        unsigned int y = round((1.0 - t) * y0 + t*y1);


        if (x < W && y < H) {
            image[3 * (y * W + x)] = r;
            image[3 * (y * W + x) + 1] = g;
            image[3 * (y * W + x) + 2] = b;
        }
    }
    
}

void loop_line3(vector<unsigned char>& image, unsigned int W, unsigned int H, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, unsigned int r, unsigned int g, unsigned int b) {
                        //Если движение по оси y больше, то будет разрежение точек
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
    }               //Двойной свап нужен для того, если y0 > y1
    }

    for (unsigned int x = x0; x <= x1; ++x) {
        double t = (double)(x - x0) / (double)(x1 - x0);
        unsigned int y = round((1.0 - t) * y0 + t * y1);

        if (xChange) {
            if (y < W && x < H) {
                image[3 * (x * W + y)] = r;
                image[3 * (x * W + y) + 1] = g;
                image[3 * (x * W + y) + 2] = b;
            }
        }
        else {

            if (x < W && y < H) {
                image[3 * (y * W + x)] = r;
                image[3 * (y * W + x) + 1] = g;
                image[3 * (y * W + x) + 2] = b;
            }
        }
    }

    
}

void loop_line4(vector<unsigned char>& image, unsigned int W, unsigned int H, 
                unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, 
                unsigned int r, unsigned int g, unsigned int b) {
    
                                //Делаем то же самое, но умножаем на 2*(х1 - х0) для того, чтобы сократить дроби и не использовать их в вычислениях
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
    double dy = 2 * (x1 - x0) * (double)abs((int)y1 - (int)y0) / (x1 - x0);
    double derror = 0.0;                //Ошибка нужна для того, чтобы понимать когда переходить на следующий пиксель вверх
    int y_update = (y1 > y0) ? 1 : -1;      //Либо вверх, либо вниз

    for (unsigned int x = x0; x <= x1; ++x) {
        if (xChange) {
            if (y < W && x < H) {
                image[3 * (x * W + y)] = r;
                image[3 * (x * W + y) + 1] = g;
                image[3 * (x * W + y) + 2] = b;
            }
        }
        else {
            if (x < W && y < H) {
                image[3 * (y * W + x)] = r;
                image[3 * (y * W + x) + 1] = g;
                image[3 * (y * W + x) + 2] = b;
            }
        }

        derror += dy;       //Как раз делаем шаг по игрик, а по х он 1
        if (derror > 2 * (x1 - x0) * 0.5) {
            derror -= 2 * (x1 - x0) * 1.0;
            y += y_update;
        }
    }
}

void alg_bresenhema(vector<unsigned char>& image, unsigned int W, unsigned int H, 
                unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, 
                unsigned int r, unsigned int g, unsigned int b) {
                            //Ну и всё сократив получается алгоритм Брезенхема, в котором нет дробей
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

    for (unsigned int x = x0; x <= x1; ++x) {
        if (xChange) {
            if (y < W && x < H) {
                image[3 * (x * W + y)] = r;
                image[3 * (x * W + y) + 1] = g;
                image[3 * (x * W + y) + 2] = b;
            }
        }
        else {
            if (x < W && y < H) {
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

void task2(unsigned int W, unsigned int H) {
    vector<unsigned char> image(H * W * 3, 255);
    unsigned int x0 = 100;
    unsigned int y0 = 100;

    /*for (int i = 0; i < 13; ++i) {
        double alpha = 2.0 * M_PI * i /  13;
        unsigned int x = (unsigned int)(x0 + 95 * cos(alpha));
        unsigned int y = (unsigned int)(y0 + 95 * sin(alpha));

        dotted_line(image, W, H, x0, y0, x, y, 200, 0, 0, 0);
    }

    lodepng::encode("star1.png", image, W, H, LCT_RGB, 8);*/

    /*for (int i = 0; i < 13; ++i) {
        double alpha = 2.0 * M_PI * i /  13;
        unsigned int x = (unsigned int)(x0 + 95 * cos(alpha));
        unsigned int y = (unsigned int)(y0 + 95 * sin(alpha));
        dotted_line2(image, W, H, x0, y0, x, y, 0, 0, 0);
    }
    lodepng::encode("star2.png", image, W, H, LCT_RGB, 8);*/


    
    /*for (int i = 0; i < 13; ++i) {
        double alpha = 2.0 * M_PI * i /  13;
        unsigned int x = (unsigned int)(x0 + 95 * cos(alpha));
        unsigned int y = (unsigned int)(y0 + 95 * sin(alpha));
        loop_line(image, W, H, x0, y0, x, y, 0, 0, 0);
    }
    lodepng::encode("star3.png", image, W, H, LCT_RGB, 8);*/

    /*for (int i = 0; i < 13; ++i) {
        double alpha = 2.0 * M_PI * i /  13;
        unsigned int x = (unsigned int)(x0 + 95 * cos(alpha));
        unsigned int y = (unsigned int)(y0 + 95 * sin(alpha));
        loop_line2(image, W, H, x0, y0, x, y, 0, 0, 0);
    }
    lodepng::encode("star4.png", image, W, H, LCT_RGB, 8);*/

    /*for (int i = 0; i < 13; ++i) {
        double alpha = 2.0 * M_PI * i /  13;
        unsigned int x = (unsigned int)(x0 + 95 * cos(alpha));
        unsigned int y = (unsigned int)(y0 + 95 * sin(alpha));
        loop_line3(image, W, H, x0, y0, x, y, 0, 0, 0);
    }
    lodepng::encode("star5.png", image, W, H, LCT_RGB, 8);*/

    /*for (int i = 0; i < 13; ++i) {
        double alpha = 2.0 * M_PI * i /  13;
        unsigned int x = (unsigned int)(x0 + 95 * cos(alpha));
        unsigned int y = (unsigned int)(y0 + 95 * sin(alpha));
        loop_line4(image, W, H, x0, y0, x, y, 0, 0, 0);
    }
    lodepng::encode("star6.png", image, W, H, LCT_RGB, 8);*/

    for (int i = 0; i < 13; ++i) {
        double alpha = 2.0 * M_PI * i /  13;
        unsigned int x = (unsigned int)(x0 + 95 * cos(alpha));
        unsigned int y = (unsigned int)(y0 + 95 * sin(alpha));
        loop_line4(image, W, H, x0, y0, x, y, 0, 0, 0);
    }
    lodepng::encode("star7.png", image, W, H, LCT_RGB, 8);
}

int main() {
    const unsigned int W = 200;
    const unsigned int H = 200;

    task2(W, H);

    return 0;
}