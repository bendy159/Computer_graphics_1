#include <iostream>
#include "lodepng.h"

using namespace std;

void black(vector<unsigned char>& image, unsigned int W, unsigned int H) {      //Чёрный фон
    for(unsigned int i = 0; i < H * W; ++i) {
        image[i] = 0;
    }

    lodepng::encode("black.png", image, W, H, LCT_GREY, 8);
}

void white(vector<unsigned char>& image, unsigned int W, unsigned int H) {      //Белый фон
    for(unsigned int i = 0; i < H * W; ++i) {
        image[i] = 255;
    }

    lodepng::encode("white.png", image, W, H, LCT_GREY, 8);
}

void red(vector<unsigned char>& image, unsigned int W, unsigned int H) {        //Красный фон
    for(unsigned int i = 0; i < H * W * 3; i += 3) {
        image[i] = 255;
        image[i + 1] = 0;
        image[i + 2] = 0;
    }

    lodepng::encode("red.png", image, W, H, LCT_RGB, 8);
}

void grad(vector<unsigned char>& image, unsigned int W, unsigned int H) {       //Что-то типа градиента
    for(unsigned int i = 0; i < H; i += 3) {
        for(unsigned int j = 0; j < W; ++j) {
            image[3 * (i * W + j)] = (i + j) % 256;
            image[3 * (i * W + j) + 1] = (i + j) % 256;
            image[3 * (i * W + j) + 2] = (i + j) % 256;
        }
    }

    lodepng::encode("grad.png", image, W, H, LCT_RGB, 8);       
}

int main() {
    const unsigned int W = 200;     //Задание границ
    const unsigned int H = 200;
    vector<unsigned char> image(H * W * 3);     //Вектор чаров, так как 0 - 255

    red(image, W, H);

    return 0;
}