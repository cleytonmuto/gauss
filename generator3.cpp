// Builds a 3x4 augmented matrix for a linear system whose solution is
// three random integers (x, y, z).
//
// Output (the augmented system only):
//   a11 a12 a13 b1
//   a21 a22 a23 b2
//   a31 a32 a33 b3
//
// Each right-hand side is bi = ai1*x + ai2*y + ai3*z. The coefficient
// matrix is kept invertible, so (x, y, z) is the unique solution.

#include <iostream>
#include <random>

namespace {

long determinant(const int a[3][4]) {
    const long a00 = a[0][0];
    const long a01 = a[0][1];
    const long a02 = a[0][2];
    const long a10 = a[1][0];
    const long a11 = a[1][1];
    const long a12 = a[1][2];
    const long a20 = a[2][0];
    const long a21 = a[2][1];
    const long a22 = a[2][2];

    return a00 * (a11 * a22 - a12 * a21) -
           a01 * (a10 * a22 - a12 * a20) +
           a02 * (a10 * a21 - a11 * a20);
}

}  // namespace

int main() {
    std::mt19937 rng{std::random_device{}()};
    std::uniform_int_distribution<int> value(-9, 9);

    const int x = value(rng);
    const int y = value(rng);
    const int z = value(rng);

    int augmented[3][4]{};
    do {
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                augmented[row][col] = value(rng);
            }
        }
    } while (determinant(augmented) == 0);

    for (int row = 0; row < 3; ++row) {
        augmented[row][3] = augmented[row][0] * x +
                            augmented[row][1] * y +
                            augmented[row][2] * z;
    }

    // std::cout << x << ' ' << y << ' ' << z << '\n';
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 4; ++col) {
            if (col > 0) {
                std::cout << ' ';
            }
            std::cout << augmented[row][col];
        }
        std::cout << '\n';
    }
}
