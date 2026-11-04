// Reads a 3x4 augmented linear system from stdin and prints its solution.
//
// Input:
//   a11 a12 a13 b1
//   a21 a22 a23 b2
//   a31 a32 a33 b3
//
// Output:
//   x y z
//
// Solves by Gaussian elimination with partial pivoting. The generator
// builds systems whose solution is three integers, so the result is
// checked against the original equations before it is printed.

#include <cmath>
#include <iostream>
#include <limits>

namespace {

constexpr int kRows = 3;
constexpr int kCols = 4;
constexpr double kPivotEpsilon = 1e-9;

bool read_system(long long original[kRows][kCols], double augmented[kRows][kCols]) {
    for (int row = 0; row < kRows; ++row) {
        for (int col = 0; col < kCols; ++col) {
            if (!(std::cin >> original[row][col])) {
                return false;
            }
            augmented[row][col] = static_cast<double>(original[row][col]);
        }
    }
    return true;
}

bool eliminate(double augmented[kRows][kCols]) {
    for (int col = 0; col < kRows; ++col) {
        int pivot = col;
        for (int row = col + 1; row < kRows; ++row) {
            if (std::fabs(augmented[row][col]) > std::fabs(augmented[pivot][col])) {
                pivot = row;
            }
        }
        if (std::fabs(augmented[pivot][col]) < kPivotEpsilon) {
            return false;
        }
        if (pivot != col) {
            for (int j = col; j < kCols; ++j) {
                const double saved = augmented[col][j];
                augmented[col][j] = augmented[pivot][j];
                augmented[pivot][j] = saved;
            }
        }

        for (int row = col + 1; row < kRows; ++row) {
            const double factor = augmented[row][col] / augmented[col][col];
            augmented[row][col] = 0.0;
            for (int j = col + 1; j < kCols; ++j) {
                augmented[row][j] -= factor * augmented[col][j];
            }
        }
    }
    return true;
}

bool back_substitute(const double augmented[kRows][kCols], long long solution[kRows]) {
    double value[kRows]{};
    for (int row = kRows - 1; row >= 0; --row) {
        double sum = augmented[row][kRows];
        for (int col = row + 1; col < kRows; ++col) {
            sum -= augmented[row][col] * value[col];
        }
        if (std::fabs(augmented[row][row]) < kPivotEpsilon) {
            return false;
        }
        value[row] = sum / augmented[row][row];

        const double rounded = std::round(value[row]);
        if (std::fabs(value[row] - rounded) > 1e-6 ||
            !std::isfinite(rounded) ||
            rounded < static_cast<double>(std::numeric_limits<long long>::min()) ||
            rounded > static_cast<double>(std::numeric_limits<long long>::max())) {
            return false;
        }
        solution[row] = static_cast<long long>(rounded);
    }
    return true;
}

bool satisfies(const long long original[kRows][kCols], const long long solution[kRows]) {
    for (int row = 0; row < kRows; ++row) {
        const long long left = original[row][0] * solution[0] +
                               original[row][1] * solution[1] +
                               original[row][2] * solution[2];
        if (left != original[row][3]) {
            return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    long long original[kRows][kCols]{};
    double augmented[kRows][kCols]{};
    if (!read_system(original, augmented)) {
        std::cerr << "expected 3 rows of 4 numbers\n";
        return 1;
    }

    long long solution[kRows]{};
    if (!eliminate(augmented) || !back_substitute(augmented, solution) ||
        !satisfies(original, solution)) {
        std::cerr << "system has no unique integer solution\n";
        return 1;
    }

    std::cout << solution[0] << ' ' << solution[1] << ' ' << solution[2] << '\n';
}
