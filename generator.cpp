// Builds an n by (n + 1) augmented matrix whose solution is n random
// integers. The dimension n is an integer greater than 1, taken from
// the command line or, if no argument is given, from stdin.
//
//   generator 4
//   echo 4 | generator
//
// Output is the augmented system only, one equation per line:
//   a11 a12 ... a1n b1
//   ...
//   an1 an2 ... ann bn
//
// Each right-hand side is the coefficient row times the hidden
// solution. The coefficient matrix is kept invertible, so that
// solution is unique.

#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

constexpr long long kMod = 1000000007LL;
constexpr int kMinDimension = 2;

long long mod_pow(long long base, long long exponent) {
    long long result = 1;
    base %= kMod;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = result * base % kMod;
        }
        base = base * base % kMod;
        exponent /= 2;
    }
    return result;
}

bool invertible(const std::vector<std::vector<long long>>& augmented) {
    const int dimension = static_cast<int>(augmented.size());
    std::vector<std::vector<long long>> matrix(dimension, std::vector<long long>(dimension));
    for (int row = 0; row < dimension; ++row) {
        for (int col = 0; col < dimension; ++col) {
            long long value = augmented[row][col] % kMod;
            if (value < 0) {
                value += kMod;
            }
            matrix[row][col] = value;
        }
    }

    for (int col = 0; col < dimension; ++col) {
        int pivot = col;
        while (pivot < dimension && matrix[pivot][col] == 0) {
            ++pivot;
        }
        if (pivot == dimension) {
            return false;
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
        }

        const long long inverse = mod_pow(matrix[col][col], kMod - 2);
        for (int row = col + 1; row < dimension; ++row) {
            if (matrix[row][col] == 0) {
                continue;
            }
            const long long factor = matrix[row][col] * inverse % kMod;
            for (int j = col; j < dimension; ++j) {
                matrix[row][j] = (matrix[row][j] - factor * matrix[col][j]) % kMod;
                if (matrix[row][j] < 0) {
                    matrix[row][j] += kMod;
                }
            }
        }
    }
    return true;
}

bool parse_dimension(const std::string& text, int& dimension) {
    if (text.empty()) {
        return false;
    }
    std::size_t end = 0;
    long long value = 0;
    try {
        value = std::stoll(text, &end);
    } catch (const std::exception&) {
        return false;
    }
    if (end != text.size() || value < kMinDimension ||
        value > std::numeric_limits<int>::max()) {
        return false;
    }
    dimension = static_cast<int>(value);
    return true;
}

bool read_dimension(int argc, char** argv, int& dimension) {
    if (argc == 2) {
        return parse_dimension(argv[1], dimension);
    }
    if (argc != 1) {
        return false;
    }

    std::string line;
    if (!std::getline(std::cin, line)) {
        return false;
    }
    const std::size_t last = line.find_last_not_of(" \t\r");
    if (last == std::string::npos) {
        return false;
    }
    return parse_dimension(line.substr(0, last + 1), dimension);
}

}  // namespace

int main(int argc, char** argv) {
    int dimension = 0;
    if (!read_dimension(argc, argv, dimension)) {
        std::cerr << "expected an integer greater than 1\n";
        return 1;
    }

    std::mt19937 rng{std::random_device{}()};
    std::uniform_int_distribution<int> value(-9, 9);

    std::vector<long long> solution(dimension);
    for (int index = 0; index < dimension; ++index) {
        solution[index] = value(rng);
    }

    std::vector<std::vector<long long>> augmented(
        dimension, std::vector<long long>(dimension + 1));
    do {
        for (int row = 0; row < dimension; ++row) {
            for (int col = 0; col < dimension; ++col) {
                augmented[row][col] = value(rng);
            }
        }
    } while (!invertible(augmented));

    for (int row = 0; row < dimension; ++row) {
        long long right_hand_side = 0;
        for (int col = 0; col < dimension; ++col) {
            right_hand_side += augmented[row][col] * solution[col];
        }
        augmented[row][dimension] = right_hand_side;
    }

    for (int row = 0; row < dimension; ++row) {
        for (int col = 0; col <= dimension; ++col) {
            if (col > 0) {
                std::cout << ' ';
            }
            std::cout << augmented[row][col];
        }
        std::cout << '\n';
    }
}
