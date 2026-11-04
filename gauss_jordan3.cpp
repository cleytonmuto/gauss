// Reads a 3x4 augmented linear system and prints its reduced row-echelon form.
//
// Input is three equations in three variables. Each row is
//   a b c d
// for a*x + b*y + c*z = d. Entries may be integers or fractions a/b.
//
// Gauss-Jordan elimination swaps a zero pivot with a lower nonzero
// entry, scales the pivot to 1, and clears every other entry in that
// column. Integers are printed as integers. Every other value is a
// reduced fraction a/b.

#include <iostream>
#include <string>
#include <utility>

namespace {

constexpr int kRows = 3;
constexpr int kCols = 4;

struct Rational {
    long long num = 0;
    long long den = 1;
};

long long gcd(long long a, long long b) {
    if (a < 0) {
        a = -a;
    }
    if (b < 0) {
        b = -b;
    }
    while (b != 0) {
        const long long remainder = a % b;
        a = b;
        b = remainder;
    }
    return a == 0 ? 1 : a;
}

void normalize(Rational& value) {
    if (value.den < 0) {
        value.num = -value.num;
        value.den = -value.den;
    }
    if (value.num == 0) {
        value.den = 1;
        return;
    }
    const long long divisor = gcd(value.num, value.den);
    value.num /= divisor;
    value.den /= divisor;
}

Rational operator+(Rational left, Rational right) {
    Rational result{left.num * right.den + right.num * left.den, left.den * right.den};
    normalize(result);
    return result;
}

Rational operator-(Rational left, Rational right) {
    right.num = -right.num;
    return left + right;
}

Rational operator*(Rational left, Rational right) {
    Rational result{left.num * right.num, left.den * right.den};
    normalize(result);
    return result;
}

Rational operator/(Rational left, Rational right) {
    Rational result{left.num * right.den, left.den * right.num};
    normalize(result);
    return result;
}

bool parse_rational(const std::string& token, Rational& value) {
    try {
        const std::size_t slash = token.find('/');
        if (slash == std::string::npos) {
            std::size_t end = 0;
            const long long num = std::stoll(token, &end);
            if (end != token.size()) {
                return false;
            }
            value = Rational{num, 1};
            return true;
        }

        std::size_t end_num = 0;
        std::size_t end_den = 0;
        const std::string numerator = token.substr(0, slash);
        const std::string denominator = token.substr(slash + 1);
        const long long num = std::stoll(numerator, &end_num);
        const long long den = std::stoll(denominator, &end_den);
        if (end_num != numerator.size() || end_den != denominator.size() || den == 0) {
            return false;
        }
        value = Rational{num, den};
        normalize(value);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool read_system(Rational augmented[kRows][kCols]) {
    for (int row = 0; row < kRows; ++row) {
        for (int col = 0; col < kCols; ++col) {
            std::string token;
            if (!(std::cin >> token) || !parse_rational(token, augmented[row][col])) {
                return false;
            }
        }
    }
    return true;
}

bool is_zero_row(const Rational row[kCols]) {
    for (int col = 0; col < kCols; ++col) {
        if (row[col].num != 0) {
            return false;
        }
    }
    return true;
}

void swap_rows(Rational augmented[kRows][kCols], int left, int right) {
    for (int col = 0; col < kCols; ++col) {
        std::swap(augmented[left][col], augmented[right][col]);
    }
}

void eliminate(Rational augmented[kRows][kCols]) {
    int pivot_row = 0;
    for (int col = 0; col < kCols && pivot_row < kRows; ++col) {
        int selected = -1;
        for (int row = pivot_row; row < kRows; ++row) {
            if (augmented[row][col].num != 0) {
                selected = row;
                break;
            }
        }
        if (selected < 0) {
            continue;
        }
        if (selected != pivot_row) {
            swap_rows(augmented, pivot_row, selected);
        }

        const Rational pivot = augmented[pivot_row][col];
        for (int col_index = col; col_index < kCols; ++col_index) {
            augmented[pivot_row][col_index] = augmented[pivot_row][col_index] / pivot;
        }
        for (int row = 0; row < kRows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const Rational factor = augmented[row][col];
            if (factor.num == 0) {
                continue;
            }
            for (int col_index = col; col_index < kCols; ++col_index) {
                augmented[row][col_index] =
                    augmented[row][col_index] - factor * augmented[pivot_row][col_index];
            }
        }
        ++pivot_row;
    }

    for (int pass = 0; pass < kRows; ++pass) {
        for (int row = 0; row < kRows - 1; ++row) {
            if (is_zero_row(augmented[row]) && !is_zero_row(augmented[row + 1])) {
                swap_rows(augmented, row, row + 1);
            }
        }
    }
}

void print_rational(Rational value) {
    if (value.den == 1) {
        std::cout << value.num;
        return;
    }
    std::cout << value.num << '/' << value.den;
}

void print_system(const Rational augmented[kRows][kCols]) {
    for (int row = 0; row < kRows; ++row) {
        for (int col = 0; col < kCols; ++col) {
            if (col > 0) {
                std::cout << ' ';
            }
            print_rational(augmented[row][col]);
        }
        std::cout << '\n';
    }
}

}  // namespace

int main() {
    Rational augmented[kRows][kCols]{};
    if (!read_system(augmented)) {
        std::cerr << "expected 3 rows of 4 rational numbers\n";
        return 1;
    }

    eliminate(augmented);
    print_system(augmented);
}
