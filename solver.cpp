// Reads an augmented linear system from stdin and prints its solution.
//
// The number of rows is the number of variables. Row i holds the
// coefficients of that equation and the right-hand side, so a system
// with n variables is n rows of n + 1 integers. Input ends at EOF or
// at a blank line. The solution is printed on one line:
//
//   x1 x2 ... xn
//
// Solves by exact Gaussian elimination with partial pivoting, and
// prints the result only when every variable is an integer.

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Rational {
    __int128 num = 0;
    __int128 den = 1;
};

__int128 gcd128(__int128 a, __int128 b) {
    if (a < 0) {
        a = -a;
    }
    if (b < 0) {
        b = -b;
    }
    while (b != 0) {
        const __int128 remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
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
    const __int128 divisor = gcd128(value.num, value.den);
    value.num /= divisor;
    value.den /= divisor;
}

Rational make_rational(long long value) {
    return Rational{value, 1};
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

double magnitude(Rational value) {
    return std::fabs(static_cast<double>(value.num) / static_cast<double>(value.den));
}

std::string decimal(__int128 value) {
    if (value == 0) {
        return "0";
    }
    const bool negative = value < 0;
    unsigned __int128 magnitude = negative
                                      ? static_cast<unsigned __int128>(-(value + 1)) + 1
                                      : static_cast<unsigned __int128>(value);
    std::string digits;
    while (magnitude > 0) {
        digits.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) {
        digits.push_back('-');
    }
    return std::string(digits.rbegin(), digits.rend());
}

bool read_system(std::vector<std::vector<long long>>& rows) {
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(std::cin, line)) {
        ++line_number;
        std::istringstream input(line);
        std::vector<long long> row;
        long long value = 0;
        while (input >> value) {
            row.push_back(value);
        }
        if (input.fail() && !input.eof()) {
            std::cerr << "row " << line_number << " contains a non-integer value\n";
            return false;
        }
        if (row.empty()) {
            if (!rows.empty()) {
                break;
            }
            continue;
        }
        rows.push_back(std::move(row));
    }
    return true;
}

bool widths_match(const std::vector<std::vector<long long>>& rows) {
    const std::size_t variables = rows.size();
    if (variables == 0) {
        std::cerr << "expected one or more rows\n";
        return false;
    }
    for (std::size_t row = 0; row < variables; ++row) {
        if (rows[row].size() != variables + 1) {
            std::cerr << "row " << (row + 1) << " has " << rows[row].size()
                      << " numbers; expected " << (variables + 1) << '\n';
            return false;
        }
    }
    return true;
}

bool eliminate(std::vector<std::vector<Rational>>& augmented) {
    const std::size_t variables = augmented.size();
    for (std::size_t col = 0; col < variables; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < variables; ++row) {
            if (magnitude(augmented[row][col]) > magnitude(augmented[pivot][col])) {
                pivot = row;
            }
        }
        if (augmented[pivot][col].num == 0) {
            return false;
        }
        if (pivot != col) {
            std::swap(augmented[col], augmented[pivot]);
        }

        for (std::size_t row = col + 1; row < variables; ++row) {
            const Rational factor = augmented[row][col] / augmented[col][col];
            for (std::size_t j = col; j < variables + 1; ++j) {
                augmented[row][j] = augmented[row][j] - factor * augmented[col][j];
            }
        }
    }
    return true;
}

bool back_substitute(const std::vector<std::vector<Rational>>& augmented,
                     std::vector<Rational>& solution) {
    const std::size_t variables = augmented.size();
    solution.assign(variables, Rational{});
    for (std::size_t row = variables; row-- > 0;) {
        if (augmented[row][row].num == 0) {
            return false;
        }
        Rational sum = augmented[row][variables];
        for (std::size_t col = row + 1; col < variables; ++col) {
            sum = sum - augmented[row][col] * solution[col];
        }
        solution[row] = sum / augmented[row][row];
        if (solution[row].den != 1) {
            return false;
        }
    }
    return true;
}

bool satisfies(const std::vector<std::vector<long long>>& original,
               const std::vector<Rational>& solution) {
    const std::size_t variables = original.size();
    for (std::size_t row = 0; row < variables; ++row) {
        __int128 left = 0;
        for (std::size_t col = 0; col < variables; ++col) {
            left += static_cast<__int128>(original[row][col]) * solution[col].num;
        }
        if (left != original[row][variables]) {
            return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    std::vector<std::vector<long long>> original;
    if (!read_system(original) || !widths_match(original)) {
        return 1;
    }

    const std::size_t variables = original.size();
    std::vector<std::vector<Rational>> augmented(variables, std::vector<Rational>(variables + 1));
    for (std::size_t row = 0; row < variables; ++row) {
        for (std::size_t col = 0; col < variables + 1; ++col) {
            augmented[row][col] = make_rational(original[row][col]);
        }
    }

    std::vector<Rational> solution;
    if (!eliminate(augmented) || !back_substitute(augmented, solution) ||
        !satisfies(original, solution)) {
        std::cerr << "system has no unique integer solution\n";
        return 1;
    }

    for (std::size_t i = 0; i < variables; ++i) {
        if (i > 0) {
            std::cout << ' ';
        }
        std::cout << decimal(solution[i].num);
    }
    std::cout << '\n';
}
