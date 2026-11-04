// Reads an augmented linear system and prints its reduced row-echelon form.
//
// The number of rows is the number of variables. A system with n
// variables is n rows of n + 1 entries, coefficients followed by the
// right-hand side. Entries may be integers or fractions a/b. Input
// ends at EOF or at a blank line.
//
// Gauss-Jordan elimination swaps a zero pivot with a lower nonzero
// entry, scales the pivot to 1, and clears every other entry in that
// column. Integers are printed as integers. Every other value is a
// reduced fraction a/b.

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
    const __int128 divisor = gcd128(value.num, value.den);
    value.num /= divisor;
    value.den /= divisor;
}

Rational operator+(Rational left, Rational right) {
    const __int128 common = gcd128(left.den, right.den);
    const __int128 left_den = left.den / common;
    const __int128 right_den = right.den / common;
    Rational result{left.num * right_den + right.num * left_den, common * left_den * right_den};
    normalize(result);
    return result;
}

Rational operator-(Rational left, Rational right) {
    right.num = -right.num;
    return left + right;
}

Rational operator*(Rational left, Rational right) {
    const __int128 left_cancel = gcd128(left.num, right.den);
    const __int128 right_cancel = gcd128(right.num, left.den);
    Rational result{(left.num / left_cancel) * (right.num / right_cancel),
                    (left.den / right_cancel) * (right.den / left_cancel)};
    normalize(result);
    return result;
}

Rational operator/(Rational left, Rational right) {
    return left * Rational{right.den, right.num};
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

bool read_system(std::vector<std::vector<Rational>>& rows) {
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(std::cin, line)) {
        ++line_number;
        std::istringstream input(line);
        std::vector<Rational> row;
        std::string token;
        while (input >> token) {
            Rational value;
            if (!parse_rational(token, value)) {
                std::cerr << "row " << line_number << " contains an invalid rational number\n";
                return false;
            }
            row.push_back(value);
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

bool widths_match(const std::vector<std::vector<Rational>>& rows) {
    if (rows.empty()) {
        std::cerr << "expected one or more rows\n";
        return false;
    }
    const std::size_t width = rows.size() + 1;
    for (std::size_t row = 0; row < rows.size(); ++row) {
        if (rows[row].size() != width) {
            std::cerr << "row " << (row + 1) << " has " << rows[row].size()
                      << " numbers; expected " << width << '\n';
            return false;
        }
    }
    return true;
}

bool is_zero_row(const std::vector<Rational>& row) {
    for (const Rational& value : row) {
        if (value.num != 0) {
            return false;
        }
    }
    return true;
}

void eliminate(std::vector<std::vector<Rational>>& augmented) {
    const std::size_t rows = augmented.size();
    const std::size_t columns = augmented[0].size();
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < columns && pivot_row < rows; ++col) {
        std::size_t selected = rows;
        for (std::size_t row = pivot_row; row < rows; ++row) {
            if (augmented[row][col].num != 0) {
                selected = row;
                break;
            }
        }
        if (selected == rows) {
            continue;
        }
        if (selected != pivot_row) {
            std::swap(augmented[pivot_row], augmented[selected]);
        }

        const Rational pivot = augmented[pivot_row][col];
        for (std::size_t col_index = col; col_index < columns; ++col_index) {
            augmented[pivot_row][col_index] = augmented[pivot_row][col_index] / pivot;
        }
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const Rational factor = augmented[row][col];
            if (factor.num == 0) {
                continue;
            }
            for (std::size_t col_index = col; col_index < columns; ++col_index) {
                augmented[row][col_index] =
                    augmented[row][col_index] - factor * augmented[pivot_row][col_index];
            }
        }
        ++pivot_row;
    }

    for (std::size_t pass = 0; pass < rows; ++pass) {
        for (std::size_t row = 0; row + 1 < rows; ++row) {
            if (is_zero_row(augmented[row]) && !is_zero_row(augmented[row + 1])) {
                std::swap(augmented[row], augmented[row + 1]);
            }
        }
    }
}

void print_rational(Rational value) {
    if (value.den == 1) {
        std::cout << decimal(value.num);
        return;
    }
    std::cout << decimal(value.num) << '/' << decimal(value.den);
}

void print_system(const std::vector<std::vector<Rational>>& augmented) {
    for (const std::vector<Rational>& row : augmented) {
        for (std::size_t col = 0; col < row.size(); ++col) {
            if (col > 0) {
                std::cout << ' ';
            }
            print_rational(row[col]);
        }
        std::cout << '\n';
    }
}

}  // namespace

int main() {
    std::vector<std::vector<Rational>> augmented;
    if (!read_system(augmented) || !widths_match(augmented)) {
        return 1;
    }

    eliminate(augmented);
    print_system(augmented);
}
