# Gauss

Programs that build random linear systems and solve them with Gaussian elimination and Gauss-Jordan elimination. Arithmetic is exact: integer entries stay integers, and every other value is a reduced fraction `a/b`.

Each program has a fixed 3-variable version and a version for any size. For `n` variables, the augmented system is `n` rows of `n + 1` entries: the coefficients, then the right-hand side.

| Program | Input | Output |
| --- | --- | --- |
| `generator3` | none | random 3×4 system |
| `generator` | integer `n` > 1 | random `n`×(`n`+1) system |
| `solver3` | 3×4 integers | the three integer unknowns |
| `solver` | `n`×(`n`+1) integers | the `n` integer unknowns |
| `gauss3` | 3×4 rationals | row-echelon form |
| `gauss` | `n`×(`n`+1) rationals | row-echelon form |
| `gauss_jordan3` | 3×4 rationals | reduced row-echelon form |
| `gauss_jordan` | `n`×(`n`+1) rationals | reduced row-echelon form |

## Build

The programs are standalone C++17 files:

```bash
g++ -std=c++17 -Wall -Wextra -o generator3 generator3.cpp
g++ -std=c++17 -Wall -Wextra -o generator generator.cpp
g++ -std=c++17 -Wall -Wextra -o solver3 solver3.cpp
g++ -std=c++17 -Wall -Wextra -o solver solver.cpp
g++ -std=c++17 -Wall -Wextra -o gauss3 gauss3.cpp
g++ -std=c++17 -Wall -Wextra -o gauss gauss.cpp
g++ -std=c++17 -Wall -Wextra -o gauss_jordan3 gauss_jordan3.cpp
g++ -std=c++17 -Wall -Wextra -o gauss_jordan gauss_jordan.cpp
```

## Generators

`generator3` prints a system whose hidden solution is three random integers. `generator` does the same for a chosen size. Pass that size as an argument, or type it on standard input:

```bash
./generator 4
echo 4 | ./generator
```

Coefficients are random integers from −9 to 9. The coefficient matrix is invertible, so the hidden integer solution is unique. The solution itself is not printed.

```text
2 -1 -9 44
3 0 -8 28
-7 5 0 -7
```

means `2x − y − 9z = 44`, and likewise for the other rows.

## Solvers

`solver3` and `solver` read a system of integers and print the unique integer solution on one line.

```bash
./generator3 | ./solver3
./generator 4 | ./solver
```

The system in the previous section has the solution

```text
-4 -7 -5
```

A system with no unique integer solution is reported on standard error, and the program exits with status 1.

## Gaussian elimination

`gauss3` and `gauss` print the row-echelon form. A zero pivot is swapped with a lower nonzero entry, the pivot is scaled to 1, and the entries below it are cleared.

```bash
./generator3 | ./gauss3
```

```text
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

```text
1 1/2 -1/2 4
0 1 1 2
0 0 1 -1
```

## Gauss-Jordan elimination

`gauss_jordan3` and `gauss_jordan` continue from that form and clear the entries above each pivot as well. The result is the reduced row-echelon form. When the solution is unique, that form is the identity matrix followed by the solution.

```bash
./generator 4 | ./gauss_jordan
```

```text
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

```text
1 0 0 2
0 1 0 3
0 0 1 -1
```

## Input

The variable-size readers (`solver`, `gauss`, and `gauss_jordan`) take one equation per line and stop at end of input or at a blank line. Every row must contain `n + 1` entries.

`gauss`, `gauss3`, `gauss_jordan`, and `gauss_jordan3` accept integers and fractions:

```text
1/2 1/3 1/6 1
0 1/4 1/2 1
0 0 1/5 1
```

Fractions are reduced, and the sign is written on the numerator (`-1/2`, not `1/-2`). An integer is printed without a denominator.
