// Implementation for Problem 092 - Resistive Network
// This header defines the resistive_network class expected by the OJ.

#ifndef SRC_HPP
#define SRC_HPP

// Bring in standard headers first so the provided fraction.hpp compiles cleanly
#include <bits/stdc++.h>

#include fraction.hpp

// We ignore the provided matrix helper and implement directly with simple vectors.
#define IGNORE_MATRIX 1

class resistive_network {
private:
    int n; // number of nodes (interfaces)
    int m; // number of connections

    // Edge lists (1-based indices as provided by input)
    std::vector<int> u; // from
    std::vector<int> v; // to
    std::vector<fraction> r; // resistances
    std::vector<fraction> g; // conductances = 1 / r

    // Laplacian matrix L of size n x n with rational entries
    std::vector<std::vector<fraction>> L;

    // Build Laplacian from edges
    void build_laplacian() {
        L.assign(n, std::vector<fraction>(n, fraction(0)));
        for (int k = 0; k < m; ++k) {
            int a = u[k] - 1;
            int b = v[k] - 1;
            const fraction &c = g[k];
            L[a][a] = L[a][a] + c;
            L[b][b] = L[b][b] + c;
            L[a][b] = L[a][b] - c;
            L[b][a] = L[b][a] - c;
        }
    }

    // Solve A x = b using Gaussian elimination over fractions.
    // A is square (size sz x sz), b is length sz. Returns x.
    static std::vector<fraction> solve_linear(std::vector<std::vector<fraction>> A,
                                              std::vector<fraction> b) {
        int sz = (int)A.size();
        // Standard Gaussian elimination with partial pivoting
        for (int k = 0; k < sz; ++k) {
            int pivot = -1;
            for (int i = k; i < sz; ++i) {
                if (!(A[i][k] == fraction(0))) { pivot = i; break; }
            }
            if (pivot == -1) {
                // Singular column; continue (system guaranteed solvable in problem constraints)
                continue;
            }
            if (pivot != k) {
                std::swap(A[pivot], A[k]);
                std::swap(b[pivot], b[k]);
            }
            fraction piv = A[k][k];
            for (int i = k + 1; i < sz; ++i) {
                if (A[i][k] == fraction(0)) continue;
                fraction factor = A[i][k] / piv;
                for (int j = k; j < sz; ++j) {
                    A[i][j] = A[i][j] - factor * A[k][j];
                }
                b[i] = b[i] - factor * b[k];
            }
        }

        // Back substitution
        std::vector<fraction> x(sz, fraction(0));
        for (int i = sz - 1; i >= 0; --i) {
            // Find leading coefficient position
            int lead = -1;
            for (int j = 0; j < sz; ++j) {
                if (!(A[i][j] == fraction(0))) { lead = j; break; }
            }
            if (lead == -1) continue; // zero row
            fraction sum = fraction(0);
            for (int j = lead + 1; j < sz; ++j) {
                sum = sum + A[i][j] * x[j];
            }
            x[lead] = (b[i] - sum) / A[i][lead];
        }
        return x;
    }

    // Build reduced matrix A by removing ground index (0-based), and corresponding b vector
    void build_reduced_system(int ground,
                              std::vector<std::vector<fraction>> &A,
                              const std::vector<fraction> &b_full,
                              std::vector<fraction> &b_red) const {
        int sz = n - 1;
        A.assign(sz, std::vector<fraction>(sz, fraction(0)));
        b_red.clear(); b_red.reserve(sz);
        // Map original indices to reduced
        auto map_index = [&](int idx) -> int {
            return (idx < ground) ? idx : idx - 1;
        };
        for (int i = 0; i < n; ++i) {
            if (i == ground) continue;
            int ri = map_index(i);
            b_red.push_back(b_full[i]);
            for (int j = 0; j < n; ++j) {
                if (j == ground) continue;
                int rj = map_index(j);
                A[ri][rj] = L[i][j];
            }
        }
    }

public:
    // Constructor: set up network with interface_size_ nodes and connection_size_ edges
    // from[k], to[k] are 1-based endpoints for edge k, resistance[k] > 0
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        n = interface_size_;
        m = connection_size_;
        u.resize(m); v.resize(m); r.resize(m); g.resize(m);
        for (int i = 0; i < m; ++i) {
            u[i] = from[i];
            v[i] = to[i];
            r[i] = resistance[i];
            g[i] = fraction(1) / r[i];
        }
        build_laplacian();
    }

    ~resistive_network() = default;

    // Equivalent resistance between two nodes (1-based). If same node, zero.
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1;
        int b = interface_id2;
        if (a == b) return fraction(0);

        int ground = b - 1; // ground at node b
        // Current injections: +1 at a, -1 at b
        std::vector<fraction> I(n, fraction(0));
        I[a - 1] = fraction(1);
        I[b - 1] = I[b - 1] - fraction(1);

        // Build reduced system excluding ground node b
        std::vector<std::vector<fraction>> A;
        std::vector<fraction> bvec;
        build_reduced_system(ground, A, I, bvec);
        // Solve for voltages (relative to ground)
        std::vector<fraction> Ured = solve_linear(A, bvec);

        // Extract voltage at node a (if a==ground handled earlier)
        int idx_a = (a - 1 < ground) ? (a - 1) : (a - 2);
        return Ured[idx_a];
    }

    // Given node currents I (length n), with u_n = 0, return voltage at node id (1-based, id < n)
    fraction get_voltage(int id, fraction current[]) {
        // Build current vector
        std::vector<fraction> I(n);
        for (int i = 0; i < n; ++i) I[i] = current[i];
        // Ground is node n-1 (0-based)
        int ground = n - 1;
        std::vector<std::vector<fraction>> A;
        std::vector<fraction> bvec;
        build_reduced_system(ground, A, I, bvec);
        std::vector<fraction> Ured = solve_linear(A, bvec);
        // id < n by contract
        return Ured[id - 1];
    }

    // Given node voltages U (length n), return total dissipated power
    fraction get_power(fraction voltage[]) {
        fraction total(0);
        for (int k = 0; k < m; ++k) {
            fraction du = voltage[u[k] - 1] - voltage[v[k] - 1];
            // Power on resistor k: (du^2) / r_k
            total = total + (du * du) / r[k];
        }
        return total;
    }
};

#endif // SRC_HPP
