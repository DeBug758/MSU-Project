#include <iostream>
#include <vector>

using namespace std;

template<typename T>
class Matrix {
private:
    vector<vector<T>> data;
    int rows;
    int columns;

public:
    Matrix(int rows, int columns, const vector<vector<T>>& initData) {
        this->rows = rows;
        this->columns = columns;
        data = initData;
    }


    vector<T> relax(const vector<T>& b, int maxiter, double eps, double omega) {
        vector<T> x(rows, 0);
        vector<T> New_x(rows);
        int iter = 0;
        T error = eps + 1.0;

        while (iter < maxiter && error > eps) {
            error = 0;
            for (int i = 0; i < rows; i++) {
                New_x[i] = function(b, x, i, omega);
                T diff = abs(New_x[i] - x[i]);
                if (diff > error) {
                    error = diff;
                }
            }
            x = New_x;
            iter++;
        }

        if (iter >= maxiter) {
            cout << "Failed to achieve the required accuracy for the specified number of iterations." << endl;
        }
        else {
            cout << "The relaxation method converges after " << iter << " iteretions" << endl;
        }

        return x;
    }

    T function(const vector<T>& b, const vector<T>& x, int i, double omega) {
        T sum = 0;
        for (int j = 0; j < columns; j++) {
            if (j != i) {
                sum += data[i][j] * x[j];
            }
        }

        T newX = (b[i] - sum) / data[i][i];
        newX = omega * newX + (1.0 - omega) * x[i];
        return newX;
    }
};

int main() {
    vector<vector<double>> data{
        {5, 2, -1},
        {-4, 7, 3},
        {2, -2, 4}
    };
    vector<double> B{ 12, 24, 9};
    Matrix<double> A(3, 3, data);
    vector<double> result = A.relax(B, 100, 0.1, 0.25);

    cout << "System solution:" << endl;
    for (int i = 0; i < result.size(); i++) {
        cout << "x" << i << "= " << result[i] << endl;
    }

    return 0;
}