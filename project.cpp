#include <iostream>
#include <vector>
#include <cmath>

template <typename T>
class Matrix {
private:
    unsigned int rows;
    unsigned int cols;
    std::vector<std::vector<T> > data;

public:
	Matrix(){
	}
	 
    Matrix(unsigned int nmRows, unsigned int nmCols, T value = 0) : rows(nmRows), cols(nmCols) {
        data.resize(rows, std::vector<T>(cols, value));
    }

    Matrix(unsigned int nmRows, unsigned int nmCols, const std::vector<std::vector<T> >& values)
        : rows(nmRows), cols(nmCols), data(values) {}

	Matrix(unsigned int r){
        data.resize(r, std::vector<T>(r, 0));
        for (unsigned int i = 0; i < r; ++i) {
            data[i][i] = 1;
        }
        this->rows = r;
        this->cols = r;
    }

int getRows() const{
    return this->rows;
}

void setRows(int value) {
    data.resize(value, std::vector<T>(cols));
    this->rows = value;
}

int getCols() const{
    return this->cols;
}

void setCols(int value){
    data.resize(rows, std::vector<T>(value));
    this->cols = value;
}

void getSize() const {
    std::cout << "Rows: " << rows << std::endl;
    std::cout << "Cols: " << cols << std::endl;
} 

void setSize(int rows, int cols){
    data.resize(rows, std::vector<T>(cols));
    this->rows = rows;
    this->cols = cols;
}

template<typename U>
Matrix<U> getRow(int value) const{
    if(value >= rows){
        std::cout << "Invalid number" << std::endl;
        return Matrix<U>();
    }
    Matrix<U> result(1, cols);
    for(int i = 0; i < cols; ++i){
        result.data[0][i] = data[value][i];
    }

    return result;
}

void setRow(int value, const std::vector<T>& values){
    if(value >= rows){
        std::cout << "Invalid number"<<std::endl;
    }
    else{
        this->data[value] = values;
    }
}

template<typename U>
Matrix<U> getCol(int value) const{
    if(value >= cols){
        std::cout << "Invalid number" << std::endl;
        return Matrix<U>();
    }
    Matrix<U> result(rows, 1);
    for(int i = 0; i < rows; ++i){
        result.data[i][0] = data[i][value];
    }
    
    return result;
}

void setCol(int value, const std::vector<T>& values){
    if(value >= cols){
        std::cout << "Invalid number"<<std::endl;
    }
    else{
        for(int i = 0; i < rows; ++i){
            data[i][value] = values[i];
        }
    }
}

void swapCols(int first, int second){
    if(first >= cols || second >= cols){
        std::cout<<"Invalid number"<<std::endl;
    }
    else{
        T tmp;
        for(int i = 0; i < rows; ++i){
            tmp = data[i][first];
            data[i][first] = data[i][second];
            data[i][second] = tmp;
        }
    }
}

void swapRows(int first, int second){
    if(first >= cols || second >= cols){
        std::cout<<"Invalid number"<<std::endl;
    }
    else{
        std::vector<T> tmp = data[first];
        data[first] = data[second];
        data[second] = tmp;
    }
}

void trans(){
    if(cols > rows){
        data.resize(cols, std::vector<T>(cols));
    }
    else{
        data.resize(rows, std::vector<T>(rows));
    }
    for(size_t i = 0; i < rows; ++i){
        for (size_t j = i; j < cols; j++)
        {
            data[j][i] = data[i][j];
        }
    }
    data.resize(cols, std::vector<T>(rows));
    int cnt = cols;
    cols = rows;
    rows = cnt;
}

void scan(){
    std::cout << "Enter number of Rows: ";
    std::cin >> this->rows;
    std::cout << "Enter number of Cols: ";
    std::cin >> this->cols;
    std::cout << "Enter your values:" << std::endl;
    for(int i = 0; i < this->rows; ++i){
        for(int j = 0; j < this->cols; ++j){
            std::cin >> this->data[i][j];
        }
    }
}

void setElement(unsigned int row, unsigned int col, T value) {
    if (row < rows && col < cols) {
         data[row][col] = value;
    }
}

T getElement(unsigned int row, unsigned int col) const {
    if (row < rows && col < cols) {
        return data[row][col];
    } 
    else {
        std::cout << "Error: Invalid indices." << std::endl;
        return T();  // Default value
    }
}

void printMatrix() const {
    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < cols; ++j) {
            std::cout << ' ' << data[i][j] << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template <typename U>
Matrix<U> operator+(const Matrix<U>& other) const {
        if (rows != other.rows || cols != other.cols) {
            std::cout << "Error: Matrices must have the same dimensions for addition." << std::endl;
            return Matrix<U>(rows, cols);
        }

    Matrix<U> result(rows, cols);
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] + other.getElement(i, j);
            }
        }
        return result;
    }

template <typename U>
Matrix<U> operator-(const Matrix<U>& other) const {
        if (rows != other.rows || cols != other.cols) {
            std::cout << "Error: Matrices must have the same dimensions for addition." << std::endl;
            return Matrix<U>(rows, cols);
        }

    Matrix<U> result(rows, cols);
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] - other.getElement(i, j);
            }
        }
        return result;
    }

template <typename U>
Matrix<U> operator*(const Matrix<U>& other) const {
        if (cols != other.rows) {
            std::cout << "Error: Invalid matrix dimensions for multiplication." << std::endl;
            return Matrix<U>(rows, other.cols);
        }

        Matrix<U> result(rows, other.cols);
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < other.cols; ++j) {
                for (unsigned int k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.getElement(k, j);
                }
            }
        }
        return result;
    }

template <typename U>
Matrix<U> operator*(U scalar) const {
 Matrix<U> result(rows, cols);
    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] * scalar;
        }
    }
    return result;
}

template <typename U>
Matrix<U> operator/(U scalar) const {
 Matrix<U> result(rows, cols);
    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] / scalar;
        }
    }
    return result;
}

T operator %(const Matrix <T>& other) {
        T scalar{};
        if (cols != other.cols  rows != other.rows)
            exit(0); 

        if ((cols != 1 && rows != 1)  (other.cols != 1 && other.rows != 1)) 
            exit(0);
        

        if (cols != 1) {
                for (int i{}; i < cols; i++) {
                    scalar += data[0][i]*other.data[0][i];
                }
        }
        
        if (rows != 1) {
            for (int i{}; i < rows; i++)
                scalar += data[i][0] * other.data[i][0];
        }
        scalar = sqrt(scalar);
        return scalar;
}

template <typename U>
void operator=(const Matrix<U>& other){
    rows = other.rows;
    cols = other.cols;
    data = other.data;
}

template <typename U>
Matrix<U> resetRow(unsigned int value){
    if(value >= rows && value != 0){
        std::cout << "Invalid number"<<std::endl;
        return Matrix<U>();
    }
    Matrix<U> result(rows - value, cols);
    for(int i = value; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
            result[i - value][j] = data[i][j];
        }
    }
    rows = value;
    return result;
}

//template <typename U>
Matrix<T> resetCol(unsigned int value){
    if(value >= cols && value != 0){
        std::cout << "Invalid number"<<std::endl;
        return Matrix<T>();
    }
    Matrix<T> result(rows, cols - value);
    for(int i = 0; i < rows; ++i){
        for(int j = value; j < cols; ++j){
            result.data[i][j - value] = data[i][j];
        }
    }
    cols = value;
    return result;
}

void appendMatrix(const Matrix<T>& other) {
    if (rows != other.rows) {
        std::cout << "Error: Matrices must have the same number of rows to append." << std::endl;
        return;
    }

    std::vector<std::vector<T> > newData(rows, std::vector<T>(cols + other.cols, 0));
    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < cols; ++j) {
            newData[i][j] = data[i][j];
        }
        for (unsigned int j = 0; j < other.cols; ++j) {
            newData[i][j + cols] = other.getElement(i, j);
        }
    }

    cols += other.cols;
    data = newData;
}

void appendMatrix(const std::vector<T>& other) {
    if(rows != other.size()){
      std::cout<<"EROR SIZE"<<std::endl;
      return;
    }
    std::vector< std::vector<T> > newData(rows, std::vector<T>(cols + 1, 0));
    for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                newData[i][j] = data[i][j];
            }
            newData[i][cols] = other[i];
        }
    data = newData;
    cols++;
  }

void prependMatrix(const Matrix<T>& other){
    if (rows != other.rows) {
        std::cout << "Error: Matrices must have the same number of rows to append." << std::endl;
        return;
    }
    std::vector<std::vector<T> > newData(rows, std::vector<T>(cols + other.cols, 0));
    for (unsigned int i = 0; i < rows; ++i) {
         for (unsigned int j = 0; j < other.cols; ++j) {
            newData[i][j] = other.getElement(i, j);
        }

        for (unsigned int j = 0; j < cols; ++j) {
            newData[i][j + other.cols] = data[i][j];
        }
       
    }
    cols += other.cols;
    data = newData;
}

void prependMatrix(const std::vector<T>& other){
    if (rows != other.rows) {
        std::cout<<"EROR SIZE"<<std::endl;
        return;
    }
    std::vector<std::vector<T> > newData(rows, std::vector<T>(cols + other.cols, 0));
    for (unsigned int i = 0; i < rows; ++i) {
        newData[i][0] = other[i];
        for (unsigned int j = 0; j < cols; ++j) {
            newData[i][j + other.cols] = data[i][j];
        }
       
    }
    cols += other.cols;
    data = newData;
}

Matrix<T> const gaussJordan() {
    Matrix<T> tmp(rows, cols, data);
    for (unsigned int i = 0; i < rows; ++i) {
        unsigned int maxRow = i;
        for (unsigned int k = i + 1; k < rows; ++k) {
            if (std::abs(tmp.data[k][i]) > std::abs(tmp.data[maxRow][i])) {
                maxRow = k;
            }
        }

       
        if (maxRow != i) {
            std::swap(tmp.data[i], tmp.data[maxRow]);
        }
        
        T factor = tmp.data[i][i];
        if (factor == 0) {
            std::cout << "Error: Matrix is singular, cannot convert to identity matrix." << std::endl;
            return *this;
        }

        for (unsigned int j = i; j < cols; ++j) {
            tmp.data[i][j] /= factor;
        }
        
        for (unsigned int k = 0; k < rows; ++k) {
            if (k != i) {
                factor = tmp.data[k][i];
                for (unsigned int j = i; j < cols; ++j) {
                    tmp.data[k][j] -= factor * tmp.data[i][j];
                }
            }
        }
    }
    return tmp;
}

const Matrix<T> inverse(){		
    if(cols != rows){
        std::cout<<"cols != rows"<<std::endl;
        return *this;
    }
    Matrix<T> tmp(rows, cols, data);
    tmp.appendMatrix(Matrix<T>(rows));
    tmp = tmp.gaussJordan();
    return tmp.resetCol(rows);
}

int cond1() const{
    return this->norm1() *  (this.inverse()).norm1();
}

int cond2() const{
    return this->norm2() *  (this.inverse()).norm2();
}

Matrix<T> LU() const {
        if (rows != cols) {
            std::cout << "Error: LU decomposition requires a square matrix." << std::endl;
            return Matrix<T>(rows, cols);
        }

        Matrix<T> lower(rows, cols);
        Matrix<T> upper(rows, cols);

        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                if (j < i) {
                    lower.setElement(j, i, 0);
                } else {
                    lower.setElement(j, i, data[j][i]);
                    for (unsigned int k = 0; k < i; ++k) {
                        lower.data[j][i] -= lower.data[j][k] * upper.data[k][i];
                    }
                }
            }

            for (unsigned int j = 0; j < cols; ++j) {
                if (j < i) {
                    upper.setElement(i, j, 0);
                } else if (j == i) {
                    upper.setElement(i, j, 1);
                } else {
                    upper.setElement(i, j, data[i][j] / lower.data[i][i]);
                    for (unsigned int k = 0; k < i; ++k) {
                        upper.data[i][j] -= (lower.data[i][k] * upper.data[k][j]) / lower.data[i][i];
                    }
                }
            }
        }

        Matrix<T> lu = lower * upper;
        std::cout << "Matrix L:" << std::endl;
        lower.printMatrix();
        std::cout << "Matrix U:" << std::endl;
        upper.printMatrix();
   
        return lu;
    }

double norm1() const{
    double max = 0.0;
    double sum;
    for (int i = 0; i < this->cols; ++i) {
        sum = 0.0;
        for(int j = 0; j < this->rows; ++j){
            sum += data[j][i];
        }
        if(sum > max){
            max = sum;
        }
    }
    return max;
}

double norm2() const{
    double sum = 0.0;
    for (int i = 0; i < this->rows; ++i) {
        for(int j = 0; j < this->cols; ++j){
            sum += data[i][j] * data[i][j]; 
        }
    }
    return std::sqrt(sum);
}

//TODO
/*void qrDecomposition(Matrix<T>& Q, Matrix<T>& R) const {
        
    }*/

  
void Gaus(){
    for (unsigned int i = 0; i < rows; ++i) {
        
        unsigned int maxRow = i;
        for (unsigned int k = i + 1; k < rows; ++k) {
            if (std::abs(data[k][i]) > std::abs(data[maxRow][i])) {
                maxRow = k;
            }
        }

        if (maxRow != i) {
            std::swap(data[i], data[maxRow]);
        }
    for (unsigned int i = 0; i < rows; ++i) {
      
      for (unsigned int k {i+1}; k < rows; ++k) {
      T factor = data[k][i] / data[i][i];
        for (unsigned int j{}; j < cols; ++j){
          //data[k-1][j] /= data[i][i];
          data[k][j] -= factor * data[i][j];
        }
    
    }
    this->printMatrix();
  }
  }
	
void relax(const vector<T>& y, int maxiter, double eps, double omega){
	vector<T> x(rows, 0);
	vector<T> New_x(rows);
	int iter = 0;
	T error = eps + 1.0;

	while (iter < maxiter && error > eps){
		error = 0;
		for (int i = 0; i < rows; i++){
			New_x[i] = function(y, x, i, omega);
			T diff = abs(New_x[i] - x[i]);
			if (diff > error){
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

        for (int i = 0; i < x.size(); i++){
		cout << "x" << i << "= " << x[i] << endl;
	}
	
}

T function(const vector<T>& y, const vector<T>& x, int i, double omega) {
        T sum = 0;
        for (int j = 0; j < columns; j++) {
            if (j != i) {
                sum += data[i][j] * x[j];
            }
        }

        T newX = (y[i] - sum) / data[i][i];
        newX = omega * newX + (1.0 - omega) * x[i];
        return newX;
}


int main() {
    std::vector<std::vector<double>> data{
        {1., -4., -2.},
        {3., 1., 1.},
        {3., -5., -6.}
    };
    std::vector<double> y {-3.,5.,-9.};
    
    Matrix<double> A(3, 3, data);
    std::cout << "Matrix A:" << std::endl;
    A.printMatrix();
    Matrix<double> C = A.inverse();
    C.printMatrix();
    A.relax(y, 100, 0.1, 0.25);
    return 0;
}
