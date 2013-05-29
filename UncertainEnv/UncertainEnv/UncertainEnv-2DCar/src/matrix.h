#ifndef __MATRIX_H__
#define __MATRIX_H__

#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <assert.h>

template <size_t _numRows, size_t _numColumns = 1> class Matrix {

private:
	double _elems[_numRows * _numColumns];

public:

	// Retrieval
	inline size_t numRows() const { 
		return _numRows; 
	}
	inline size_t numColumns() const { 
		return _numColumns; 
	}

	// Subscript operator
	inline double& operator () (size_t row, size_t column) { 
		return _elems[row * _numColumns + column]; 
	}
	inline double  operator () (size_t row, size_t column) const { 
		return _elems[row * _numColumns + column]; 
	}

	inline double& operator [] (size_t elt) { 
		return _elems[elt]; 
	}
	inline double  operator [] (size_t elt) const { 
		return _elems[elt]; 
	}

	// Reset to zeros
	inline void reset() { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] = double(0);
		}
	}

	// Submatrix
	template <size_t numRows, size_t numColumns>
	inline Matrix<numRows, numColumns> subMatrix(size_t row, size_t column) const {
		Matrix<numRows, numColumns> m;
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				m(i, j) = (*this)(row + i, column + j);
			}
		}
		return m;
	}

	// Insert
	template <size_t numRows, size_t numColumns>
	inline void insert(size_t row, size_t column, const Matrix<numRows, numColumns>& q) {
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				(*this)(row + i, column + j) = q(i, j);
			}
		}
	}

	// Unary minus
	inline Matrix<_numRows, _numColumns> operator-() const {
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = -_elems[i];
		}
		return m;
	}

	// Unary plus
	inline const Matrix<_numRows, _numColumns>& operator+() const { 
		return *this; 
	}

	// Equality
	inline bool operator==(const Matrix<_numRows, _numColumns>& q) const { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			if (_elems[i] < q._elems[i] || _elems[i] > q._elems[i]) {
				return false;
			}
		}
		return true;
	}

	// Inequality
	inline bool operator!=(const Matrix<_numRows, _numColumns>& q) const { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			if (_elems[i] < q._elems[i] || _elems[i] > q._elems[i]) {
				return true;
			}
		}
		return false;
	}

	// Matrix addition
	inline Matrix<_numRows, _numColumns> operator+(const Matrix<_numRows, _numColumns>& q) const { 
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = _elems[i] + q._elems[i];
		}
		return m;
	}
	inline const Matrix<_numRows, _numColumns>& operator+=(const Matrix<_numRows, _numColumns>& q) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] += q._elems[i];
		}
		return *this; 
	}

	// Matrix subtraction
	inline Matrix<_numRows, _numColumns> operator-(const Matrix<_numRows, _numColumns>& q) const { 
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = _elems[i] - q._elems[i];
		}
		return m;
	}
	inline const Matrix<_numRows, _numColumns>& operator-=(const Matrix<_numRows, _numColumns>& q) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] -= q._elems[i];
		}
		return *this; 
	}

	// Scalar multiplication
	inline Matrix<_numRows, _numColumns> operator*(double a) const { 
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = _elems[i] * a;
		}
		return m;
	}
	inline const Matrix<_numRows, _numColumns>& operator*=(double a) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] *= a;
		}
		return *this;
	}

	// Scalar division
	inline Matrix<_numRows, _numColumns> operator/(double a) const { 
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = _elems[i] / a;
		}
		return m;
	}
	inline const Matrix<_numRows, _numColumns>& operator/=(double a) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] /= a;
		}
		return *this;
	}

	// Matrix multiplication
	template <size_t numColumns>
	inline Matrix<_numRows, numColumns> operator*(const Matrix<_numColumns, numColumns>& q) const {
		Matrix<_numRows, numColumns> m;
		for (size_t i = 0; i < _numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				double temp = double(0);
				for (size_t k = 0; k < _numColumns; ++k) {
					temp += (*this)(i, k) * q(k, j);
				}
				m(i, j) = temp;
			}
		}
		return m;
	}

	inline const Matrix<_numRows, _numColumns>& operator*=(const Matrix<_numColumns, _numColumns>& q) { 
		return ((*this) = (*this) * q); 
	}

	// Matrix transpose
	inline Matrix<_numColumns, _numRows> operator~() const {
		Matrix<_numColumns, _numRows> m;
		for (size_t i = 0; i < _numColumns; ++i) {
			for (size_t j = 0; j < _numRows; ++j) {
				m(i, j) = (*this)(j, i);
			}
		}
		return m;
	}
};

/*template < >
class Matrix<1,1>: public Matrix
{
// Casting to double for 1x1 matrix
inline operator double() const {
return _elems[0];
}
};*/


// Scalar multiplication 
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator*(double a, const Matrix<_numRows, _numColumns>& q) { return q*a; }

// Matrix trace
template <size_t _size>
inline double tr(const Matrix<_size, _size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i){
		trace += q(i, i);
	}
	return trace;
}

// Matrix 1-norm
template <size_t _numRows, size_t _numColumns>
inline double norm(const Matrix<_numRows, _numColumns>& q) {
	double norm1 = double(0);
	for (size_t j = 0; j < _numColumns; ++j) {
		double colabssum = double(0);
		for (size_t i = 0; i < _numRows; ++i) {
			colabssum += abs(q(i,j));
		}
		if (colabssum > norm1) {
			norm1 = colabssum;
		}
	}
	return norm1;
}


// Identity matrix
template <size_t _size>
inline Matrix<_size, _size> identity() {
	Matrix<_size, _size> m;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(i, j) = (i == j ? double(1) : double(0));
		}
	}
	return m;
}

// Zero matrix
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> zeros() {
	Matrix<_numRows, _numColumns> m;
	m.reset();
	return m;
}

// Matrix determinant
template <size_t _size>
inline double det(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> m(q);
	double D = double(1);

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = std::abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// swap rows and columns
		if (k != max_row) {
			size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
			D = -D;
		}
		if (k != max_col) {
			size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
			D = -D;
		}

		D *= m(row_p[k], col_p[k]);
		if (D == double(0)) {
			return double(0);
		}

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
		}  
	}

	return D;
}



// Matrix inverse
template <size_t _size>
inline Matrix<_size, _size> operator!(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> m(q);
	Matrix<_size, _size> inv = identity<_size>();

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// swap rows and columns
		if (k != max_row) {
			size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
		}
		if (k != max_col) {
			size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
		}

		// eliminate column
		assert(m(row_p[k], col_p[k]) != double(0));
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
			for (size_t j = 0; j < k; ++j) {
				inv(row_p[i], row_p[j]) -= factor * inv(row_p[k], row_p[j]);
			}
			inv(row_p[i], row_p[k]) = -factor;
		} 
	}

	// Backward substitution
	for (size_t k = _size - 1; k != -1; --k) {
		double quotient = m(row_p[k], col_p[k]);
		for (size_t j = 0; j < _size; ++j) {
			inv(row_p[k], j) /= quotient;
		}

		for (size_t i = 0; i < k; ++i) {
			double factor = m(row_p[i], col_p[k]);
			for (size_t j = 0; j < _size; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		} 
	}

	// reshuffle result
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(col_p[i], j) = inv(row_p[i], j);
		}
	}

	return m; 
}

template <size_t _size>
inline void jacobi(const Matrix<_size, _size>& q, Matrix<_size, _size>& V, Matrix<_size, _size>& D) {
	D = q;
	V = identity<_size>();

	while (true) {
		double maximum = 0; size_t max_row = 0; size_t max_col = 0;
		for (size_t i = 0; i < _size; ++i) {
			for (size_t j = i + 1; j < _size; ++j) {
				if (abs(D(i,j)) > maximum) {
					maximum = abs(D(i,j));
					max_row = i;
					max_col = j;
				}
			}
		}

		if (maximum < 1e-10) {
			break;
		}

		double theta = (D(max_col, max_col) - D(max_row, max_row)) / (2 * D(max_row, max_col));
		double t = 1 / (abs(theta) + sqrt(theta*theta+1));
		if (theta < 0) t = -t;
		double c = 1 / sqrt(t*t+1); 
		double s = c*t;

		Matrix<_size, _size> R = identity<_size>();
		R(max_row,max_row) = c;
		R(max_col,max_col) = c;
		R(max_row,max_col) = s;
		R(max_col,max_row) = -s;

		// update D // 
		//std::cout << D << std::endl;
		//D = ~R * D * R;

		double temp1 = c*c*D(max_row, max_row) + s*s*D(max_col, max_col) - 2*c*s*D(max_row, max_col);
		double temp2 = s*s*D(max_row, max_row) + c*c*D(max_col, max_col) + 2*c*s*D(max_row, max_col);
		D(max_row, max_col) = 0;
		D(max_col, max_row) = 0;
		D(max_row, max_row) = temp1;
		D(max_col, max_col) = temp2;
		for (int j = 0; j < _size; ++j) {
			if ((j != max_row) && (j != max_col)) {
				temp1 = c * D(j, max_row) - s * D(j, max_col);
				temp2 = c * D(j, max_col) + s * D(j, max_row);
				D(j, max_row) = (D(max_row, j) = temp1);
				D(j, max_col) = (D(max_col, j) = temp2);
			}
		}
		//std::cout << D << std::endl << std::endl;


		V = V * R;
	} 
}


// Matrix exponentiation
#define _B0 1729728e1
#define _B1 864864e1
#define _B2 199584e1
#define _B3 2772e2
#define _B4 252e2
#define _B5 1512e0
#define _B6 56e0
#define _B7 1e0
#define _NORMLIM 9.504178996162932e-1

template <size_t _size>
inline Matrix<_size, _size> exp(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> A(q);
	int s = (int) std::max(double(0), ceil(log(norm(A)/_NORMLIM)*M_LOG2E)); 

	A /= pow(2.0,s);
	Matrix<_size, _size> A2(A*A);
	Matrix<_size, _size> A4(A2*A2);
	Matrix<_size, _size> A6(A2*A4);
	Matrix<_size, _size> U( A*(A6*_B7 + A4*_B5 + A2*_B3 + identity<_size>()*_B1) );
	Matrix<_size, _size> V( A6*_B6 + A4*_B4 + A2*_B2 + identity<_size>()*_B0 ); 
	Matrix<_size, _size> R7 = !(V - U)*(V + U);

	for (int i = 0; i < s; ++i) {
		R7 *= R7;
	}
	return R7;
}

// Input stream
template <size_t _numRows, size_t _numColumns>
inline std::istream& operator>>(std::istream& is, Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			is >> q(i,j);
		}
	}
	return is;
}

// Output stream
template <size_t _numRows, size_t _numColumns>
inline std::ostream& operator<<(std::ostream& os, const Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			//os << q(i,j) << "\t";
			printf_s("%24.24g ", q(i,j));
		}
		//printf_s("\n");
		os << std::endl;
	}
	return os;
}

// Output to file
template <size_t _numRows, size_t _numColumns>
inline void serialize(std::ofstream & ofs, const Matrix<_numRows, _numColumns>& q) 
{
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			ofs << std::setprecision(12) << " " << q(i,j);
		}
	}
}

template <size_t _numRows, size_t _numColumns>
inline void deserialize(std::ifstream& ifs, Matrix<_numRows, _numColumns>& q) 
{
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			ifs >> q(i,j);
		}
	}
}

// Hilbert matrix
template <size_t _size>
inline Matrix<_size, _size> hilbert() {
	Matrix<_size, _size> m;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(i, j) = double(1) / (double) (i + j + 1);
		}
	}
	return m;
}

#endif