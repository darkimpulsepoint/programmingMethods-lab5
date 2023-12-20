#include <iostream>
#include <ctime>
#include <cmath>

#define MATRIX_TYPE Fraction
/*
dont use int!
preferable to use double | Fraction
*/

using namespace std;

template <typename T>
class Matrix
{

private:
    T **elements;
    int rows, cols;

public:
    Matrix(int rows, int cols);

    template <typename C>
    Matrix(C **els, int rows, int cols);

    void randomFill();

    int getRows() const;

    int getCols() const;

    T getValue(int row, int col) const;

    template <typename C>
    friend ostream &operator<<(ostream &os, const Matrix<C> &myMatrix);

    void insertValue(int row, int col, T value);
};

template <typename T>
class QuadricMatrix : public Matrix<T>
{
public:
    T determinant() const;

    QuadricMatrix(int rank);
};

class Fraction
{

private:
    int m;
    int n;

    void normalize();

    int gcd(int m, int n) const;

public:
    Fraction(int m, int n);

    Fraction(int m);

    Fraction();

    Fraction operator+(const Fraction &second) const;

    Fraction operator-(const Fraction &second) const;

    Fraction operator*(const Fraction &second) const;

    Fraction operator/(const Fraction &second) const;

    bool operator==(const Fraction &second) const;

    friend ostream &operator<<(ostream &, const Fraction &fr);

    friend istream &operator>>(istream &, Fraction &fr);

    operator bool() const;
};

template <typename T>
class QuadrSysLinAlgEq : public Matrix<T>
{

public:
    QuadrSysLinAlgEq(int size);

    template <typename C>
    QuadrSysLinAlgEq(C **els, int size);

    QuadrSysLinAlgEq();

    Matrix<T> getSolution() const;

    bool isSolutionExist() const;
};

template <typename T>
Matrix<T> QuadrSysLinAlgEq<T>::getSolution() const
{
    if (!isSolutionExist())
    {
        return Matrix<T>(0, 0);
    }

    T **temp;

    temp = new T *[this->getRows()];

    for (size_t i = 0; i < this->getRows(); i++)
    {
        temp[i] = new T[this->getCols()];
        for (size_t j = 0; j < this->getCols(); j++)
        {
            temp[i][j] = this->getValue(i, j);
        }
    }

    for (size_t k = 1; k < this->getRows(); k++)
    {
        for (size_t j = k; j < this->getRows(); j++)
        {
            if (temp[k - 1][k - 1] == T(0))
            {
                continue;
            }
            T m = temp[j][k - 1] / temp[k - 1][k - 1];

            for (size_t i = 0; i < this->getRows() + 1; i++)
            {
                temp[j][i] = temp[j][i] - m * temp[k - 1][i];
            }
        }
    }
    Matrix<T> answer = Matrix<T>(1, this->getRows());

    for (int i = this->getRows() - 1; i >= 0; i--)
    {
        answer.insertValue(0, i, temp[i][this->getRows()] / temp[i][i]);
        for (int c = this->getRows() - 1; c > i; c--)
        {
            answer.insertValue(0, i, answer.getValue(0, i) - temp[i][c] * answer.getValue(0, c) / temp[i][i]);
        }
    }

    return answer;
}

template <typename T>
template <typename C>
Matrix<T>::Matrix(C **els, int rows, int cols)
{
    this->rows = rows;
    this->cols = cols;
    this->elements = new T *[rows];
    for (size_t i = 0; i < rows; i++)
    {
        this->elements[i] = new T[cols];
    }

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            this->elements[i][j] = els[i][j];
        }
    }
}

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &myMatrix)
{
    for (int i = 0; i < myMatrix.rows; ++i)
    {
        for (int j = 0; j < myMatrix.cols; ++j)
        {
            cout << myMatrix.getValue(i, j) << " ";
        }
        cout << endl;
    }
    return os;
}

template <typename T>
T Matrix<T>::getValue(int row, int col) const
{
    return this->elements[row][col];
}

template <typename T>
int Matrix<T>::getCols() const
{
    return this->cols;
}

template <typename T>
int Matrix<T>::getRows() const
{
    return this->rows;
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols)
{
    this->rows = rows;
    this->cols = cols;
    elements = new T *[rows];

    for (size_t i = 0; i < rows; i++)
    {
        elements[i] = new T[cols];
        for (size_t j = 0; j < cols; j++)
        {
            elements[i][j] = T();
        }
    }
}

template <typename T>
void Matrix<T>::randomFill()
{
    srand(time(NULL));

    for (size_t row = 0; row < this->rows; row++)
    {
        for (size_t col = 0; col < this->cols; col++)
        {
            T el = T(rand() % (2 * this->rows + 1) - this->rows);

            if (el == T(0))
            {
                this->elements[row][col] = el + T(1);
            }
            else
            {
                this->elements[row][col] = el;
            }
        }
    }
}

template <typename T>
void Matrix<T>::insertValue(int row, int col, T value)
{
    if (row > rows - 1 || row < 0 || col > cols - 1 || col < 0)
    {
        return;
    }

    this->elements[row][col] = value;
}

template <typename T>
QuadricMatrix<T>::QuadricMatrix(int size) : Matrix<T>(size, size) {}

template <typename T>
T QuadricMatrix<T>::determinant() const
{
    T det = 0;

    if (this->getRows() == 1)
    {
        return this->getValue(0, 0);
    }

    for (size_t i = 0; i < this->getCols(); i++)
    {
        QuadricMatrix temp = QuadricMatrix(this->getRows() - 1);

        for (size_t j = 1; j < this->getRows(); j++)
        {
            for (size_t k = 0; k < i; k++)
            {
                temp.insertValue(j - 1, k, this->getValue(j, k));
            }

            for (size_t k = i + 1; k < this->getCols(); k++)
            {
                temp.insertValue(j - 1, k - 1, this->getValue(j, k));
            }
        }
        det = det + T(pow(-1, i)) * this->getValue(0, i) * temp.determinant();
    }
    return det;
}

template <typename T>
QuadrSysLinAlgEq<T>::QuadrSysLinAlgEq(int size) : Matrix<T>(size, size + 1) {}

template <typename T>
QuadrSysLinAlgEq<T>::QuadrSysLinAlgEq() : Matrix<T>(0, 0) {}

template <typename T>
template <typename C>
QuadrSysLinAlgEq<T>::QuadrSysLinAlgEq(C **els, int size) : Matrix<C>(els, size, size + 1) {}

template <typename T>
bool QuadrSysLinAlgEq<T>::isSolutionExist() const
{
    QuadricMatrix<T> temp = QuadricMatrix<T>(this->getRows());
    for (size_t row = 0; row < this->getRows(); row++)
    {
        for (size_t col = 0; col < this->getRows(); col++)
        {
            temp.insertValue(row, col, this->getValue(row, col));
        }
    }

    return temp.determinant();
}

QuadrSysLinAlgEq<MATRIX_TYPE> getMatrixByDefaultPattern()
{
    MATRIX_TYPE **els;
    els = new MATRIX_TYPE *[3];
    for (int i = 0; i < 3; i++)
    {
        els[i] = new MATRIX_TYPE[4];
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            els[i][j] = MATRIX_TYPE(10 - (i * 2 + j * j) * i + 2 * j);
        }
    }

    QuadrSysLinAlgEq<MATRIX_TYPE> pregeneratedMatrix = QuadrSysLinAlgEq<MATRIX_TYPE>(els, 3);

    return pregeneratedMatrix;
}

Fraction::Fraction()
{
    this->m = 0;
    this->n = 1;
}

Fraction::Fraction(int m, int n)
{
    this->m = m;
    this->n = n;

    this->normalize();
}

Fraction::Fraction(int m)
{
    this->m = m;
    this->n = 1;
}

int Fraction::gcd(int x, int y) const
{
    while (y != 0)
    {
        int c = x % y;
        x = y;
        y = c;
    }
    return x;
}

void Fraction::normalize()
{
    int reducer = gcd(abs(this->m), abs(this->n));

    this->m = m / reducer;
    this->n = n / reducer;

    if (n < 0)
    {
        this->m = m * (-1);
        this->n = n * (-1);
    }
}

Fraction::operator bool() const
{
    return (this->m != 0);
}

Fraction Fraction::operator+(const Fraction &second) const
{
    return Fraction(this->m * second.n + this->n * second.m, this->n * second.n);
}

Fraction Fraction::operator/(const Fraction &second) const
{
    return Fraction(this->m * second.n, this->n * second.m);
}

Fraction Fraction::operator-(const Fraction &second) const
{
    return Fraction(this->m * second.n - this->n * second.m, this->n * second.n);
}

Fraction Fraction::operator*(const Fraction &second) const
{
    return Fraction(this->m * second.m, this->n * second.n);
}

bool Fraction::operator==(const Fraction &second) const
{
    return this->m * second.n == this->n * second.m;
}

ostream &operator<<(ostream &os, const Fraction &fr)
{
    if (fr.n != 1)
    {
        return os << fr.m << "/" << fr.n;
    }
    else
    {
        return os << fr.m;
    }
}

istream &operator>>(istream &in, Fraction &fr)
{
    int numerator;

    in >> numerator;

    fr = Fraction(numerator, 1);

    return in;
}

int askSize()
{
    int size = 0;

    cout << "Enter matrix size (non zero)" << endl;
    cin >> size;

    return size > 0 ? size : 1;
}

int main()
{
    QuadrSysLinAlgEq<MATRIX_TYPE> q, pregenerated = getMatrixByDefaultPattern();

    int choice = 0;
    while (true)
    {

        cout << "1 - use pregenerated matrix\n2 - generate matrix\n3 - fill yourself\n0 - exit\n";
        cin >> choice;

        Fraction temp;
        switch (choice)
        {

        default:
            cout << "[-] No match case! Choosing pregenerated matrix..." << endl;
        case 1:
            q = pregenerated;
            break;

        case 2:
            q = QuadrSysLinAlgEq<MATRIX_TYPE>(askSize());
            q.randomFill();
            break;

        case 3:
            q = QuadrSysLinAlgEq<MATRIX_TYPE>(askSize());
            for (size_t i = 0; i < q.getRows(); i++)
            {
                for (size_t j = 0; j < q.getCols(); j++)
                {
                    cout << "Enter [" << i + 1 << "][" << j + 1 << "] element" << endl;
                    cin >> temp;
                    q.insertValue(i, j, temp);
                }
            }
            break;

        case 0:
            exit(0);
            break;
        }

        cout << "choosen matrix:" << endl;
        cout << q << endl;

        Matrix<MATRIX_TYPE> solution = q.getSolution();

        cout << "solution:" << endl;

        if (!solution.getRows())
        {
            cout << "the system has no or an infinite number of solutions" << endl;
        }
        else
        {
            cout << solution << endl;
        }
    }
}