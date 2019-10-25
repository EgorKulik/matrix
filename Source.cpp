#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

ifstream fin("input.txt");

class Matrix
{
public:
	Matrix()
	{
		cout << "Inputing data...";
		fin >> _size;
		cout << endl << endl << "X:	";
		for (int i = 0; i < _size; i++)
		{
			double current_input_x;
			fin >> current_input_x;
			cout.width(3);
			cout << current_input_x << " ";
			_input_x.push_back(current_input_x);
		}
		cout << endl << endl << "F(X):	";
		for (int i = 0; i < _size; i++)
		{
			double current_input_fx;
			fin >> current_input_fx;
			cout.width(3);
			cout << current_input_fx << " ";
			_input_fx.push_back(current_input_fx);
		}

		cout << endl << endl << "System of linear algebraic equations: " << endl << endl;
		for (int i = 0; i < _size; i++)
		{
			vector<double>equation;
			cout.width(3);
			cout << "1*a0 + ";
			equation.push_back(1);
			for (int j = 1; j < _size; j++)
			{
				double coef_x;
				coef_x = _input_x[i];
				cout.width(3);
				cout << pow(coef_x, j) << "*a" << j;
				if (j != _size - 1)
				{
					cout << " + ";
				}
				else
				{
					cout << " = ";
				}
				equation.push_back(pow(coef_x, j));
			}
			double current_fx;
			current_fx = _input_fx[i];
			cout.width(3);
			cout << current_fx << endl;
			_input_matrix_x.push_back(equation);
			_input_fx.push_back(current_fx);
		}

		cout << endl << endl << "Filling matrix... " << endl << endl;
		for (int i = 0; i < _size; i++)
		{
			vector<double>equation;
			cout.width(3);
			cout << "1 ";
			equation.push_back(1);
			for (int j = 1; j < _size; j++)
			{
				double coef_x;
				coef_x = _input_x[i];
				cout.width(3);
				cout << pow(coef_x, j) << " ";
				equation.push_back(pow(coef_x, j));
			}
			double current_fx;
			current_fx = _input_fx[i];
			cout.width(3);
			cout << "  " << current_fx << endl;
			_input_matrix_x.push_back(equation);
			_input_fx.push_back(current_fx);
		}
		cout << endl << "Coefficient matrix is full!" << endl << endl;
	}

	~Matrix()
	{
		/**/
	}

	void GetTriangleMatrix()
	{
		cout << endl;
		for (int k = 0; k < _size - 1; k++)
		{
			SwapStrings(k);
			ZeroingMatrix(k);
		}

		cout << "Triangle matrix: " << endl;
		for (int i = 0; i < _size; i++)
		{
			for (int j = 0; j < _size; j++)
			{
				cout.width(3);
				cout << _input_matrix_x[i][j] << " ";
			}
			cout.width(3);
			cout << "	  " << _input_fx[i] << endl;
		}
		cout << endl;
	}

	void SwapStrings(int x)
	{
		if (_input_matrix_x[x][x] == 0)
		{
			for (int i = x + 1; i < _size; i++)
			{
				if (_input_matrix_x[i][x] != 0)
				{
					vector<double>tmp_vector = _input_matrix_x[i];
					double tmp_value = _input_fx[i];
					_input_matrix_x[i] = _input_matrix_x[x];
					_input_fx[i] = _input_fx[x];
					_input_matrix_x[x] = tmp_vector;
					_input_fx[x] = tmp_value;
				}
			}
		}
	}

	void ZeroingMatrix(int x)
	{
		for (int i = x + 1; i < _size; i++)
		{
			double coefficient = (double)(_input_matrix_x[i][x] / _input_matrix_x[x][x]);
			for (int j = x; j < _size; j++)
			{
				_input_matrix_x[i][j] -= coefficient * _input_matrix_x[x][j];
			}
			_input_fx[i] -= coefficient*_input_fx[x];
		}
	}

	double GetDeterminant()
	{
		int l;
		double d;
		double sum11 = 1;
		double sum12 = 0;
		double sum21 = 1;
		double sum22 = 0;
		// находим детерминант
		for (int i = 0; i<_size; i++)
		{
			sum11 = 1;
			l = 2 * _size - 1 - i;
			sum21 = 1;
			for (int j = 0; j<_size; j++)
			{
				sum21 *= _input_matrix_x[j][l%_size];
				l--;
				sum11 *= _input_matrix_x[j][(j + i) % (_size)];
			}
			sum22 += sum21;
			sum12 += sum11;
		}
		d = sum12 - sum22;
		return d;
	}

	int GetSize()
	{
		return _size;
	}

	double GetCoef(int x, int y)
	{
		return _input_matrix_x[x][y];
	}

	double GetAns(int x)
	{

		return _input_fx[x];
	}

private:
	int _size;
	vector<vector<double>> _input_matrix_x;
	vector<double> _input_fx;
	vector<double> _input_x;
};

int main()
{
	Matrix M;
	vector<double>answers;
	M.GetTriangleMatrix();
	cout << "The solution:" << endl;
	answers.push_back(M.GetAns(M.GetSize() - 1) / M.GetCoef(M.GetSize() - 1, M.GetSize() - 1));
	for (int i = (M.GetSize() - 1) - 1; i >= 0; i--)
	{
		//здесь будет храниться (b)
		double sum_tmp = 0;
		//подставляем найденные уже корни, приводим уравнение к виду (y=kx+b)
		for (int j = 0; j < answers.size(); j++)
		{
			//находим коэффициент уже вычисленного корня
			double coef_local = M.GetCoef(i, M.GetSize() - 1 - j);
			//считаем значение одного слагаемого
			sum_tmp += answers[j] * coef_local;
		}
		double y_local = M.GetAns(i); // (y)
		double y_b_local = (M.GetAns(i) - sum_tmp); // (y-b)
		double k_local = M.GetCoef(i, i); // (k)
		double answer_local = y_b_local / k_local; // ( x=(y-b)/k )
		answers.push_back(answer_local);
	}
	for (int l = 0; l < answers.size(); l++)
	{
		cout << "a" << l << " = ";
		cout.width(3);
		cout << answers[answers.size() - l - 1] << endl;
	}
	cout << endl << "L = ";

	for (int l = 0; l < answers.size(); l++)
	{
		cout << answers[answers.size() - l - 1] << " * x^" << l;
		if (l != answers.size() - 1)
		{
			cout << " + ";
		}
	}
	cout << endl << endl << endl;
	return 0;
}
