#pragma warning(disable : 4996)
#include <iostream>
#include <vector>
#include <math.h> 
#include <ctime>
using namespace std;
void qrFactorisation(vector<float> Q, vector<float> R, vector<float> matrix, int a);

void main()
{
	setlocale(LC_ALL, "Russian");
	int a = 0;
	cout << "Input matrix rank" << endl;
	cin >> a;
	srand(time(0));
	vector<float> matrix(a*a);// ({ 1,2,4,3,3,2,4,1,3 }); // базовая матрица, она будет подаваться на вход
	vector<float> Q(matrix.size());
	vector<float> R(matrix.size());
	cout << "Input matrix elements" << endl;
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < a; j++) {
			cin >> matrix[i*a + j];
		}
	}
	/*for (int i = 0; i < a; i++) {
		for (int j = 0; j < a; j++) {
			matrix[i*a + j] = rand() % 100 - rand() % 100;
		}
	}
	//cout << "Матрица, которую требуется разложить" << endl;
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < a; j++) {
			cout << matrix[i*a + j] << "\t";
		}
		cout << endl;
	}
	*/
	qrFactorisation(Q, R, matrix, a);
	cout << "runtime = " << clock() / 1000.0 << "ms" << endl; // время работы программы   
	system("pause");
}

void qrFactorisation(vector<float> Q, vector<float> mat, vector<float> matrix, int a) {
	vector<float> P(matrix.size());
	//Заполняем единичную матрицу
	for (int i = 0; i < a; i++) {
		Q[i * a + i] = 1;
	}

	for (int i = 0; i < a; i++) {
		double norm = 0; // норма темп вэлью
		vector<double> tempValue(a); // для записи взаимно ортогональных векторов
		for (int j = i; j < a; j++) {
			tempValue[j - i] = -matrix[j * a + i];
			norm += tempValue[j - i] * tempValue[j - i];
		}
		norm = sqrt(norm);

		if (tempValue[0] > 0) norm = -norm;
		tempValue[0] = tempValue[0] + norm;
		norm = 0;
		for (int j = 0; j < a - i; j++) {
			norm += tempValue[j] * tempValue[j];
		}
		norm = sqrt(norm);

		//Нормальзируем
		if (norm > 0) {
			for (int k = 0; k < tempValue.size(); k++) {
				tempValue[k] /= norm;
			}
		}
		//Заполняем P
		for (int k = 0; k < a - i; k++) {
			for (int l = 0; l < a - i; l++) {
				if (k == l) P[k * (a)+k] = 1 - 2 * tempValue[k] * tempValue[l];
				else P[k * a + l] = -2 * tempValue[k] * tempValue[l];
			}
		}



		//Заполняем R
		for (int k = i; k < a; k++) {
			for (int l = i; l < a; l++) {
				float tm = 0;
				for (int m = i; m < a; m++) {
					tm += P[(k - i) * a + m - i] * matrix[m * a + l];
				}
				mat[k * a + l] = tm;

			}
		}
		for (int k = i; k < a; k++) {
			for (int l = i; l < a; l++) {
				matrix[k * a + l] = mat[k * a + l];
			}
		}


		//Заполняем Q

		for (int k = 0; k < a; k++) {
			for (int l = i; l < a; l++) {
				float tm = 0;
				for (int m = i; m < a; m++) {
					tm += Q[k * a + m] * P[(m - i) * a + l - i];
				}
				mat[k * a + l] = tm;

			}
		}
		for (int k = 0; k < a; k++) {
			for (int l = i; l < a; l++) {
				Q[k * a + l] = mat[k * a + l];
			}
		}
	}

	cout << "P" << endl;
	for (int k = 0; k < a; k++) {
		for (int l = 0; l < a; l++) {
			cout << P[k*a + l] << "\t";
		}
		cout << endl;
	}
	cout << "Q" << endl;
	for (int k = 0; k < a; k++) {
		for (int l = 0; l < a; l++) {
			cout << Q[k*a + l] << "\t";
		}
		cout << endl;
	}
	cout << "R" << endl;
	for (int k = 0; k < a; k++) {
		for (int l = 0; l < a; l++) {
			cout << matrix[k*a + l] << "\t";
		}
		cout << endl;
	}


}
