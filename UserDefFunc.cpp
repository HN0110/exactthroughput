

#include <iostream>
#include <random>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Parameter.h"
using namespace std;

//��_�Ԃ̋���
double TwoPdistance(double x1, double y1, double x2, double y2)
{
	double dist;

	dist = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

	return dist;
}

//�p�X���X
double Pathloss(double dist, double wavelength)
{
	return 20 * log10(4 * PI * RefDistance / wavelength) + 10 * PathlossExponent * log10(dist / RefDistance);
}

//�V���h�E�C���O
double Shadowing()
{
	random_device rd;

	mt19937 mt(rd());

	normal_distribution<> norm(0.0, 8.0);

	return norm(mt);

}

//���C���[�t�F�[�W���O
double RayleighFading()
{
	double H = 0.0;

	random_device rd;

	mt19937 mt(rd());

	//�͈�x1�ȏ�x2�����̈�l����
	uniform_real_distribution<> rand100(0.0, 1.0000000000000000000000000000000000000000001);


	return -log(rand100(mt));
}
