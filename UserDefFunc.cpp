

#include <iostream>
#include <random>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Parameter.h"
using namespace std;

//二点間の距離
double TwoPdistance(double x1, double y1, double x2, double y2)
{
	double dist;

	dist = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

	return dist;
}

//パスロス
double Pathloss(double dist, double wavelength)
{
	return 20 * log10(4 * PI * RefDistance / wavelength) + 10 * PathlossExponent * log10(dist / RefDistance);
}

//シャドウイング
double Shadowing()
{
	random_device rd;

	mt19937 mt(rd());

	normal_distribution<> norm(0.0, 8.0);

	return norm(mt);

}

//レイリーフェージング
double RayleighFading()
{
	double H = 0.0;

	random_device rd;

	mt19937 mt(rd());

	//範囲x1以上x2未満の一様乱数
	uniform_real_distribution<> rand100(0.0, 1.0000000000000000000000000000000000000000001);


	return -log(rand100(mt));
}
