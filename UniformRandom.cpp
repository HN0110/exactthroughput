#include <iostream>
#include <random>
#include <stdio.h>
#include <math.h>
#include <time.h>
using namespace std;

double UniformRandom(double x1, double x2)
{
	random_device rd;

	mt19937 mt(rd());

	//????x1????x2?????????l????
	uniform_real_distribution<> rand100(x1, x2);

	return rand100(mt);
}