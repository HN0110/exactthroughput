#ifndef _INITIAL_CONSTRUCTION_HPP_INCLUDED2_
#define _INITIAL_CONSTRUCTION_HPP_INCLUDED2_

#pragma once
//��l����
double UniformRandom(double x1, double x2);

//2�_�ԋ���
double TwoPdistance(double x1, double y1, double x2, double y2);

//Channel
double Pathloss(double dist, double wavelength);

double Shadowing();

double RayleighFading();

#endif