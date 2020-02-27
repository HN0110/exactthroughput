
#ifndef _INITIAL_CONSTRUCTION_HPP_INCLUDED_
#define _INITIAL_CONSTRUCTION_HPP_INCLUDED_

#pragma once

//#define  oppai  Hatada 0.72;
#define SmallNum 3//スモールセルの数
#define InitialSmallCellRadius 150.0// スモールセルの半径

#define DeltaRadius 50.0 //電力制御時に狭める半径

#define MapSize 1000 //マップの最大距離
#define GridNum 100  //グリッドの数

#define SensorNum 10 //センサーの数

#define DeltaTx_P 0.1 //電力制御時に下げる

//#define Tx_P 30.0			/* 送信電力[dBm] */
#define Pout 0.01			/* アウテージ確率 */

#define Margin 0			/* 干渉到達エリアのマージン */

#define PathlossExponent 3.5	/* パスロス減衰係数 */
#define RefDistance	 10		/* パスロス参照距離10[m] */

#define MinTale		 -70.0 /* シャドウイングのpdfの横軸最小値 */
#define MaxTale		 100.0	/* シャドウイングのpdfの横軸最大値 */
#define ShadowMean	 0.0	/* シャドウイングの対数正規分布の平均 */
#define ShadowSD	 8.0	/* シャドウイングの対数正規分布の標準偏差 (standard deviation : SD) */
#define DeltaD		 0.1    /* 微小幅 */

#define DesiredSINR	 5.0	/* 所望SINR[dB] */

#define MinTale1	 -40.0  /* SIRのpdfの横軸最小値 */
#define MaxTale1	 60.0   /* SIRのpdfの横軸最大値 */

#define PI			 3.141592653589793
#define loop		 1000

#define SNR			 30.0		 /* 所望SNR */
#define Boltzmann	 1.38064852  /*ボルツマン定数 10^(-23)なし */
#define Kelvin		 293.15		 /* 導体の温度 20℃ */
#define FB			 100.0		 /* 周波数帯域幅 100MHz 10^(6)なし */

#define NOISE_DB -95.0
#define INITIAL_BANDWIDTH 100000000.0
#define BAND_MAX 4.7
#define BAND_MIN 4.6
#define DIVIDED_NUM 2

#define GB 5 /* ガードバンド 5MHz(MHzはMain上でたす)*/

#define DivideNumMax 10
#define DivideNumMin 2
typedef struct
{
	double x;//mapのx座標
	double y;//mapのy座標
	double dist[SmallNum];//スモールセル間の距離
	int InterferenceArrivalAreaFlag[SmallNum];
	int PowerControlDoneFlag;
	int AllSpectrumBandFlag;
	double UseBand;//使用中心周波数
	double AlreadyAllocatedBandFlag; //すでに割り当てたバンドを割り当てた場合はフラグ立て
	double Throughput;
	double Throughput_GB;
	double PowerControlThroughput;
	double BandDivisionThroughput;
}SmallCell;

typedef struct
{
	double x;//最悪点のx座標
	double y;//最悪点のy座標
	double difference_x; //最悪点と着目スモールセルのx座標の差分
	double difference_y; //最悪点と着目スモールセルのy座標の差分

}InterferenceWorstPoint;

typedef struct
{
	double x; //センサーのx座標
	double y; //センサーのy座標

}Sensor;

typedef struct
{
	int Num;
	int flag;
	double dist[4];
	double AveSINR[SmallNum];
	double tmpAveSINR[SmallNum];
	double AveOwn[SmallNum];
	double AveInter[SmallNum][SmallNum];
}Grid;

#endif