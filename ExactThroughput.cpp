/*************

c++ -o Main Main.cpp UserDefFunc.cpp -std=c++0x -g
***********/




#include "UserDefFunc.h"
#include "Parameter.h"

#include <iomanip>
#include <fstream>  //ifstream, ofstream
#include <algorithm>
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <map>
#include <sstream>

#define SIZE_OF_ARRAY(array) (sizeof(array)/sizeof(array[0]))

using namespace std;

namespace {
	void printElems(const double* begin, const double* end)
	{
		for (const double* p = begin; p != end; ++p) {
			cout << *p << ends << "\n";
		}
		cout << endl;
	}
}

int main()
{

    //csv出力
    string filename1 = "cell4_3-1-1.csv";
    ofstream writing_file1;
    writing_file1.open(filename1, ios::out);


    //スループットが最大の場合のグラフ導出

    int i;
    double c_max = 0.0, c, cf, snr, B, noise_mw, noise_db, wavelength, rx;

    string filename = "otameshi_th.csv";

    ofstream writing_file;
    writing_file.open(filename, ios::out);

    cf = (BAND_MAX * pow(10.0, 9.0) + BAND_MIN * pow(10.0, 9.0)) / 2.0;

    wavelength = 3.0 * pow(10.0, 8.0) / cf;
    rx = 41.9540253107965 - Pathloss(150.0, wavelength);
    snr = pow(10.0, rx / 10.0) / pow(10.0, NOISE_DB / 10.0);

    c = INITIAL_BANDWIDTH * log2(1 + snr);

    writing_file1 << "wavelength," << wavelength << ",centerfrequency," << cf << ",rx," << rx << ",snr," << snr << ",c," << endl;
    //writing_file << "number_of_cells" << "," << "max_total_througput" << endl;

    for (i = 0; i < SmallNum; i++)
    {
        c_max += c / pow(10.0, 6.0);
        writing_file << i << "," << c_max << endl;
    }

    ifstream ifs(filename);
    string str;
    while(getline(ifs,str))
    {
        string tmp;
        istringstream stream(str);
        while(getline(stream,tmp, ','))
        {
            cout<< tmp << endl;
        }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //電力制御なし・帯域分割なしの周波数共用可能かの判定

	int j, k, l, m, n, x, y, z;
	int o = 0;

	int num1, num2;
	double* pT;

	double Rx_P; //着目スモールセルの基地局からの受信信号電力
	//double Rx_I; //スモールセルからの受信信号電力
	//double SIR[SmallNum];  //SIR


	map<string, int> MicroDist;
	string names[] = { "Micro01", "Micro23", "Micro31" };


	SmallCell SmallCell[SmallNum];

	//着目スモールセル内の最も干渉が大きくなる点

	InterferenceWorstPoint WorstPoint[SmallNum];


	Grid Grid[GridNum][GridNum];

	//スモールセルの基地局位置

	//SmallCell[0].x = 125.0;
	//SmallCell[0].y = 125.0;

	//SmallCell[1].x = 125.0;
	//SmallCell[1].y = 285.0;

	//SmallCell[2].x = 255.0;
	//SmallCell[2].y = 285.0;



	//送信電力の制御
	double Tx_P[SmallNum];

	//セル半径
	double SmallCellRadius[SmallNum] = { 0.0 };

	//SIRthと所望SIRの差分
	double DeltaP[SmallNum] = { 0.0 };

	//センサー
	Sensor Sensor[SmallNum][SensorNum];

	//干渉スモールセルと着目スモールセル間の距離の中で最も短い距離
	double MinDist[SmallNum] = { 0.0 };

	//干渉スモールセルに最も近い干渉スモールセルのナンバー
	int MinDistNum[SmallNum] = { 0 };

	//最悪点が含まれるGrid座標
	int WorstGridNum_x[GridNum];
	int WorstGridNum_y[GridNum];

	//最悪点が含まれるメッシュの始まりと終わりのx,y座標
	double WorstMeshStart_x[SmallNum];
	double WorstMeshFinish_x[SmallNum];

	double WorstMeshStart_y[SmallNum];
	double WorstMeshFinish_y[SmallNum];

	//センサと着目スモールセル基地局間の距離
	double SensorSmallCellDist[SmallNum][SensorNum];

	//ある微小区間で発生した回数
	double NumofTimes[3410] = { 0.0 };
	double CDFofTimes[3410] = { 0.0 };
	int NumofTime = 0;

	//自局信号電力，他局信号電力
	double OwnPower = 0.0;
	double OtherPower_mW = 0.0;
	double OtherPower_dBm = 0.0;
	double tmp = 0.0;
	double SINR = 0.0;
	double MinGridAveSINR = 0.0;
	int MinGridx = 0;
	int MinGridy = 0;
	double AveSINR = 0.0;
	//double AveOwn = 0.0;
	//double AveInter1 = 0.0;
	//double AveInter2 = 0.0;
	int Aveflag = 0;

	double tmpRSSI = 0.0;
	double tmpOwn = 0.0;
	double tmp2 = 0.0;
	double tmpOther_mw = 0.0;
	double tmpInterNoise = 0.0;
	double tmpOther_dBm = 0.0;
	double tmpSINR = 0.0;
	double tmpMinGridAveSINR = 0.0;
	int tmpMinGridx[SmallNum] = { 0 };
	int tmpMinGridy[SmallNum] = { 0 };

	//double RayleighFadingValue[loop * SensorNum * SmallNum] = { 0.0 };

	double CDF = 0.0;
	double SINRth = 0.0;
	int SINRthFlag = 0;

	int DeltaPFlag = 0;

	//周波数共用失敗カウンター
	int SpectrumSharingImpossible_nonpc[SmallNum] = { 0 };
	int SpectrumSharingImpossible_pc[SmallNum] = { 0 };
	int SpectrumSharingImpossible_nonpc_initial[SmallNum] = { 0 };

	//アウテージ確率
	double pout[SmallNum] = { 0.0 };
	double pout_nonpc[SmallNum] = { 0.0 };
	double pout_nonpc_initial[SmallNum] = { 0.0 };
	int PoutInitialFlag = 0;

	//平均アウテージ確率
	double Avepout[SmallNum] = { 0.0 };

	//SNR関連 熱雑音導出など
	double Noise_dBm;
	double Noise_mW;

	//干渉電力足す雑音電力
	double InterferenceNoise = 0.0;

	//平均干渉電力
	double AveIsum = 0.0;

	int numsum = 0;
	double instantaneousSINR[loop * SensorNum] = { 0.0 };

	//アウテージ出すようの瞬時SINR導出用電力配列
	double WorstMeshNumX[SmallNum] = { 0.0 };
	double WorstMeshNumY[SmallNum] = { 0.0 };
	double CellAveOwn[SmallNum] = { 0.0 };
	double CellAveInter[SmallNum][SmallNum] = { 0.0 };

	double InterferenceArrivalAreaRadius = 0.0;
	int AddFlag[SmallNum] = { 0 };

	int PossibleSharingCounter = 0;
	int PossibleSharingCounter_NonPC = 0;
	double WorstMeshCenterx[SmallNum] = { 0.0 };
	double WorstMeshCentery[SmallNum] = { 0.0 };
	double WorstMeshCenterInterCellDist[SmallNum] = { 0.0 };
	int PowerControlDoneFlag[SmallNum] = { 0 };
	double WorMeshCentMinDist = 0.0;
	int WorMeshCentMinDistNum = 0;
	double InitialTx_P = 0.0;
	int NumofPowerControlDoneFlag = 0;
	int AllCellPowerControlCounter = 0;
	int NonPowerControlCounter = 0;
	int NotSharingCellFlag[SmallNum] = { 0 };
	double InitialDeltaP[SmallNum] = { 0.0 };
	double ComparisonDeltaP = 0.0;
	int DoNotSharingCellFlag[SmallNum] = { 0 };
	double InitialCellAveInter[SmallNum][SmallNum] = { 0.0 };
	double TotalCellAveInter[SmallNum] = { 0.0 };
	double ComparisonTotalCellAveInter = 0.0;
	int MaxTotalCellAveInterNum = 0;
	int FirstTotalCellAveInterFlag = 0;
	int FirstTotalCellAveInterNum = 0;
	int DoNotSharingCellNum = 0;
	int SharingAchievementFlag = 0;
	int DoNotSharingFlag = 0;
	double beforeTx_P = 0.0;
	double beforeCellRadius = 0.0;

	double cmin = 0.0;

    double rx_edge[SmallNum] = {0.0}; //セル端のスループットを算出する際に用いる
    double snr_edge[SmallNum] = {0.0};
    double c_edge[SmallNum] = {0.0};
    double c_sum = 0.0;
    double GB_c_sum = 0.0;

    double f1 = 0.0, f2 = 0.0;  //帯域分割後の中心周波数
    double f[DivideNumMax] = { 0.0 };   //各分割バンドの中心周波数
    double bandwitdh = 0.0; //スループット計算用の帯域幅
    double SumAveInter = 0.0;
    double numf1 = 0.0;
    double numf2 = 0.0;
    int    numf[DivideNumMax] = { 0 };
    double numcf = 0.0;
    int    DividedNum = 0;
    double CellAveInter_nofade = 0.0;
    double AveSum_nofade = 0.0;
    double AveSINR_nofade = 0.0;
    double InterferenceNoise_nofade = 0.0;
    double OwnPower_nofade = 0.0;
    int pc_count = 0;
    int breakflag = 0;

    int OUTAGE = 0;
    double c_sum_dp = 0.0; //spectrumbanddivideallocationとpowercontrol両方使った場合の合計スループット


    OUTAGE = Pout * loop * SensorNum - 1;
    cout << "OUTAGE," << OUTAGE << endl;

	//初期設定---------------------------------------------------------------------------------------------------------------------------------------------------------------
	//セルの設定
	for (i = 0; i < SmallNum; i++)
	{
		//セルの配置
		SmallCell[i].x = UniformRandom(0.00000000000000000000000001, 1000.0);
		SmallCell[i].y = UniformRandom(0.00000000000000000000000001, 1000.0);
        SmallCell[i].AlreadyAllocatedBandFlag = 0;

        writing_file1 << setprecision(20) << "cellx," << SmallCell[i].x << ",celly," << SmallCell[i].y << endl;
	}

    SmallCell[0].x = 234.7835635;
    SmallCell[0].y = 77.2975318;

    SmallCell[1].x = 625.094498599999;
    SmallCell[1].y = 892.554598199999;

    SmallCell[2].x = 884.5234791;
    SmallCell[2].y = 785.4373726;
/*
    SmallCell[3].x = 884.5234791;
    SmallCell[3].y = 785.4373726;

    SmallCell[4].x = 787.345887;
    SmallCell[4].y = 459.865465;
*/
	for (i = 0; i < SmallNum; i++)
	{
        writing_file1 << setprecision(20) << "saihaichi,cellx," << SmallCell[i].x << ",celly," << SmallCell[i].y << endl;
	}

    //初期のセル半径決め
    for (i = 0; i < SmallNum; i++)
    {
        SmallCellRadius[i] = InitialSmallCellRadius;
    }

    //雑音電力
    Noise_dBm = -95.0;
	Noise_mW = pow(10, Noise_dBm / 10.0);

    //cf = (BAND_MAX * pow(10.0, 9.0) + BAND_MIN * pow(10.0, 9.0)) / 2.0;
    //cf = 4625000000;
    //wavelength = 3.0 * pow(10.0, 8.0) / cf;

    //各スモールセルの初期送信電力決め
    for (i = 0; i < SmallNum; i++)
    {
        Tx_P[i] = SNR + Pathloss(SmallCellRadius[i], wavelength) + Noise_dBm;
		SmallCellRadius[i] = RefDistance * pow(10.0, (Tx_P[i] - Noise_dBm - SNR - 20 * log10(4 * PI * RefDistance / wavelength))/ (RefDistance * PathlossExponent));

        //初期送信電力をとっておく
        InitialTx_P = Tx_P[i];
        writing_file1 << "SmallCellRadius," << SmallCellRadius[i] << endl;
        writing_file1 << "Tx_P," << i << "," << Tx_P[i] << endl;

    }

    //干渉到達エリアの半径決め
    InterferenceArrivalAreaRadius = RefDistance * pow(10.0, (InitialTx_P - Noise_dBm - 20.0 * log10(4.0 * PI * RefDistance / wavelength)) / (10.0 * PathlossExponent));
    writing_file1 << "InterferenceArrivalAreaRadius," << InterferenceArrivalAreaRadius << endl;

    //セルに属するメッシュの判定フラグの初期化
    for (i = 0; i < GridNum; i++)
    {
        for (j = 0; j < GridNum; j++)
        {
            for (k = 0; k < SmallNum; k++)
            {
                Grid[i][j].flag = 0;
                Grid[i][j].AveSINR[k] = 0.0;
                Grid[i][j].tmpAveSINR[k] = 0.0;
                Grid[i][j].AveOwn[k] = 0.0;
            }
        }
    }

    for (i = 0; i < GridNum; i++)
    {
        for (j = 0; j < GridNum; j++)
        {
            for (k = 0; k < SmallNum; k++)
            {
                for (l = 0; l < SmallNum; l++)
                {
                    Grid[i][j].AveInter[l][k] = 0.0;
                }
            }
        }
    }

    for (i = 0; i < SmallNum; i++)
    {
        CellAveOwn[i] = 0.0;
        WorstMeshCenterInterCellDist[i] = 0.0;
        SmallCell[i].PowerControlDoneFlag = 0;

        for (j = 0; j < SmallNum; j++)
        {
            CellAveInter[i][j] = 0.0;
            InitialCellAveInter[i][j] = 0.0;

            SmallCell[i].InterferenceArrivalAreaFlag[j] = 0;
        }
    }

    //セルごとに最悪メッシュの導出 それらの情報を保持---------------------------------------------------------------------------------------------------------
    for (i = 0; i < SmallNum; i++)
    {
        writing_file1 << "\nSmallNum," << i << endl;

        //干渉信号到達エリア内に着目スモールセル以外のセルがあるか判定
        for (j = 0; j < SmallNum; j++)
        {
            SmallCell[i].dist[j] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, SmallCell[j].x, SmallCell[j].y);

            writing_file1 << "inter-CellDist," << SmallCell[i].dist[j] << endl;

            //cout << "Radius" << InterferenceArrivalAreaRadius - SmallCellRadius[i] << "\n";
            //着目セルと着目セル以外のセルの距離と干渉信号到達エリアの半径を比較
            if ((SmallCell[i].dist[j] < (InterferenceArrivalAreaRadius + Margin - SmallCellRadius[i])) && (i != j))
            {
                SmallCell[i].InterferenceArrivalAreaFlag[j] = 1;

                //cout << "InterferenceArrivalArea" << j << "\n";
                writing_file1 << "InterferenxeSinalArea_naino_smallcell," << j << ",SmallCell[i].InterferenceArrivalAreaFlag[j]," << SmallCell[i].InterferenceArrivalAreaFlag[j] << endl;
                //cout << SmallCell[i].InterferenceArrivalAreaFlag[j] << "\n";

            }
        }


        //着目スモールセルの最悪点を求める
        //「全メッシュの四隅の点とセルの中心点間の距離」，「スモールセルの半径」を比較して半径よりも距離が短い点が1つでもあればスモールセル内にあるメッシュとして旗を立てる
        for (j = 0; j < GridNum; j++)
        {
            for (k = 0; k < GridNum; k++)
            {
                //四隅分の距離を導出
                Grid[j][k].dist[0] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, j  * (MapSize / GridNum), k * (MapSize / GridNum));
                Grid[j][k].dist[1] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, (j + 1) * (MapSize / GridNum), k * (MapSize / GridNum));
                Grid[j][k].dist[2] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, j * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));
                Grid[j][k].dist[3] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, (j + 1) * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));
                //fprintf(fpA, "GridNum, %d, %d, SmallCellNum, %d, dist, %f\n", j, k, i, Grid[j][k].dist[0]);
                //fprintf(fpA, "GridNum, %d, %d, SmallCellNum, %d, dist, %f\n", j, k, i, Grid[j][k].dist[1]);
                //fprintf(fpA, "GridNum, %d, %d, SmallCellNum, %d, dist, %f\n", j, k, i, Grid[j][k].dist[2]);
                //fprintf(fpA, "GridNum, %d, %d, SmallCellNum, %d, dist, %f\n", j, k, i, Grid[j][k].dist[3]);


                //四隅の点と中心点間距離，スモールセルの半径を比較，半径以下だったら旗を立てる
                if ((SmallCellRadius[i] >= Grid[j][k].dist[0]) || (SmallCellRadius[i] >= Grid[j][k].dist[1]) || (SmallCellRadius[i] >= Grid[j][k].dist[2]) || (SmallCellRadius[i] >= Grid[j][k].dist[3]))
                {
                    Grid[j][k].flag = 1;
                    //cout << "Gridcoordinates(" << j << "," << k << ")," << i << "numberSmallCellCenter(" << SmallCell[i].x << "," << SmallCell[i].y << ")distance" << Grid[j][k].dist[0] << "\n";
                    //fprintf(fpA, "Gridcoordinates, %d, %d, %d, numberSmallCellCenter, %f, %f, distance, %f\n", j, k, i, SmallCell[i].x, SmallCell[i].y, Grid[j][k].dist[0]);
                }


                //旗が1となる（最悪グリッドとなりうる可能性がある）グリッドでセンサをばらまき試行を重ね，SINRの平均が最も低いグリッドを探す

                //条件式（旗が1ならば実行）
                if (Grid[j][k].flag == 1)
                {
                    //センサをSensorNum個だけばらまく
                    for (l = 0; l < SensorNum; l++)
                    {
                        //センサはグリッド座標(j, k)を使って一様乱数でまく
                        Sensor[i][l].x = UniformRandom(j * (MapSize / GridNum), (j + 1) * (MapSize / GridNum));
                        Sensor[i][l].y = UniformRandom(k * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));
                        //fprintf(fpA, "SensorCoordinates, %f, %f\n", Sensor[i][l].x, Sensor[i][l].y);
                        //cout << "sensor_x" << Sensor[i][l].x << "sensor_y" << Sensor[i][l].y << "\n";
                    }

                    //cout << "センサーx, センサーy" << Sensor[i][0].x << "," << Sensor[i][0].y << "\n";

                    //試行回数分回す
                    for (l = 0; l < loop; l++)
                    {

                        //センサの数分受信電力を求める
                        for (n = 0; n < SensorNum; n++)
                        {

                            //SINRの分布を作るためSmallNum分の受信信号電力が必要
                            for (m = 0; m < SmallNum; m++)
                            {

                                //センサとスモールセル間の距離
                                SensorSmallCellDist[m][n] = TwoPdistance(Sensor[i][n].x, Sensor[i][n].y, SmallCell[m].x, SmallCell[m].y);

                                //受信電力の計算
                                Rx_P = Tx_P[m] - Pathloss(SensorSmallCellDist[m][n], wavelength); //+Shadowing();

                                tmpRSSI = pow(10.0, Rx_P / 10.0);

                                //RayleighFadingValue[o] = tmpRSSI - Rx_P;

                                //cout << "rayleigh" << RayleighFadingValue[o] << "\n";
                                //fprintf(fpA, "Rayleigh, %f\n", RayleighFadingValue[o]);

                                if (i == m)
                                {
                                    OwnPower = tmpRSSI;
                                    tmpOwn = tmpRSSI * RayleighFading();

                                    //自局セルの平均受信電力出し
                                    Grid[j][k].AveOwn[i] += tmpOwn;
                                }
                                else if (SmallCell[i].InterferenceArrivalAreaFlag[m] == 1)
                                {
                                    //dBmをmWに変換
                                    OtherPower_mW += tmpRSSI;

                                    tmp = tmpRSSI * RayleighFading();
                                    tmpOther_mw += tmp;

                                    //SmallCellの番号ごとの干渉電力を集める
                                    Grid[j][k].AveInter[i][m] += tmp;
                                }


                                o++;

                            }

                            InterferenceNoise = OtherPower_mW + Noise_mW;

                            tmpInterNoise = tmpOther_mw + Noise_mW;

                            //fprintf(fpA, "OtherPower_mW, %9.9f,  Noise_mW, %f, InterferenceNoise, %f\n", OtherPower_mW, Noise_mW, InterferenceNoise);
                            //cout << OtherPower_mW << Noise_mW << InterferenceNoise << "\n";

                            //OtherPower_dBm = 10 * log10(InterferenceNoise);

                            //tmpOther_dBm = 10 * log10(tmpInterNoise);

                            //fprintf(fpA, "OwnPower, %f, OtherPower_dBm, %f, OtherPower_noisenasi, %f\n", OwnPower, OtherPower_dBm, 10*log10(OtherPower_mW));
                            //cout << OtherPower_dBm << "\n";

                            SINR = OwnPower / (InterferenceNoise);

                            tmpSINR = tmpOwn / (tmpInterNoise);

                            //cout << tmpInterNoise << "\n";
                            //cout << InterferenceNoise << "\n";
                            //fprintf(fpA, "SINR, %f\n", SINR);

                            Grid[j][k].AveSINR[i] += SINR;
                            Grid[j][k].tmpAveSINR[i] += tmpSINR;


                            //fprintf(fpA, "自局電力, %f, 他局電力, %f, SIR, %f, SumSIR, %f\n", OwnPower, OtherPower_dBm, SIR, Grid[j][k].AveSIR);

                            OwnPower = 0.0;
                            OtherPower_mW = 0.0;
                            OtherPower_dBm = 0.0;
                            SINR = 0.0;

                            tmpOwn = 0.0;
                            tmpOther_mw = 0.0;
                            tmpOther_dBm = 0.0;
                            tmpSINR = 0.0;

                        }

                    }


                    Grid[j][k].AveSINR[i] /= (loop * SensorNum);
                    Grid[j][k].tmpAveSINR[i] /= (loop * SensorNum);
                    Grid[j][k].AveOwn[i] /= (loop * SensorNum);

                    //メッシュの各干渉信号電力の平均値
                    for (l = 0; l < SmallNum; l++)
                    {
                        if (DoNotSharingCellFlag[l] == 0)
                        {
                            Grid[j][k].AveInter[i][l] /= (loop * SensorNum);
                        }
                    }

                    //cout << "fadingnashi" << 10 * log10(Grid[j][k].AveSINR) << "fadingari" << 10 * log10(Grid[j][k].tmpAveSINR) << "\n";

                    //fprintf(fpA, "Gridcoordinates, %d, %d, AVESINR, %f\n", j, k, Grid[j][k].AveSINR);
                    //cout << tmpMinGridx << tmpMinGridy << "\n";
                    if (MinGridAveSINR == 0.0)
                    {
                        MinGridAveSINR = Grid[j][k].AveSINR[i];
                        MinGridx = j;
                        MinGridy = k;
                    }
                    else if (MinGridAveSINR > Grid[j][k].AveSINR[i])
                    {
                        MinGridAveSINR = Grid[j][k].AveSINR[i];
                        MinGridx = j;
                        MinGridy = k;
                    }

                    //fadingなしとありの比較
                    if (tmpMinGridAveSINR == 0.0)
                    {
                        tmpMinGridAveSINR = Grid[j][k].tmpAveSINR[i];
                        tmpMinGridx[i] = j;
                        tmpMinGridy[i] = k;
                    }
                    else if (tmpMinGridAveSINR > Grid[j][k].tmpAveSINR[i])
                    {
                        tmpMinGridAveSINR = Grid[j][k].tmpAveSINR[i];
                        tmpMinGridx[i] = j;
                        tmpMinGridy[i] = k;
                        //cout << tmpMinGridAveSINR << "grid" << j << "grid" << k <<"\n";
                    }
                }
            }
        }


        WorstMeshNumX[i] = tmpMinGridx[i];
        WorstMeshNumY[i] = tmpMinGridy[i];

        CellAveOwn[i] = Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveOwn[i];

        for (j = 0; j < SmallNum; j++)
        {
            if (DoNotSharingCellFlag[j] == 0)
            {
                //各着目セルの最悪メッシュの各干渉受信電力
                CellAveInter[i][j] = Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j];
                InitialCellAveInter[i][j] = Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j];
            }
        }

        //cout << "AveOwn" << 10 * log10(Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveOwn[i]) << "\n";
        writing_file1 << "AveOwn," << 10 * log10(Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveOwn[i]) << endl;


        for (j = 0; j < SmallNum; j++)
        {
            writing_file1 << "AveInter_mW," << Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j] << endl;
            writing_file1 << "AveInter_dBm," << 10.0* log10(Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j]) << endl;

            //cout << "AveInter_dBm" << 10 * log10(Grid[tmpMinGridx][tmpMinGridy].AveInter[j]) << "\n";
        }
        //cout << "MinAveSINR" << 10 * log10(MinGridAveSINR) << "MinGridCoordinates(" << MinGridx << "," << MinGridy << ")" << "\n";
        //cout << "MinAveSINR" << 10 * log10(tmpMinGridAveSINR) << "MinGridCoordinates(" << tmpMinGridx[i] << "," << tmpMinGridy[i] << ")" << "\n";

        writing_file1 << "Non_Fading_MinimumAverageSINR," << 10 * log10(MinGridAveSINR) << ",MinimumGrid_zahyou," << MinGridx << "," << MinGridy << endl;
        writing_file1 << "Fading_MinimumAverageSINR," << 10 * log10(tmpMinGridAveSINR) << ",MinimumGrid_zahyou," << tmpMinGridx[i] << "," << tmpMinGridy[i] << endl;

        //fprintf(fpA, "Fadingあり最小平均自セル受信電力, %f, 他セル1, %f, 他セル2, %f,最小グリッド座標, %d, %d\n", Grid[tmpMinGridx][tmpMinGridy].AveOwn, Grid[tmpMinGridx][tmpMinGridy].AveInter1, Grid[tmpMinGridx][tmpMinGridy].AveInter2, tmpMinGridx, tmpMinGridy);
        //fprintf(fpA, "\n");

        //cout << WorstMeshNumX[i] << "\n";
        //cout << WorstMeshNumY[i] << "\n";
        //cout << "CellAveOwn" << CellAveOwn[i] << "\n";
        writing_file1 << "CellAveOwn," << CellAveOwn[i] << endl;

        for (j = 0; j < SmallNum; j++)
        {
            //cout << "CellAveInter" << 10 * log10(CellAveInter[i][j]) << "\n";
            writing_file1 << "CellAveInter," << 10 * log10(CellAveInter[i][j]) << endl;
        }

        m = 0;
        o = 0;

        MinGridAveSINR = 0.0;
        tmpMinGridAveSINR = 0.0;

        //Gridの初期化
        for (j = 0; j < GridNum; j++)
        {
            for (k = 0; k < GridNum; k++)
            {
                Grid[j][k].flag = 0;

            }
        }
	}


    //Gridの初期化
    for (j = 0; j < SmallNum; j++)
    {
        for (k = 0; k < SmallNum; k++)
        {
            Grid[j][k].flag = 0;

            Grid[j][k].dist[0] = 0.0;
            Grid[j][k].dist[1] = 0.0;
            Grid[j][k].dist[2] = 0.0;
            Grid[j][k].dist[3] = 0.0;
        }
    }

    //電力制御なしで共用できるかの確認
    writing_file1 << endl;
    writing_file1 << "Confirm_spectrum_sharing_without_powercontrol" << endl;

    for (i = 0; i < SmallNum; i++)
    {
        SpectrumSharingImpossible_nonpc[i] = 0;
    }

    for (i = 0; i < SmallNum; i++)
    {
        //cout << "\n";
        //cout << "check point" << i << "\n";
        //cout << "\n";

        //cout << "ADDFlag" << AddFlag[i] << "\n";

        if (DoNotSharingCellFlag[i] == 0)
        {

            //試行回数分回す
            for (j = 0; j < loop * SensorNum; j++)
            {

                OwnPower = CellAveOwn[i] * RayleighFading();

                //干渉電力の総和を求める
                for (l = 0; l < SmallNum; l++)
                {
                    if (DoNotSharingCellFlag[l] == 0)
                    {
                        OtherPower_mW += CellAveInter[i][l] * RayleighFading();
                    }
                }

                AveIsum += OtherPower_mW;

                //(干渉電力の和)と雑音電力の和
                InterferenceNoise = OtherPower_mW + Noise_mW;

                AveSINR += OwnPower / InterferenceNoise;

                SINR = 10 * log10(OwnPower / InterferenceNoise);

                //	cout << "SINR" << SINR << "\n";

                instantaneousSINR[m] = SINR;
                m++;

                if (SINR < DesiredSINR)
                {
                    SpectrumSharingImpossible_nonpc[i]++;
                    if (PoutInitialFlag == 0)
                    {
                        SpectrumSharingImpossible_nonpc_initial[i]++;


                    }
                }

                //諸々初期化
                OwnPower = 0.0;
                OtherPower_mW = 0.0;
                OtherPower_dBm = 0.0;
                SINR = 0.0;

                //cout << "check point4" << "\n";

            }

            m = 0;

            //cout << "check point4" << "\n";

            sort(instantaneousSINR, instantaneousSINR + SIZE_OF_ARRAY(instantaneousSINR));

            //cout << "instantaneousSINR" << instantaneousSINR[1299] << "\n";


            AveIsum /= (loop * SensorNum);
            AveSINR /= (loop * SensorNum);

            //cout << "AveIsum,AveSINR" << AveIsum << "," << AveSINR << "\n";

            writing_file1 << "AverageInterferencePower," << 10 * log10(AveIsum) << endl;
            writing_file1 << "AverageSINR," << 10 * log10(AveSINR) << endl;

            //cout << CDFofTimes[1300] << "\n";



            //SIRthと所望SIRの差分の導出
            DeltaP[i] = instantaneousSINR[OUTAGE] - DesiredSINR;
            InitialDeltaP[i] = instantaneousSINR[OUTAGE] - DesiredSINR;

            writing_file1 << "SINRth,DeltaP," << instantaneousSINR[OUTAGE] << "," <<  DeltaP[i] << endl;
            //cout << "SINRth" << instantaneousSINR[4999] << "DeltaP" << DeltaP[i] << "\n";
            //cout << "AveIsum" << AveIsum << "\n";

            //所望SINRとSINRthを比較しSINRthが大きいと共用可能と判定
            if (DeltaP[i] >= 0.0)
            {
                PossibleSharingCounter_NonPC++;
                //cout << "SharingAchievement!!!" << "\n";
                //cout << "NomPCCellNumber" << i << "\n";

                if (PossibleSharingCounter_NonPC == SmallNum)
                {
                    cout << "Sharing_Achievement(NonPowerControl)!!!" << endl;
                    writing_file1 << "\n**************************************************************************************************************************************************************************" << endl;
                    writing_file1 << "Sharing_Achievement(NonPowerControl)!!!" << endl;
                    writing_file1 << "**************************************************************************************************************************************************************************\n" << endl;
                }

                ////電力制御なしで共用できた数
                //if (NonPowerControlCounter == 0)
                //{
                //	PossibleSharingCounter_NonPC++;
                //}
                //cout << "SharingAchievement!!!" << "\n";
            }
        }


        AveIsum = 0.0;
        AveSINR = 0.0;
    }

    for (i = 0; i < SmallNum; i++)
    {
        writing_file1 << SpectrumSharingImpossible_nonpc[i] << endl;
        writing_file1 << "outage," << (double)SpectrumSharingImpossible_nonpc[i] / (double(loop) * double(SensorNum)) << endl;
    }

    if (PossibleSharingCounter_NonPC != SmallNum)
    {
        cout << "Sharing_Failure(NonPowerControl)!!!" << endl;
        writing_file1 << "\n**************************************************************************************************************************************************************************" << endl;
        writing_file1 << "Sharing_Failure(NonPowerControl)!!!" << endl;
        writing_file1 << "**************************************************************************************************************************************************************************\n" << endl;
    }
/*
    for (i = 0; i < SmallNum; i++)
    {
        pout_nonpc[i] = (double)SpectrumSharingImpossible_nonpc[i] / (double)(loop * SensorNum);

        writing_file1 << "pout_nonpc[i]" << pout_nonpc[i] << endl;
    }
*/

    //初期のpoutは1回だけとる
    PoutInitialFlag = 1;

    //電力制御なしで共用できた数
    writing_file1 << "PossibleSharingCounter_NonPC," << PossibleSharingCounter_NonPC << endl;

    for (i = 0; i < SmallNum; i++)
    {
        SmallCell[i].AllSpectrumBandFlag = 0;
    }

    //共用できなかった場合，電力制御か帯域分割を実行
    if (PossibleSharingCounter_NonPC != SmallNum)
    {
        writing_file1 << endl;
        writing_file1 << "\n**************************************************************************************************************************************************************************" << endl;
        writing_file1 << "Execute_powercontrol_and_banddivision" << endl;
        writing_file1 << "**************************************************************************************************************************************************************************\n" << endl;

        //smallcellに干渉がないため，使用可能周波数の全帯域を割り振る
        for (i = 0; i < SmallNum; i++)
        {
            for (j = 0; j < SmallNum; j++)
            {
                writing_file1 << "AveInter," << Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j] << endl;

                SumAveInter += Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j];

            }
            if (SumAveInter == 0.0)
            {
                //全帯域割り当て可能かのフラッグ立て
                SmallCell[i].AllSpectrumBandFlag = 1;
                SmallCell[i].AlreadyAllocatedBandFlag = 1;
            }

            writing_file1 << "smallcellnum," << i << ",SmallCell[i]AllSpectrumBandFlag," << SmallCell[i].AllSpectrumBandFlag << endl;
            writing_file1 << "smallcellnum," << i << ",SmallCell[i]AlreadyAllocatedBandFlag," << SmallCell[i].AlreadyAllocatedBandFlag << endl;
            SumAveInter = 0.0;
        }

        //先に電力制御をした場合のスループットを導出して帯域分割の手法と比較する
		//電力制御なしで共用できたかを判定
		//できなかった場合に電力制御を行っていく
		if (PossibleSharingCounter_NonPC != (SmallNum - DoNotSharingCellNum))
		{

			//PossibleSharingCounterの初期化
			PossibleSharingCounter = 0;


			////除外したセル以外で共用する
			//if (DoNotSharingCellFlag[i] == 0)
			//{

			//総干渉量が大きいセルから電力制御
			//半径が最小値になるまで繰り返す
			//すべてのセルが制御するまで無限ループ
            pc_count = 1;
            writing_file1 << endl << "Execute_powercontrol(while)" << endl;

            cout << "while no naka" << endl;

			while (1)
			{
                writing_file1 << "\n************************************************************************************************************************" << endl;
                writing_file1 << "PowerControl_count" << pc_count << endl;
                writing_file1 << "************************************************************************************************************************\n" << endl;
                pc_count += 1;
				//総干渉量の導出
				for (i = 0; i < SmallNum; i++)
				{
                    for (j = 0; j < SmallNum; j++)
                    {
                        TotalCellAveInter[j] += CellAveInter[i][j];
                    }
				}
                //他のセルから受ける総干渉量
				for (i = 0; i < SmallNum; i++)
				{
                    writing_file1 << "SmallCellNum," << i << ",OtherCell_karaukeru_soukannsyouryou," << TotalCellAveInter[i] << endl;
				}


				//総干渉量を比較し最大のセルを導出
				for (i = 0; i < SmallNum; i++)
				{

					//
					//cout << "\n";
					//cout << "check point" << i << "\n";
					//cout << "\n";
					//cout << "ADDFlag" << AddFlag[i] << "\n";


					//総干渉量を比較しセル半径を下げるセルの選出
					//電力制御していないかつ除外していないセルにおいて総干渉量を調べる
					if (SmallCell[i].PowerControlDoneFlag == 0)
					{
						if (ComparisonTotalCellAveInter == 0.0)
						{
							ComparisonTotalCellAveInter = TotalCellAveInter[i];
							MaxTotalCellAveInterNum = i;
						}
						else if (ComparisonTotalCellAveInter < TotalCellAveInter[i])
						{
							ComparisonTotalCellAveInter = TotalCellAveInter[i];
							MaxTotalCellAveInterNum = i;
						}
						//cout << "MaxTotalCellAveInterNum" << MaxTotalCellAveInterNum << "\n";
						//cout << "ComparisonTotalCellAveInter" << ComparisonTotalCellAveInter << "\n";

                        writing_file1 << "MaxTotalCellAveInterNum," << MaxTotalCellAveInterNum << endl;
                        writing_file1 << "ComparisonTotalCellAveInter," << ComparisonTotalCellAveInter << endl;
					}

				}
				//もし初めて電力制御するときの最大総干渉量を持つセル番号は除外するときに必要なため保持しておく
				if (FirstTotalCellAveInterFlag == 0)
				{
					FirstTotalCellAveInterFlag = 1;
					FirstTotalCellAveInterNum = MaxTotalCellAveInterNum;

					//cout << "FirstTotalCellAveInterNum" << FirstTotalCellAveInterNum << "\n";
                    writing_file1 << "FirstTotalCellAveInterNum," << FirstTotalCellAveInterNum << endl;
				}

                //cout << "unchi" << endl;
				//最も総干渉量が多かったスモールセルの電力を下げる
				if (SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag == 0)
				{
					//総干渉量が全セルの中最大で電力制御したというflagを立てる

					SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag = 1;
					//cout << "PowerControlDoneFlag" << SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag << "\n";
                    writing_file1 << "PowerControlDoneFlag," << SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag << endl;

    				//対象のセルの送信電力をΔ分下げる
					//cout << "Tx_P[MaxTotalCellAveInterNum]" << Tx_P[MaxTotalCellAveInterNum] << "\n";
                    writing_file1 << "Tx_P[MaxTotalCellAveInterNum]," << Tx_P[MaxTotalCellAveInterNum] << endl;

					//もし総干渉量がゼロのセルが選ばれていたらスルーする
					if (ComparisonTotalCellAveInter != 0.0)
					{
						beforeTx_P = Tx_P[MaxTotalCellAveInterNum];
						Tx_P[MaxTotalCellAveInterNum] -= DeltaTx_P;
						//cout << "seigyogoTx_P[MaxTotalCellAveInterNum]" << Tx_P[MaxTotalCellAveInterNum] << "\n";
                        writing_file1 << "seigyogoTx_P[MaxTotalCellAveInterNum]," << Tx_P[MaxTotalCellAveInterNum] << endl;
					}
					else
					{
                        writing_file1 << "soukannsyouryougazeronoserugaerabareteitanodede_dennryokuseigyohaokonawanai" << endl;
					}


					if (ComparisonTotalCellAveInter != 0.0)
					{
						//対象の干渉セルの送信電力の再設計
						//対象の干渉セルの半径の再設計
						beforeCellRadius = SmallCellRadius[MaxTotalCellAveInterNum];

                        writing_file1 << "beforecellradius," << SmallCellRadius[MaxTotalCellAveInterNum] << endl;

						//beforeTx_P = Tx_P[MaxTotalCellAveInterNum];
						SmallCellRadius[MaxTotalCellAveInterNum] = RefDistance * pow(10.0, (Tx_P[MaxTotalCellAveInterNum] - Noise_dBm - SNR - 20 * log10(4 * PI * RefDistance / wavelength))/ (RefDistance * PathlossExponent));
						//Tx_P[MaxTotalCellAveInterNum] = SNR + Pathloss(SmallCellRadius[MaxTotalCellAveInterNum]) + Noise_dBm;

						//cout << "beforTx_P,Tx_P" << beforeTx_P << Tx_P[MaxTotalCellAveInterNum] << "\n";
						//cout << "beforeCellRadius, SmallCellRadius" << beforeCellRadius << SmallCellRadius[MaxTotalCellAveInterNum] << "\n";

						//fprintf(fpA, "beforTx_P,Tx_P[MaxTotalCellAveInterNum], %20.20f, %20.20f\n", beforeTx_P, Tx_P[MaxTotalCellAveInterNum]);
                        writing_file1 << endl << "Execute PowerControl" << endl;
                        writing_file1 << "BeforeCellRadius," << beforeCellRadius << ",SmallCellRadius[MaxTotalCellAveInterNum]," << SmallCellRadius[MaxTotalCellAveInterNum] << endl;
                        writing_file1 << "InitialTx_P," << InitialTx_P << ",Tx_P[MaxTotalCellAveInterNum]," << Tx_P[MaxTotalCellAveInterNum] << ",10 * log10(CellAveOwn[MaxTotalCellAveInterNum])," << 10 * log10(CellAveOwn[MaxTotalCellAveInterNum]) << ",CellAveOwn[MaxTotalCellAveInterNum]," << CellAveOwn[MaxTotalCellAveInterNum] << endl;

						//初期送信電力と設計後の送信電力の差だけ平均電力を下げる
						for (i = 0; i < SmallNum; i++)
						{
							if (CellAveInter[i][MaxTotalCellAveInterNum] != 0.0)
							{
								CellAveInter[i][MaxTotalCellAveInterNum] = pow(10.0, (10 * log10(CellAveInter[i][MaxTotalCellAveInterNum]) - (beforeTx_P - Tx_P[MaxTotalCellAveInterNum])) / 10.0);

								//cout << "CellAveInter[" << i << "][" << MaxTotalCellAveInterNum << "]" << CellAveInter[i][MaxTotalCellAveInterNum] << "\n";
							}
						}
					}

					NumofPowerControlDoneFlag++;

					//cout << "NumofPowerControlDoneFlag" << NumofPowerControlDoneFlag << "\n";
                    writing_file1 << "NumofPowerControlDoneFlag," << NumofPowerControlDoneFlag << endl;

					//電力制御したというflagのリセット

				}

				//縮小したセルの最悪メッシュを距離と平均SINRより変更しそこの平均自局信号電力と各平均干渉電力を求める
				if ((ComparisonTotalCellAveInter != 0.0) && (SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag == 1))
				{

					//初期化
					tmpMinGridAveSINR = 0.0;

					for (j = 0; j < GridNum; j++)
					{
						for (k = 0; k < GridNum; k++)
						{
							//初期化
							Grid[j][k].flag = 0;


							//四隅分の距離を導出
							Grid[j][k].dist[0] = TwoPdistance(SmallCell[MaxTotalCellAveInterNum].x, SmallCell[MaxTotalCellAveInterNum].y, j  * (MapSize / GridNum), k * (MapSize / GridNum));
							Grid[j][k].dist[1] = TwoPdistance(SmallCell[MaxTotalCellAveInterNum].x, SmallCell[MaxTotalCellAveInterNum].y, (j + 1) * (MapSize / GridNum), k * (MapSize / GridNum));
							Grid[j][k].dist[2] = TwoPdistance(SmallCell[MaxTotalCellAveInterNum].x, SmallCell[MaxTotalCellAveInterNum].y, j * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));
							Grid[j][k].dist[3] = TwoPdistance(SmallCell[MaxTotalCellAveInterNum].x, SmallCell[MaxTotalCellAveInterNum].y, (j + 1) * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));


							//四隅の点と中心点間距離，スモールセルの半径を比較，半径以下だったら旗を立てる
							if ((SmallCellRadius[MaxTotalCellAveInterNum] >= Grid[j][k].dist[0]) || (SmallCellRadius[i] >= Grid[j][k].dist[1]) || (SmallCellRadius[i] >= Grid[j][k].dist[2]) || (SmallCellRadius[i] >= Grid[j][k].dist[3]))
							{
								Grid[j][k].flag = 1;
							}


							//旗が1となる（最悪グリッドとなりうる可能性がある）グリッドでセンサをばらまき試行を重ね，SINRの平均が最も低いグリッドを探す

							//条件式（旗が1ならば実行）
							if (Grid[j][k].flag == 1)
							{

								//平均SINRを比較し最悪メッシュの導出
								if (tmpMinGridAveSINR == 0.0)
								{
									tmpMinGridAveSINR = Grid[j][k].tmpAveSINR[MaxTotalCellAveInterNum];
									tmpMinGridx[MaxTotalCellAveInterNum] = j;
									tmpMinGridy[MaxTotalCellAveInterNum] = k;
								}
								else if (tmpMinGridAveSINR > Grid[j][k].tmpAveSINR[MaxTotalCellAveInterNum])
								{
									tmpMinGridAveSINR = Grid[j][k].tmpAveSINR[MaxTotalCellAveInterNum];
									tmpMinGridx[MaxTotalCellAveInterNum] = j;
									tmpMinGridy[MaxTotalCellAveInterNum] = k;
								}
							}
						}
					}

					//cout << "worstmesh_change" << tmpMinGridx[MaxTotalCellAveInterNum] << "," << tmpMinGridy[MaxTotalCellAveInterNum] << "\n";
                    writing_file1 << "worstmesh_change," <<tmpMinGridx[MaxTotalCellAveInterNum] << "," << tmpMinGridy[MaxTotalCellAveInterNum] << endl;

					//最悪メッシュの平均自局受信電力と平均他局信号電力を格納しなおす

					//平均受信電力にセルを縮小する前とした後の電力の差分を引く
					CellAveOwn[MaxTotalCellAveInterNum] = pow(10.0, (10 * log10(Grid[tmpMinGridx[MaxTotalCellAveInterNum]][tmpMinGridy[MaxTotalCellAveInterNum]].AveOwn[MaxTotalCellAveInterNum]) - (InitialTx_P - Tx_P[MaxTotalCellAveInterNum])) / 10.0);
					//cout << "CellAveOwn[MaxTotalCellAveInterNum]" << CellAveOwn[MaxTotalCellAveInterNum] << "\n";

                    writing_file1 << "CellAveOwn[MaxTotalCellAveInterNum]," << CellAveOwn[MaxTotalCellAveInterNum] << endl;

					for (j = 0; j < SmallNum; j++)
					{
						if ((MaxTotalCellAveInterNum != j) && (Grid[tmpMinGridx[MaxTotalCellAveInterNum]][tmpMinGridy[MaxTotalCellAveInterNum]].AveInter[MaxTotalCellAveInterNum][j] != 0.0))
						{
							CellAveInter[MaxTotalCellAveInterNum][j] = pow(10.0, (10 * log10(Grid[tmpMinGridx[MaxTotalCellAveInterNum]][tmpMinGridy[MaxTotalCellAveInterNum]].AveInter[MaxTotalCellAveInterNum][j]) - (InitialTx_P - Tx_P[j])) / 10.0);

							//cout << "CellAveInter[MaxTotalCellAveInterNum][j] " << CellAveInter[MaxTotalCellAveInterNum][j] << "\n";
                            writing_file1 << "CellAveInter[MaxTotalCellAveInterNum][j]," << CellAveInter[MaxTotalCellAveInterNum][j] << endl;
						}

					}
				}

				for (i = 0; i < SmallNum; i++)
				{
					SpectrumSharingImpossible_pc[i] = 0;
				}
                writing_file1 << "OtherPower_mW," << OtherPower_mW << ",AveIsum," << AveIsum << endl;

				for (i = 0; i < SmallNum; i++)
				{
                    writing_file1 << endl << "SmallNum," << i << endl;
                    writing_file1 << "CellAveOwn[i]," << CellAveOwn[i] << endl;


                    //試行回数分の瞬時受信電力を求める（メッシュ内の受信電力平均値にフェージングによる瞬時変動を加えたもの）
                    for (j = 0; j < loop * SensorNum; j++)
                    {
                        OwnPower = CellAveOwn[i] * RayleighFading();
                        OwnPower_nofade = CellAveOwn[i];

                        //干渉電力の総和を求める
                        for (l = 0; l < SmallNum; l++)
                        {
                            OtherPower_mW += CellAveInter[i][l] * RayleighFading();
                            CellAveInter_nofade += CellAveInter[i][l];
                        }

                        AveIsum += OtherPower_mW;
                        AveSum_nofade += CellAveInter_nofade;

                        //(干渉電力の和)と雑音電力の和
                        InterferenceNoise = OtherPower_mW + Noise_mW;
                        InterferenceNoise_nofade = CellAveInter_nofade + Noise_mW;

                        AveSINR += OwnPower / InterferenceNoise;
                        AveSINR_nofade += OwnPower_nofade / InterferenceNoise_nofade;

                        SINR = 10 * log10(OwnPower / InterferenceNoise);

                        //	cout << "SINR" << SINR << "\n";

                        instantaneousSINR[m] = SINR;
                        m++;

                        if (SINR < DesiredSINR)
                        {
                            SpectrumSharingImpossible_pc[i]++;
                        }

                        //諸々初期化
                        OwnPower = 0.0;
                        OwnPower_nofade = 0.0;
                        OtherPower_mW = 0.0;
                        CellAveInter_nofade = 0.0;
                        OtherPower_dBm = 0.0;
                        SINR = 0.0;
                    }

                    m = 0;

                    //cout << "check point4" << "\n";

                    sort(instantaneousSINR, instantaneousSINR + SIZE_OF_ARRAY(instantaneousSINR));

                    //cout << "instantaneousSINR" << instantaneousSINR[1299] << "\n";

                    AveIsum /= (loop * SensorNum);
                    AveSum_nofade /= (loop * SensorNum);
                    AveSINR /= (loop * SensorNum);
                    AveSINR_nofade /= (loop * SensorNum);

                    writing_file1 << "AveIsum," << AveIsum << ",AveSINR," << AveSINR << "," << 10 * log10(AveSINR) << endl;
                    writing_file1 << "AveSum_nofade," << AveSum_nofade << ",AveSINR_nofade," << AveSINR_nofade << endl;
                    //cout << "AveIsum,AveSINR" << AveIsum << "," << AveSINR << "\n";

                    //writing_file1 << "平均SINR," << 10 * log10(AveSINR) << endl;
                    //cout << CDFofTimes[1300] << "\n";



                    //SIRthと所望SIRの差分の導出
                    DeltaP[i] = instantaneousSINR[OUTAGE] - DesiredSINR;
                    InitialDeltaP[i] = instantaneousSINR[OUTAGE] - DesiredSINR;

                    writing_file1 << "SINRth," << instantaneousSINR[OUTAGE] << ",DeltaP," << DeltaP[i] << endl;
                    //cout << "SINRth" << instantaneousSINR[4999] << "DeltaP" << DeltaP[i] << "\n";
                    //cout << "AveIsum" << AveIsum << "\n";


                    //共用可否判定とFlag立て-------------------------------------------------------------------------------------------------------------------------------------------------
                    if (DeltaP[i] >= 0.0)
                    {
                        PossibleSharingCounter++;

                        ////電力制御なしで共用できた数
                        //if (NonPowerControlCounter == 0)
                        //{
                        //	PossibleSharingCounter_NonPC++;
                        //}
                        //cout << "SharingAchievement!!!" << "\n";

                    }
                    AveIsum = 0.0;
                    AveSum_nofade = 0.0;
					AveSINR = 0.0;
                    AveSINR_nofade = 0.0;
				}


                writing_file1 << endl;
				//通信路容量の導出
				//cmin = log(1 + (pow(10.0,-75.0 / 10.0)) / pow(10, Noise_dBm / 10));
				//c = log(1 + (pow(10.0, Tx_P[] / 10.0)) / pow(10, Noise_dBm / 10));


                //すべてのセルで電力制御を行った場合でもスループットが最小値より大きく，さらに一から電力制御を繰り返す場合

                //スループットの計算
                //cf = (BAND_MAX * pow(10.0, 9.0) + BAND_MIN * pow(10.0, 9.0)) / 2.0;
                writing_file1 << "wavelength," << wavelength << ",cf," << cf << endl;

                //セル端のスループットの計算を行う
                for (k= 0; k < SmallNum; k++)
                {
                    rx_edge[k] = Tx_P[k] - Pathloss(InitialSmallCellRadius, wavelength);
                    writing_file1 << "Tx_P[k]," << Tx_P[k] << ",SmallCellNum," << k << ",rx_edge," << rx_edge[k] << ",InitialSmallCellRadius," << InitialSmallCellRadius << ",SmallCellUseBand," << SmallCell[k].UseBand << endl;
                    writing_file1 << "rx_edge[k]," << rx_edge[k] << endl;

                    snr_edge[k] = pow(10.0, rx_edge[k] / 10.0) / pow(10.0, NOISE_DB / 10.0);

                    writing_file1 << "SmallCellNum," << k << ",snr_edge," << snr_edge[k] << endl;

                    SmallCell[k].Throughput = (INITIAL_BANDWIDTH) * log2(1 + snr_edge[k]);

                    if (pc_count != 1)
                    {
                        SmallCell[k].PowerControlThroughput = SmallCell[k].Throughput;
                    }

                    writing_file1 << "SmallCellNum," << k << ",Throughput," << SmallCell[k].Throughput << ",Throughput[Mbps]," << SmallCell[k].Throughput / pow(10.0, 6.0) << endl;
                    c_sum += SmallCell[k].Throughput;
                }

                writing_file1 << "ThroughputSum," << c_sum << "," << c_sum / pow(10.0, 6.0) << endl;
                writing_file1 << "AverageThroughput," << c_sum / pow(10.0, 6.0) / SmallNum << endl;

				//すべて共用できていたらbreak(抜く必要あり)
				if (PossibleSharingCounter == SmallNum)
				{
					SharingAchievementFlag = 1;
					break;
				}

                for (k = 0; k < SmallNum; k++)
                {
                    //cout << SmallCell[k].Throughput / pow(10.0, 6.0) << endl;
                    if ((SmallCell[k].Throughput / pow(10.0, 6.0)) < 10)
                    {
                        //cout << "lower than 10Mpbs" << endl;
                        writing_file1 << "SmallCell[k].Throughput_is_lower_than_10Mbps," << SmallCell[k].Throughput / pow(10.0, 6.0) << endl;
                        breakflag = 1;
                        break;
                    }
                }

                if (breakflag == 1)
                {
                    break;
                }

                if (NumofPowerControlDoneFlag == SmallNum)
                {
					for (i = 0; i < SmallNum; i++)
					{
						SmallCell[i].PowerControlDoneFlag = 0;
					}
					NumofPowerControlDoneFlag = 0;
                }

				//諸々初期化
				ComparisonTotalCellAveInter = 0.0;
				MaxTotalCellAveInterNum = 0;
				PossibleSharingCounter = 0;
                c_sum = 0.0;

				for (i = 0; i < SmallNum; i++)
				{
					TotalCellAveInter[i] = 0.0;
                    rx_edge[i] = 0.0;
                    snr_edge[i] = 0.0;
                    SmallCell[i].Throughput = 0.0;

				}
			}

			//諸々初期化
			NumofPowerControlDoneFlag = 0;

			FirstTotalCellAveInterFlag = 0;


			//諸々リセット
			for (i = 0; i < SmallNum; i++)
			{
				for (j = 0; j < SmallNum; j++)
				{
					CellAveInter[i][j] = InitialCellAveInter[i][j];
				}
			}
        }

        cout << "whiel no soto" << endl;

		//諸々初期化
		PossibleSharingCounter_NonPC = 0;


        //cout << "SharingAchievement" << "\n";
        writing_file1 << "\n**************************************************************************************************************************************************************************" << endl;
        writing_file1 << "SharingAchievement(PowerControlOnly)!!!" << endl;
        writing_file1 << "**************************************************************************************************************************************************************************\n" << endl;
        //cout << "\n";

        for (i = 0; i < SmallNum; i++)
        {
            //pout_nonpcとpout_nonpc_initialは多分同じで，電力制御しない場合のアウテージ確率
            pout_nonpc[i] = (double)SpectrumSharingImpossible_nonpc[i] / (double)(loop * SensorNum);
            //poutは電力制御をしている時の最悪メッシュのアウテージ確率
            pout[i] = (double)SpectrumSharingImpossible_pc[i] / (double)(loop * SensorNum);
            pout_nonpc_initial[i] = (double)SpectrumSharingImpossible_nonpc_initial[i] / (double)(loop * SensorNum);

            writing_file1 << "pout_nonpc_initial[i]," << pout_nonpc_initial[i] << ",pout_nonpc[i]," << pout_nonpc[i] << ",pout[i]," << pout[i] << ",DoNotSharingCellFlag[i]," << DoNotSharingCellFlag[i] << endl;
        }

        for (i = 0; i < SmallNum; i++)
        {
            //アウテージ確率の初期化
            pout_nonpc[i] = 0.0;
            pout[i] = 0.0;
            pout_nonpc_initial[i] = 0.0;
            SpectrumSharingImpossible_nonpc[i] = 0;
            SpectrumSharingImpossible_pc[i] = 0;
            SpectrumSharingImpossible_nonpc_initial[i] = 0;

            writing_file1 << "pout_nonpc_initial[i]," << pout_nonpc_initial[i] << ",pout_nonpc[i]," << pout_nonpc[i] << ",pout[i]," << pout[i] << ",DoNotSharingCellFlag[i]," << DoNotSharingCellFlag[i] << endl;
        }

        //cout << "SharingAchievement" << "\n";
        writing_file1 << "\n**************************************************************************************************************************************************************************" << endl;
        writing_file1 << "Execute SpectrumBandDivideAllocation and PowerControl!!!" << endl;
        writing_file1 << "**************************************************************************************************************************************************************************\n" << endl;
        //cout << "\n";


        //new program sakuseit now
        //ループして分割数を増やす構造にして，分割数の少ない状態から電力制御のみのスループットと比較and共用可否の判定をし続ける
        for (DividedNum = DivideNumMin; DividedNum < DivideNumMax; DividedNum++)
        {

            //各バンドの中心周波数の作成
            for (i = 0; i < DividedNum; i++)
            {
                //f1 = ((BAND_MAX + BAND_MIN) / 2.0 + BAND_MIN) / 2.0;
                //f2 = ((BAND_MAX + BAND_MIN) / 2.0 + BAND_MAX) / 2.0;
                f[0] = ((BAND_MAX + BAND_MIN) / 2.0 + BAND_MIN) / 2.0;
                f[1] = ((BAND_MAX + BAND_MIN) / 2.0 + BAND_MAX) / 2.0;
            }

            cout << "f1:" << f1 << endl;
            cout << "f2:" << f2 << endl;
            cout << "f[0]:" << f[0] << endl;
            cout << "f[1]:" << f[1] << endl;

            for (i  = 0; i < DividedNum; i++)
            {
                for (j = 0; j < DividedNum; j++)
                {
                    for (k = 0; k < DividedNum; k++)
                    {
                        //全て同じバンドは割り当てない
                        if (i == j == k)
                        {
                            writing_file1 << "allcell is same band" << endl;
                        }
                        else
                        {
                            writing_file1 << i << "," << j << "," << k << endl;


                            //各セルにバンドを適切に割り当て
                            //SmallCell[0].UseBand = i *
                            //SmallCell[1].UseBand =
                            //SmallCell[2].UseBand =

                            if (c_sum_dp > c_sum)
                            {
                                break;
                            }
                        }
                    }
                }

            }

        }



        //全てのセルとバンドの組み合わせてスループットを導出
        //for (i < 0; i < SmallNum * DividedNum; i++)
        //{
            //全帯域で割り当てられなかったセルの帯域割り当て
            //i番目のセルがf1で，他のセルがf2だった時に共用の可否・スループットがどうなるか

            for (z = 0; z < SmallNum; z++)
            {
                writing_file1 << "\n\nThroughput_\n\n" << endl;

                if ((SmallCell[z].AllSpectrumBandFlag == 0) && (SmallCell[z].AlreadyAllocatedBandFlag == 0))
                {

                    writing_file1 << "SmallCellNum," << z << ",target_allocation_cell" << endl;

                    SmallCell[z].AlreadyAllocatedBandFlag = 1;

                    SmallCell[z].UseBand = f[0] * pow(10.0, 9.0);

                    numf[0] += 1;
                    //numf1 += 1;

                    for (k = 0; k < SmallNum; k++)
                    {
                        if ((z != k) && (SmallCell[k].AlreadyAllocatedBandFlag == 0))
                        {
                            SmallCell[k].UseBand = f[1] * pow(10.0, 9.0);
                            numf[1] += 1;
                            //numf2 += 1;
                        }
                        else if((z != k) && (SmallCell[k].AllSpectrumBandFlag == 1))
                        {
                            SmallCell[k].UseBand = cf;
                            numcf += 1;
                        }
                    }

                    for (k = 0; k < SmallNum; k++)
                    {
                        writing_file1 << "SmallCellNum," << k << ",UseBand," << SmallCell[k].UseBand << endl;
                    }

                    writing_file1 << "numf[1]," << SmallNum - numf[1] - numcf << ",numf2," << numf[1] << ",numcf," << numcf << endl;

                    //各バンドでセルが一つ以下だったら共用可否の判定はいらない
                    if ((numf[0] <= 1) && (numf[1] <= 1))
                    {
                        writing_file1 << "\n\nkakubandode_cell_ga_hitotuika\n\n" << endl;

                        //スループットの計算を行う
                        for (k= 0; k < SmallNum; k++)
                        {
                            wavelength = 3.0 * pow(10.0, 8.0) / SmallCell[k].UseBand;
                            writing_file1 << "wavelength," << wavelength << ",SmallCell[k].UseBand," << SmallCell[k].UseBand << endl;
                            Tx_P[k] = SNR + Pathloss(InitialSmallCellRadius, wavelength) + Noise_dBm;
                            rx_edge[k] = Tx_P[k] - Pathloss(InitialSmallCellRadius, wavelength);
                            writing_file1 << "Tx_P[k]," << Tx_P[k] << ",SmallCellNum," << k << ",rx_edge," << rx_edge[k] << ",InitialSmallCellRadius," << InitialSmallCellRadius << ",SmallCellUseBand," << SmallCell[k].UseBand << endl;
                            writing_file1 << "SmallCellNum," << k << ",rx_edge," << rx_edge[k] << ",InitialSmallCellRadius," << InitialSmallCellRadius << ",SmallCellUseBand," << SmallCell[k].UseBand << endl;
                            //cout << "rx_edge[k]" << rx_edge[k] << endl;

                            snr_edge[k] = pow(10.0, rx_edge[k] / 10.0) / pow(10.0, NOISE_DB / 10.0);

                            writing_file1 << "SmallCellNum," << k << ",snr_edge," << snr_edge[k] << endl;
                            //cout << "snr_edge[k]" << snr_edge[k] << endl;

                            if (SmallCell[k].AllSpectrumBandFlag == 1)
                            {
                                SmallCell[k].Throughput = (INITIAL_BANDWIDTH) * log2(1 + snr_edge[k]);
                                SmallCell[k].Throughput_GB = (INITIAL_BANDWIDTH) * log2(1 + snr_edge[k]);

                            }
                            else
                            {
                                SmallCell[k].Throughput = (INITIAL_BANDWIDTH / DividedNum) * log2(1 + snr_edge[k]);
                                SmallCell[k].Throughput_GB = ((INITIAL_BANDWIDTH - (GB * pow(10.0,6.0) * (DividedNum-1))) / DividedNum) * log2(1 + snr_edge[k]);

                                writing_file1 << "GB_BandWidth," << ((INITIAL_BANDWIDTH - (GB * pow(10.0,6.0) * (DividedNum-1))) / DividedNum) << endl;
                            }


                            writing_file1 << "SmallCellNum," << k << ",Throughput," << SmallCell[k].Throughput << endl;
                            writing_file1 << "SmallCellNum," << k << ",GB_Throughput," << SmallCell[k].Throughput_GB << endl;
                            c_sum += SmallCell[k].Throughput / pow(10.0, 6.0);
                            GB_c_sum += SmallCell[k].Throughput_GB / pow(10.0, 6.0);
                        }

                        writing_file1 << "ThroughputSum," << c_sum << ",AverageThroughput," << c_sum / SmallNum << endl;
                        writing_file1 << "GB_ThroughputSum," << GB_c_sum << ",AverageThroughput," << GB_c_sum / SmallNum << endl;
                        c_sum = 0.0;
                        GB_c_sum = 0.0;
                    }
                    else
                    {
                        writing_file1 << "\n\nkakubandode_cell_ga_hutatuizyou\n\n" << endl;

                        writing_file1 << "f1f2de_powercontrol\n\n" << endl;

                        for (y = 0; y < DIVIDED_NUM; y++)
                        {
                            writing_file1 << "f" << f[y] << "\n\n" << endl;


                            if (numf[y] >= 2)
                            {


                                if (SmallCell[z].UseBand != (f[y] * pow(10.0, 9.0)))
                                {
                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        if (SmallCell[i].UseBand != (f[y] * pow(10.0, 9.0)))
                                        {
                                            writing_file1 << "saihaichi,cellx," << SmallCell[i].x << ",celly," << SmallCell[i].y << endl;
                                        }
                                    }

                                    //初期のセル半径決め
                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        if (SmallCell[i].UseBand != (f[y] * pow(10.0, 9.0)))
                                        {
                                            SmallCellRadius[i] = InitialSmallCellRadius;
                                        }
                                    }

                                    //雑音電力
                                    Noise_dBm = -95.0;
                                    Noise_mW = pow(10, Noise_dBm / 10.0);

                                    //cf = (BAND_MAX * pow(10.0, 9.0) + BAND_MIN * pow(10.0, 9.0)) / 2.0;
                                    //cf = 4625000000;
                                    //wavelength = 3.0 * pow(10.0, 8.0) / cf;

                                    //各スモールセルの初期送信電力決め
                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        if (SmallCell[i].UseBand != (f[y] * pow(10.0, 9.0)))
                                        {
                                            Tx_P[i] = SNR + Pathloss(SmallCellRadius[i], wavelength) + Noise_dBm;
                                            SmallCellRadius[i] = RefDistance * pow(10.0, (Tx_P[i] - Noise_dBm - SNR - 20 * log10(4 * PI * RefDistance / wavelength))/ (RefDistance * PathlossExponent));

                                            //初期送信電力をとっておく
                                            InitialTx_P = Tx_P[i];
                                            writing_file1 << "SmallCellRadius," << SmallCellRadius[i] << endl;
                                            writing_file1 << "Tx_P,SmallCellnum," << i << "," << Tx_P[i] << endl;
                                        }
                                    }

                                    //干渉到達エリアの半径決め
                                    InterferenceArrivalAreaRadius = RefDistance * pow(10.0, (InitialTx_P - Noise_dBm - 20.0 * log10(4.0 * PI * RefDistance / wavelength)) / (10.0 * PathlossExponent));
                                    writing_file1 << "InterferenceArrivalAreaRadius," << InterferenceArrivalAreaRadius << endl;

                                    //セルに属するメッシュの判定フラグの初期化
                                    for (i = 0; i < GridNum; i++)
                                    {
                                        for (j = 0; j < GridNum; j++)
                                        {
                                            for (k = 0; k < SmallNum; k++)
                                            {
                                                Grid[i][j].flag = 0;
                                                Grid[i][j].AveSINR[k] = 0.0;
                                                Grid[i][j].tmpAveSINR[k] = 0.0;
                                                Grid[i][j].AveOwn[k] = 0.0;
                                            }
                                        }
                                    }

                                    for (i = 0; i < GridNum; i++)
                                    {
                                        for (j = 0; j < GridNum; j++)
                                        {
                                            for (k = 0; k < SmallNum; k++)
                                            {
                                                for (l = 0; l < SmallNum; l++)
                                                {
                                                    Grid[i][j].AveInter[l][k] = 0.0;
                                                }
                                            }
                                        }
                                    }

                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        if (SmallCell[i].UseBand != (f[y] * pow(10.0, 9.0)))
                                        {
                                            CellAveOwn[i] = 0.0;
                                            WorstMeshCenterInterCellDist[i] = 0.0;
                                            SmallCell[i].PowerControlDoneFlag = 0;

                                            for (j = 0; j < SmallNum; j++)
                                            {
                                                if (SmallCell[j].UseBand != (f[y] * pow(10.0, 9.0)))
                                                {
                                                    CellAveInter[i][j] = 0.0;
                                                    InitialCellAveInter[i][j] = 0.0;

                                                    SmallCell[i].InterferenceArrivalAreaFlag[j] = 0;
                                                }
                                            }
                                        }
                                    }

                                    //セルごとに最悪メッシュの導出 それらの情報を保持---------------------------------------------------------------------------------------------------------
                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        if (SmallCell[i].UseBand != (f[y] * pow(10.0, 9.0)))
                                        {
                                            writing_file1 << "\nSmallNum," << i << endl;

                                            //干渉信号到達エリア内に着目スモールセル以外のセルがあるか判定
                                            for (j = 0; j < SmallNum; j++)
                                            {
                                                if (SmallCell[j].UseBand != (f[y] * pow(10.0, 9.0)))
                                                {
                                                    SmallCell[i].dist[j] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, SmallCell[j].x, SmallCell[j].y);

                                                    writing_file1 << "inter-CellDist," << SmallCell[i].dist[j] << endl;

                                                    //cout << "Radius" << InterferenceArrivalAreaRadius - SmallCellRadius[i] << "\n";
                                                    //着目セルと着目セル以外のセルの距離と干渉信号到達エリアの半径を比較
                                                    if ((SmallCell[i].dist[j] < (InterferenceArrivalAreaRadius + Margin - SmallCellRadius[i])) && (i != j))
                                                    {
                                                        SmallCell[i].InterferenceArrivalAreaFlag[j] = 1;

                                                        //cout << "InterferenceArrivalArea" << j << "\n";
                                                        writing_file1 << "InterferenceSignalArea_naino_smallcell," << j << ",SmallCell[i].InterferenceArrivalAreaFlag[j]," << SmallCell[i].InterferenceArrivalAreaFlag[j] << endl;
                                                        //cout << SmallCell[i].InterferenceArrivalAreaFlag[j] << "\n";
                                                    }
                                                }
                                            }


                                            //着目スモールセルの最悪点を求める
                                            //「全メッシュの四隅の点とセルの中心点間の距離」，「スモールセルの半径」を比較して半径よりも距離が短い点が1つでもあればスモールセル内にあるメッシュとして旗を立てる
                                            for (j = 0; j < GridNum; j++)
                                            {
                                                for (k = 0; k < GridNum; k++)
                                                {
                                                    //四隅分の距離を導出
                                                    Grid[j][k].dist[0] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, j  * (MapSize / GridNum), k * (MapSize / GridNum));
                                                    Grid[j][k].dist[1] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, (j + 1) * (MapSize / GridNum), k * (MapSize / GridNum));
                                                    Grid[j][k].dist[2] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, j * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));
                                                    Grid[j][k].dist[3] = TwoPdistance(SmallCell[i].x, SmallCell[i].y, (j + 1) * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));
                                                    //fprintf(fpA, "GridNum, %d, %d, SmallCellNum, %d, dist, %f\n", j, k, i, Grid[j][k].dist[0]);
                                                    //fprintf(fpA, "GridNum, %d, %d, SmallCellNum, %d, dist, %f\n", j, k, i, Grid[j][k].dist[1]);
                                                    //fprintf(fpA, "GridNum, %d, %d, SmallCellNum, %d, dist, %f\n", j, k, i, Grid[j][k].dist[2]);
                                                    //fprintf(fpA, "GridNum, %d, %d, SmallCellNum, %d, dist, %f\n", j, k, i, Grid[j][k].dist[3]);


                                                    //四隅の点と中心点間距離，スモールセルの半径を比較，半径以下だったら旗を立てる
                                                    if ((SmallCellRadius[i] >= Grid[j][k].dist[0]) || (SmallCellRadius[i] >= Grid[j][k].dist[1]) || (SmallCellRadius[i] >= Grid[j][k].dist[2]) || (SmallCellRadius[i] >= Grid[j][k].dist[3]))
                                                    {
                                                        Grid[j][k].flag = 1;
                                                        //cout << "Gridcoordinates(" << j << "," << k << ")," << i << "numberSmallCellCenter(" << SmallCell[i].x << "," << SmallCell[i].y << ")distance" << Grid[j][k].dist[0] << "\n";
                                                        //fprintf(fpA, "Gridcoordinates, %d, %d, %d, numberSmallCellCenter, %f, %f, distance, %f\n", j, k, i, SmallCell[i].x, SmallCell[i].y, Grid[j][k].dist[0]);
                                                    }


                                                    //旗が1となる（最悪グリッドとなりうる可能性がある）グリッドでセンサをばらまき試行を重ね，SINRの平均が最も低いグリッドを探す

                                                    //条件式（旗が1ならば実行）
                                                    if (Grid[j][k].flag == 1)
                                                    {
                                                        //センサをSensorNum個だけばらまく
                                                        for (l = 0; l < SensorNum; l++)
                                                        {
                                                            //センサはグリッド座標(j, k)を使って一様乱数でまく
                                                            Sensor[i][l].x = UniformRandom(j * (MapSize / GridNum), (j + 1) * (MapSize / GridNum));
                                                            Sensor[i][l].y = UniformRandom(k * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));
                                                            //fprintf(fpA, "SensorCoordinates, %f, %f\n", Sensor[i][l].x, Sensor[i][l].y);
                                                            //cout << "sensor_x" << Sensor[i][l].x << "sensor_y" << Sensor[i][l].y << "\n";
                                                        }

                                                        //cout << "センサーx, センサーy" << Sensor[i][0].x << "," << Sensor[i][0].y << "\n";

                                                        //試行回数分回す
                                                        for (l = 0; l < loop; l++)
                                                        {

                                                            //センサの数分受信電力を求める
                                                            for (n = 0; n < SensorNum; n++)
                                                            {

                                                                //SINRの分布を作るためSmallNum分の受信信号電力が必要
                                                                for (m = 0; m < SmallNum; m++)
                                                                {
                                                                    if (SmallCell[m].UseBand != (f[y] * pow(10.0, 9.0)))
                                                                    {
                                                                        //センサとスモールセル間の距離
                                                                        SensorSmallCellDist[m][n] = TwoPdistance(Sensor[i][n].x, Sensor[i][n].y, SmallCell[m].x, SmallCell[m].y);

                                                                        //受信電力の計算
                                                                        Rx_P = Tx_P[m] - Pathloss(SensorSmallCellDist[m][n], wavelength); //+Shadowing();

                                                                        tmpRSSI = pow(10.0, Rx_P / 10.0);

                                                                        //RayleighFadingValue[o] = tmpRSSI - Rx_P;

                                                                        //cout << "rayleigh" << RayleighFadingValue[o] << "\n";
                                                                        //fprintf(fpA, "Rayleigh, %f\n", RayleighFadingValue[o]);

                                                                        if (i == m)
                                                                        {
                                                                            OwnPower = tmpRSSI;
                                                                            tmpOwn = tmpRSSI * RayleighFading();

                                                                            //自局セルの平均受信電力出し
                                                                            Grid[j][k].AveOwn[i] += tmpOwn;
                                                                        }
                                                                        else if (SmallCell[i].InterferenceArrivalAreaFlag[m] == 1)
                                                                        {
                                                                            //dBmをmWに変換
                                                                            OtherPower_mW += tmpRSSI;

                                                                            tmp = tmpRSSI * RayleighFading();
                                                                            tmpOther_mw += tmp;

                                                                            //SmallCellの番号ごとの干渉電力を集める
                                                                            Grid[j][k].AveInter[i][m] += tmp;
                                                                        }


                                                                        o++;
                                                                    }
                                                                }

                                                                InterferenceNoise = OtherPower_mW + Noise_mW;

                                                                tmpInterNoise = tmpOther_mw + Noise_mW;

                                                                //fprintf(fpA, "OtherPower_mW, %9.9f,  Noise_mW, %f, InterferenceNoise, %f\n", OtherPower_mW, Noise_mW, InterferenceNoise);
                                                                //cout << OtherPower_mW << Noise_mW << InterferenceNoise << "\n";

                                                                //OtherPower_dBm = 10 * log10(InterferenceNoise);

                                                                //tmpOther_dBm = 10 * log10(tmpInterNoise);

                                                                //fprintf(fpA, "OwnPower, %f, OtherPower_dBm, %f, OtherPower_noisenasi, %f\n", OwnPower, OtherPower_dBm, 10*log10(OtherPower_mW));
                                                                //cout << OtherPower_dBm << "\n";

                                                                SINR = OwnPower / (InterferenceNoise);

                                                                tmpSINR = tmpOwn / (tmpInterNoise);

                                                                //cout << tmpInterNoise << "\n";
                                                                //cout << InterferenceNoise << "\n";
                                                                //fprintf(fpA, "SINR, %f\n", SINR);

                                                                Grid[j][k].AveSINR[i] += SINR;
                                                                Grid[j][k].tmpAveSINR[i] += tmpSINR;


                                                                //fprintf(fpA, "自局電力, %f, 他局電力, %f, SIR, %f, SumSIR, %f\n", OwnPower, OtherPower_dBm, SIR, Grid[j][k].AveSIR);

                                                                OwnPower = 0.0;
                                                                OtherPower_mW = 0.0;
                                                                OtherPower_dBm = 0.0;
                                                                SINR = 0.0;

                                                                tmpOwn = 0.0;
                                                                tmpOther_mw = 0.0;
                                                                tmpOther_dBm = 0.0;
                                                                tmpSINR = 0.0;

                                                            }

                                                        }


                                                        Grid[j][k].AveSINR[i] /= (loop * SensorNum);
                                                        Grid[j][k].tmpAveSINR[i] /= (loop * SensorNum);
                                                        Grid[j][k].AveOwn[i] /= (loop * SensorNum);

                                                        //メッシュの各干渉信号電力の平均値
                                                        for (l = 0; l < SmallNum; l++)
                                                        {
                                                            if (SmallCell[l].UseBand != (f[y] * pow(10.0, 9.0)))
                                                            {
                                                                if (DoNotSharingCellFlag[l] == 0)
                                                                {
                                                                    Grid[j][k].AveInter[i][l] /= (loop * SensorNum);
                                                                }
                                                            }
                                                        }

                                                        //cout << "fadingnashi" << 10 * log10(Grid[j][k].AveSINR) << "fadingari" << 10 * log10(Grid[j][k].tmpAveSINR) << "\n";

                                                        //fprintf(fpA, "Gridcoordinates, %d, %d, AVESINR, %f\n", j, k, Grid[j][k].AveSINR);
                                                        //cout << tmpMinGridx << tmpMinGridy << "\n";
                                                        if (MinGridAveSINR == 0.0)
                                                        {
                                                            MinGridAveSINR = Grid[j][k].AveSINR[i];
                                                            MinGridx = j;
                                                            MinGridy = k;
                                                        }
                                                        else if (MinGridAveSINR > Grid[j][k].AveSINR[i])
                                                        {
                                                            MinGridAveSINR = Grid[j][k].AveSINR[i];
                                                            MinGridx = j;
                                                            MinGridy = k;
                                                        }

                                                        //fadingなしとありの比較
                                                        if (tmpMinGridAveSINR == 0.0)
                                                        {
                                                            tmpMinGridAveSINR = Grid[j][k].tmpAveSINR[i];
                                                            tmpMinGridx[i] = j;
                                                            tmpMinGridy[i] = k;
                                                        }
                                                        else if (tmpMinGridAveSINR > Grid[j][k].tmpAveSINR[i])
                                                        {
                                                            tmpMinGridAveSINR = Grid[j][k].tmpAveSINR[i];
                                                            tmpMinGridx[i] = j;
                                                            tmpMinGridy[i] = k;
                                                            //cout << tmpMinGridAveSINR << "grid" << j << "grid" << k <<"\n";
                                                        }
                                                    }
                                                }
                                            }


                                            WorstMeshNumX[i] = tmpMinGridx[i];
                                            WorstMeshNumY[i] = tmpMinGridy[i];

                                            CellAveOwn[i] = Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveOwn[i];

                                            for (j = 0; j < SmallNum; j++)
                                            {
                                                if (SmallCell[j].UseBand != (f[y] * pow(10.0, 9.0)))
                                                {
                                                    if (DoNotSharingCellFlag[j] == 0)
                                                    {
                                                        //各着目セルの最悪メッシュの各干渉受信電力
                                                        CellAveInter[i][j] = Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j];
                                                        InitialCellAveInter[i][j] = Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j];
                                                    }
                                                }
                                            }

                                            //cout << "AveOwn" << 10 * log10(Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveOwn[i]) << "\n";
                                            writing_file1 << "AveOwn," << 10 * log10(Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveOwn[i]) << endl;


                                            for (j = 0; j < SmallNum; j++)
                                            {
                                                if (SmallCell[j].UseBand != (f[y] * pow(10.0, 9.0)))
                                                {
                                                    writing_file1 << "AveInter_mW," << Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j] << endl;
                                                    writing_file1 << "AveInter_dBm," << 10.0* log10(Grid[tmpMinGridx[i]][tmpMinGridy[i]].AveInter[i][j]) << endl;

                                                    //cout << "AveInter_dBm" << 10 * log10(Grid[tmpMinGridx][tmpMinGridy].AveInter[j]) << "\n";
                                                }
                                            }
                                            //cout << "MinAveSINR" << 10 * log10(MinGridAveSINR) << "MinGridCoordinates(" << MinGridx << "," << MinGridy << ")" << "\n";
                                            //cout << "MinAveSINR" << 10 * log10(tmpMinGridAveSINR) << "MinGridCoordinates(" << tmpMinGridx[i] << "," << tmpMinGridy[i] << ")" << "\n";

                                            writing_file1 << "Non_Fading_MinimumAverageSINR," << 10 * log10(MinGridAveSINR) << ",MinimumGrid_zahyou最小," << MinGridx << "," << MinGridy << endl;
                                            writing_file1 << "Fading_MinimumAverageSINR," << 10 * log10(tmpMinGridAveSINR) << ",MinimumGrid_zahyou最小," << tmpMinGridx[i] << "," << tmpMinGridy[i] << endl;

                                            //fprintf(fpA, "Fadingあり最小平均自セル受信電力, %f, 他セル1, %f, 他セル2, %f,最小グリッド座標, %d, %d\n", Grid[tmpMinGridx][tmpMinGridy].AveOwn, Grid[tmpMinGridx][tmpMinGridy].AveInter1, Grid[tmpMinGridx][tmpMinGridy].AveInter2, tmpMinGridx, tmpMinGridy);
                                            //fprintf(fpA, "\n");

                                            //cout << WorstMeshNumX[i] << "\n";
                                            //cout << WorstMeshNumY[i] << "\n";
                                            //cout << "CellAveOwn" << CellAveOwn[i] << "\n";
                                            writing_file1 << "CellAveOwn," << CellAveOwn[i] << endl;

                                            for (j = 0; j < SmallNum; j++)
                                            {
                                                if (SmallCell[j].UseBand != (f[y] * pow(10.0, 9.0)))
                                                {
                                                    //cout << "CellAveInter" << 10 * log10(CellAveInter[i][j]) << "\n";
                                                    writing_file1 << "CellAveInter," << 10 * log10(CellAveInter[i][j]) << endl;
                                                }
                                            }

                                            m = 0;
                                            o = 0;

                                            MinGridAveSINR = 0.0;
                                            tmpMinGridAveSINR = 0.0;

                                            //Gridの初期化
                                            for (j = 0; j < GridNum; j++)
                                            {
                                                for (k = 0; k < GridNum; k++)
                                                {
                                                    Grid[j][k].flag = 0;

                                                }
                                            }
                                        }
                                    }


                                    //Gridの初期化
                                    for (j = 0; j < SmallNum; j++)
                                    {
                                        for (k = 0; k < SmallNum; k++)
                                        {
                                            Grid[j][k].flag = 0;

                                            Grid[j][k].dist[0] = 0.0;
                                            Grid[j][k].dist[1] = 0.0;
                                            Grid[j][k].dist[2] = 0.0;
                                            Grid[j][k].dist[3] = 0.0;
                                        }
                                    }

                                    //電力制御なしで共用できるかの確認
                                    writing_file1 << endl;
                                    writing_file1 << "Confirm_spectrum_sharing_without_powercontrol" << endl;

                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        if (SmallCell[i].UseBand != (f[y] * pow(10.0, 9.0)))
                                        {
                                            SpectrumSharingImpossible_nonpc[i] = 0;
                                        }
                                    }

                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        if (SmallCell[i].UseBand != (f[y] * pow(10.0, 9.0)))
                                        {
                                            //cout << "\n";
                                            //cout << "check point" << i << "\n";
                                            //cout << "\n";

                                            //cout << "ADDFlag" << AddFlag[i] << "\n";

                                            if (DoNotSharingCellFlag[i] == 0)
                                            {

                                                //試行回数分回す
                                                for (j = 0; j < loop * SensorNum; j++)
                                                {

                                                    OwnPower = CellAveOwn[i] * RayleighFading();

                                                    //干渉電力の総和を求める
                                                    for (l = 0; l < SmallNum; l++)
                                                    {
                                                        if (SmallCell[l].UseBand != (f[y] * pow(10.0, 9.0)))
                                                        {
                                                            if (DoNotSharingCellFlag[l] == 0)
                                                            {
                                                                OtherPower_mW += CellAveInter[i][l] * RayleighFading();
                                                            }
                                                        }
                                                    }

                                                    AveIsum += OtherPower_mW;

                                                    //(干渉電力の和)と雑音電力の和
                                                    InterferenceNoise = OtherPower_mW + Noise_mW;

                                                    AveSINR += OwnPower / InterferenceNoise;

                                                    SINR = 10 * log10(OwnPower / InterferenceNoise);

                                                    //	cout << "SINR" << SINR << "\n";

                                                    instantaneousSINR[m] = SINR;
                                                    m++;

                                                    if (SINR < DesiredSINR)
                                                    {
                                                        SpectrumSharingImpossible_nonpc[i]++;
                                                        if (PoutInitialFlag == 0)
                                                        {
                                                            SpectrumSharingImpossible_nonpc_initial[i]++;


                                                        }
                                                    }

                                                    //諸々初期化
                                                    OwnPower = 0.0;
                                                    OtherPower_mW = 0.0;
                                                    OtherPower_dBm = 0.0;
                                                    SINR = 0.0;

                                                    //cout << "check point4" << "\n";

                                                }

                                                m = 0;

                                                //cout << "check point4" << "\n";

                                                sort(instantaneousSINR, instantaneousSINR + SIZE_OF_ARRAY(instantaneousSINR));

                                                //cout << "instantaneousSINR" << instantaneousSINR[1299] << "\n";


                                                AveIsum /= (loop * SensorNum);
                                                AveSINR /= (loop * SensorNum);

                                                //cout << "AveIsum,AveSINR" << AveIsum << "," << AveSINR << "\n";

                                                writing_file1 << "AverageInterferencePower," << 10 * log10(AveIsum) << endl;
                                                writing_file1 << "AverageSINR," << 10 * log10(AveSINR) << endl;

                                                //cout << CDFofTimes[1300] << "\n";



                                                //SIRthと所望SIRの差分の導出
                                                DeltaP[i] = instantaneousSINR[OUTAGE] - DesiredSINR;
                                                InitialDeltaP[i] = instantaneousSINR[OUTAGE] - DesiredSINR;

                                                writing_file1 << "SINRth,DeltaP," << instantaneousSINR[OUTAGE] << "," <<  DeltaP[i] << endl;
                                                //cout << "SINRth" << instantaneousSINR[4999] << "DeltaP" << DeltaP[i] << "\n";
                                                //cout << "AveIsum" << AveIsum << "\n";

                                                //所望SINRとSINRthを比較しSINRthが大きいと共用可能と判定
                                                if (DeltaP[i] >= 0.0)
                                                {
                                                    PossibleSharingCounter_NonPC++;
                                                    //cout << "SharingAchievement!!!" << "\n";
                                                    //cout << "NomPCCellNumber" << i << "\n";

                                                    writing_file1 << "\n**************************************************************************************************************************************************************************" << endl;
                                                    writing_file1 << "SharingAchievement!!!" << endl;
                                                    writing_file1 << "**************************************************************************************************************************************************************************\n" << endl;
                                                    writing_file1 << "NomPCCellNumber," << i << endl;
                                                    ////電力制御なしで共用できた数
                                                    //if (NonPowerControlCounter == 0)
                                                    //{
                                                    //	PossibleSharingCounter_NonPC++;
                                                    //}
                                                    //cout << "SharingAchievement!!!" << "\n";
                                                }
                                            }
                                            AveIsum = 0.0;
                                            AveSINR = 0.0;
                                        }
                                    }

                                /*
                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        pout_nonpc[i] = (double)SpectrumSharingImpossible_nonpc[i] / (double)(loop * SensorNum);

                                        writing_file1 << "pout_nonpc[i]" << pout_nonpc[i] << endl;
                                    }
                                */
                                /*
                                    //初期のpoutは1回だけとる
                                    PoutInitialFlag = 1;

                                    //電力制御なしで共用できた数
                                    writing_file1 << "PossibleSharingCounter_NonPC," << PossibleSharingCounter_NonPC << endl;
                                */
                                /*
                                    //電力制御（スペクトラムバンド分割した状態で）
                                    while (1)
                                    {
                                        writing_file1 << "\n************************************************************************************************************************" << endl;
                                        writing_file1 << "PowerControl_count" << pc_count << endl;
                                        writing_file1 << "************************************************************************************************************************\n" << endl;
                                        pc_count += 1;
                                        //総干渉量の導出
                                        for (i = 0; i < SmallNum; i++)
                                        {
                                            for (j = 0; j < SmallNum; j++)
                                            {
                                                TotalCellAveInter[j] += CellAveInter[i][j];
                                            }
                                        }
                                        //他のセルから受ける総干渉量
                                        for (i = 0; i < SmallNum; i++)
                                        {
                                            writing_file1 << "SmallCellNum," << i << ",OtherCell_karaukeru_soukannsyouryou," << TotalCellAveInter[i] << endl;
                                        }


                                        //総干渉量を比較し最大のセルを導出
                                        for (i = 0; i < SmallNum; i++)
                                        {

                                            //
                                            //cout << "\n";
                                            //cout << "check point" << i << "\n";
                                            //cout << "\n";
                                            //cout << "ADDFlag" << AddFlag[i] << "\n";


                                            //総干渉量を比較しセル半径を下げるセルの選出
                                            //電力制御していないかつ除外していないセルにおいて総干渉量を調べる
                                            if (SmallCell[i].PowerControlDoneFlag == 0)
                                            {
                                                if (ComparisonTotalCellAveInter == 0.0)
                                                {
                                                    ComparisonTotalCellAveInter = TotalCellAveInter[i];
                                                    MaxTotalCellAveInterNum = i;
                                                }
                                                else if (ComparisonTotalCellAveInter < TotalCellAveInter[i])
                                                {
                                                    ComparisonTotalCellAveInter = TotalCellAveInter[i];
                                                    MaxTotalCellAveInterNum = i;
                                                }
                                                //cout << "MaxTotalCellAveInterNum" << MaxTotalCellAveInterNum << "\n";
                                                //cout << "ComparisonTotalCellAveInter" << ComparisonTotalCellAveInter << "\n";

                                                writing_file1 << "MaxTotalCellAveInterNum," << MaxTotalCellAveInterNum << endl;
                                                writing_file1 << "ComparisonTotalCellAveInter," << ComparisonTotalCellAveInter << endl;
                                            }

                                        }
                                        //もし初めて電力制御するときの最大総干渉量を持つセル番号は除外するときに必要なため保持しておく
                                        if (FirstTotalCellAveInterFlag == 0)
                                        {
                                            FirstTotalCellAveInterFlag = 1;
                                            FirstTotalCellAveInterNum = MaxTotalCellAveInterNum;

                                            //cout << "FirstTotalCellAveInterNum" << FirstTotalCellAveInterNum << "\n";
                                            writing_file1 << "FirstTotalCellAveInterNum," << FirstTotalCellAveInterNum << endl;
                                        }

                                        cout << "unchi" << endl;
                                        //最も総干渉量が多かったスモールセルの電力を下げる
                                        if (SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag == 0)
                                        {
                                            //総干渉量が全セルの中最大で電力制御したというflagを立てる

                                            SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag = 1;
                                            //cout << "PowerControlDoneFlag" << SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag << "\n";
                                            writing_file1 << "PowerControlDoneFlag," << SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag << endl;

                                            //対象のセルの送信電力をΔ分下げる
                                            //cout << "Tx_P[MaxTotalCellAveInterNum]" << Tx_P[MaxTotalCellAveInterNum] << "\n";
                                            writing_file1 << "Tx_P[MaxTotalCellAveInterNum]," << Tx_P[MaxTotalCellAveInterNum] << endl;

                                            //もし総干渉量がゼロのセルが選ばれていたらスルーする
                                            if (ComparisonTotalCellAveInter != 0.0)
                                            {
                                                beforeTx_P = Tx_P[MaxTotalCellAveInterNum];
                                                Tx_P[MaxTotalCellAveInterNum] -= DeltaTx_P;
                                                //cout << "seigyogoTx_P[MaxTotalCellAveInterNum]" << Tx_P[MaxTotalCellAveInterNum] << "\n";
                                                writing_file1 << "seigyogoTx_P[MaxTotalCellAveInterNum]," << Tx_P[MaxTotalCellAveInterNum] << endl;
                                            }
                                            else
                                            {
                                                writing_file1 << "soukannsyouryougazeronoserugaerabareteitanodede_dennryokuseigyohaokonawanai" << endl;
                                            }


                                            if (ComparisonTotalCellAveInter != 0.0)
                                            {
                                                //対象の干渉セルの送信電力の再設計
                                                //対象の干渉セルの半径の再設計
                                                beforeCellRadius = SmallCellRadius[MaxTotalCellAveInterNum];

                                                writing_file1 << "beforecellradius," << SmallCellRadius[MaxTotalCellAveInterNum] << endl;

                                                //beforeTx_P = Tx_P[MaxTotalCellAveInterNum];
                                                SmallCellRadius[MaxTotalCellAveInterNum] = RefDistance * pow(10.0, (Tx_P[MaxTotalCellAveInterNum] - Noise_dBm - SNR - 20 * log10(4 * PI * RefDistance / wavelength))/ (RefDistance * PathlossExponent));
                                                //Tx_P[MaxTotalCellAveInterNum] = SNR + Pathloss(SmallCellRadius[MaxTotalCellAveInterNum]) + Noise_dBm;

                                                //cout << "beforTx_P,Tx_P" << beforeTx_P << Tx_P[MaxTotalCellAveInterNum] << "\n";
                                                //cout << "beforeCellRadius, SmallCellRadius" << beforeCellRadius << SmallCellRadius[MaxTotalCellAveInterNum] << "\n";

                                                //fprintf(fpA, "beforTx_P,Tx_P[MaxTotalCellAveInterNum], %20.20f, %20.20f\n", beforeTx_P, Tx_P[MaxTotalCellAveInterNum]);
                                                writing_file1 << endl << "Execute PowerControl" << endl;
                                                writing_file1 << "BeforeCellRadius," << beforeCellRadius << ",SmallCellRadius[MaxTotalCellAveInterNum]," << SmallCellRadius[MaxTotalCellAveInterNum] << endl;
                                                writing_file1 << "InitialTx_P," << InitialTx_P << ",Tx_P[MaxTotalCellAveInterNum]," << Tx_P[MaxTotalCellAveInterNum] << ",10 * log10(CellAveOwn[MaxTotalCellAveInterNum])," << 10 * log10(CellAveOwn[MaxTotalCellAveInterNum]) << ",CellAveOwn[MaxTotalCellAveInterNum]," << CellAveOwn[MaxTotalCellAveInterNum] << endl;

                                                //初期送信電力と設計後の送信電力の差だけ平均電力を下げる
                                                for (i = 0; i < SmallNum; i++)
                                                {
                                                    if (CellAveInter[i][MaxTotalCellAveInterNum] != 0.0)
                                                    {
                                                        CellAveInter[i][MaxTotalCellAveInterNum] = pow(10.0, (10 * log10(CellAveInter[i][MaxTotalCellAveInterNum]) - (beforeTx_P - Tx_P[MaxTotalCellAveInterNum])) / 10.0);

                                                        //cout << "CellAveInter[" << i << "][" << MaxTotalCellAveInterNum << "]" << CellAveInter[i][MaxTotalCellAveInterNum] << "\n";
                                                    }
                                                }
                                            }

                                            NumofPowerControlDoneFlag++;

                                            //cout << "NumofPowerControlDoneFlag" << NumofPowerControlDoneFlag << "\n";
                                            writing_file1 << "NumofPowerControlDoneFlag," << NumofPowerControlDoneFlag << endl;

                                            //電力制御したというflagのリセット

                                        }

                                        //縮小したセルの最悪メッシュを距離と平均SINRより変更しそこの平均自局信号電力と各平均干渉電力を求める
                                        if ((ComparisonTotalCellAveInter != 0.0) && (SmallCell[MaxTotalCellAveInterNum].PowerControlDoneFlag == 1))
                                        {

                                            //初期化
                                            tmpMinGridAveSINR = 0.0;

                                            for (j = 0; j < GridNum; j++)
                                            {
                                                for (k = 0; k < GridNum; k++)
                                                {
                                                    //初期化
                                                    Grid[j][k].flag = 0;


                                                    //四隅分の距離を導出
                                                    Grid[j][k].dist[0] = TwoPdistance(SmallCell[MaxTotalCellAveInterNum].x, SmallCell[MaxTotalCellAveInterNum].y, j  * (MapSize / GridNum), k * (MapSize / GridNum));
                                                    Grid[j][k].dist[1] = TwoPdistance(SmallCell[MaxTotalCellAveInterNum].x, SmallCell[MaxTotalCellAveInterNum].y, (j + 1) * (MapSize / GridNum), k * (MapSize / GridNum));
                                                    Grid[j][k].dist[2] = TwoPdistance(SmallCell[MaxTotalCellAveInterNum].x, SmallCell[MaxTotalCellAveInterNum].y, j * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));
                                                    Grid[j][k].dist[3] = TwoPdistance(SmallCell[MaxTotalCellAveInterNum].x, SmallCell[MaxTotalCellAveInterNum].y, (j + 1) * (MapSize / GridNum), (k + 1) * (MapSize / GridNum));


                                                    //四隅の点と中心点間距離，スモールセルの半径を比較，半径以下だったら旗を立てる
                                                    if ((SmallCellRadius[MaxTotalCellAveInterNum] >= Grid[j][k].dist[0]) || (SmallCellRadius[i] >= Grid[j][k].dist[1]) || (SmallCellRadius[i] >= Grid[j][k].dist[2]) || (SmallCellRadius[i] >= Grid[j][k].dist[3]))
                                                    {
                                                        Grid[j][k].flag = 1;
                                                    }


                                                    //旗が1となる（最悪グリッドとなりうる可能性がある）グリッドでセンサをばらまき試行を重ね，SINRの平均が最も低いグリッドを探す

                                                    //条件式（旗が1ならば実行）
                                                    if (Grid[j][k].flag == 1)
                                                    {

                                                        //平均SINRを比較し最悪メッシュの導出
                                                        if (tmpMinGridAveSINR == 0.0)
                                                        {
                                                            tmpMinGridAveSINR = Grid[j][k].tmpAveSINR[MaxTotalCellAveInterNum];
                                                            tmpMinGridx[MaxTotalCellAveInterNum] = j;
                                                            tmpMinGridy[MaxTotalCellAveInterNum] = k;
                                                        }
                                                        else if (tmpMinGridAveSINR > Grid[j][k].tmpAveSINR[MaxTotalCellAveInterNum])
                                                        {
                                                            tmpMinGridAveSINR = Grid[j][k].tmpAveSINR[MaxTotalCellAveInterNum];
                                                            tmpMinGridx[MaxTotalCellAveInterNum] = j;
                                                            tmpMinGridy[MaxTotalCellAveInterNum] = k;
                                                        }
                                                    }
                                                }
                                            }

                                            //cout << "worstmesh_change" << tmpMinGridx[MaxTotalCellAveInterNum] << "," << tmpMinGridy[MaxTotalCellAveInterNum] << "\n";
                                            writing_file1 << "worstmesh_change," <<tmpMinGridx[MaxTotalCellAveInterNum] << "," << tmpMinGridy[MaxTotalCellAveInterNum] << endl;

                                            //最悪メッシュの平均自局受信電力と平均他局信号電力を格納しなおす

                                            //平均受信電力にセルを縮小する前とした後の電力の差分を引く
                                            CellAveOwn[MaxTotalCellAveInterNum] = pow(10.0, (10 * log10(Grid[tmpMinGridx[MaxTotalCellAveInterNum]][tmpMinGridy[MaxTotalCellAveInterNum]].AveOwn[MaxTotalCellAveInterNum]) - (InitialTx_P - Tx_P[MaxTotalCellAveInterNum])) / 10.0);
                                            //cout << "CellAveOwn[MaxTotalCellAveInterNum]" << CellAveOwn[MaxTotalCellAveInterNum] << "\n";

                                            writing_file1 << "CellAveOwn[MaxTotalCellAveInterNum]," << CellAveOwn[MaxTotalCellAveInterNum] << endl;

                                            for (j = 0; j < SmallNum; j++)
                                            {
                                                if ((MaxTotalCellAveInterNum != j) && (Grid[tmpMinGridx[MaxTotalCellAveInterNum]][tmpMinGridy[MaxTotalCellAveInterNum]].AveInter[MaxTotalCellAveInterNum][j] != 0.0))
                                                {
                                                    CellAveInter[MaxTotalCellAveInterNum][j] = pow(10.0, (10 * log10(Grid[tmpMinGridx[MaxTotalCellAveInterNum]][tmpMinGridy[MaxTotalCellAveInterNum]].AveInter[MaxTotalCellAveInterNum][j]) - (InitialTx_P - Tx_P[j])) / 10.0);

                                                    //cout << "CellAveInter[MaxTotalCellAveInterNum][j] " << CellAveInter[MaxTotalCellAveInterNum][j] << "\n";
                                                    writing_file1 << "CellAveInter[MaxTotalCellAveInterNum][j]," << CellAveInter[MaxTotalCellAveInterNum][j] << endl;
                                                }

                                            }
                                        }

                                        for (i = 0; i < SmallNum; i++)
                                        {
                                            SpectrumSharingImpossible_pc[i] = 0;
                                        }
                                        writing_file1 << "OtherPower_mW," << OtherPower_mW << ",AveIsum," << AveIsum << endl;

                                        for (i = 0; i < SmallNum; i++)
                                        {
                                            writing_file1 << endl << "SmallNum," << i << endl;
                                            writing_file1 << "CellAveOwn[i]," << CellAveOwn[i] << endl;


                                            //試行回数分の瞬時受信電力を求める（メッシュ内の受信電力平均値にフェージングによる瞬時変動を加えたもの）
                                            for (j = 0; j < loop * SensorNum; j++)
                                            {
                                                OwnPower = CellAveOwn[i] * RayleighFading();
                                                OwnPower_nofade = CellAveOwn[i];

                                                //干渉電力の総和を求める
                                                for (l = 0; l < SmallNum; l++)
                                                {
                                                    OtherPower_mW += CellAveInter[i][l] * RayleighFading();
                                                    CellAveInter_nofade += CellAveInter[i][l];
                                                }

                                                AveIsum += OtherPower_mW;
                                                AveSum_nofade += CellAveInter_nofade;

                                                //(干渉電力の和)と雑音電力の和
                                                InterferenceNoise = OtherPower_mW + Noise_mW;
                                                InterferenceNoise_nofade = CellAveInter_nofade + Noise_mW;

                                                AveSINR += OwnPower / InterferenceNoise;
                                                AveSINR_nofade += OwnPower_nofade / InterferenceNoise_nofade;

                                                SINR = 10 * log10(OwnPower / InterferenceNoise);

                                                //	cout << "SINR" << SINR << "\n";

                                                instantaneousSINR[m] = SINR;
                                                m++;

                                                if (SINR < DesiredSINR)
                                                {
                                                    SpectrumSharingImpossible_pc[i]++;
                                                }

                                                //諸々初期化
                                                OwnPower = 0.0;
                                                OwnPower_nofade = 0.0;
                                                OtherPower_mW = 0.0;
                                                CellAveInter_nofade = 0.0;
                                                OtherPower_dBm = 0.0;
                                                SINR = 0.0;
                                            }

                                            m = 0;

                                            //cout << "check point4" << "\n";

                                            sort(instantaneousSINR, instantaneousSINR + SIZE_OF_ARRAY(instantaneousSINR));

                                            //cout << "instantaneousSINR" << instantaneousSINR[1299] << "\n";

                                            AveIsum /= (loop * SensorNum);
                                            AveSum_nofade /= (loop * SensorNum);
                                            AveSINR /= (loop * SensorNum);
                                            AveSINR_nofade /= (loop * SensorNum);

                                            writing_file1 << "AveIsum," << AveIsum << ",AveSINR," << AveSINR << "," << 10 * log10(AveSINR) << endl;
                                            writing_file1 << "AveSum_nofade," << AveSum_nofade << ",AveSINR_nofade," << AveSINR_nofade << endl;
                                            //cout << "AveIsum,AveSINR" << AveIsum << "," << AveSINR << "\n";

                                            //writing_file1 << "平均SINR," << 10 * log10(AveSINR) << endl;
                                            //cout << CDFofTimes[1300] << "\n";



                                            //SIRthと所望SIRの差分の導出
                                            DeltaP[i] = instantaneousSINR[OUTAGE] - DesiredSINR;
                                            InitialDeltaP[i] = instantaneousSINR[OUTAGE] - DesiredSINR;

                                            writing_file1 << "SINRth," << instantaneousSINR[OUTAGE] << ",DeltaP," << DeltaP[i] << endl;
                                            //cout << "SINRth" << instantaneousSINR[4999] << "DeltaP" << DeltaP[i] << "\n";
                                            //cout << "AveIsum" << AveIsum << "\n";


                                            //共用可否判定とFlag立て-------------------------------------------------------------------------------------------------------------------------------------------------
                                            if (DeltaP[i] >= 0.0)
                                            {
                                                PossibleSharingCounter++;

                                                ////電力制御なしで共用できた数
                                                //if (NonPowerControlCounter == 0)
                                                //{
                                                //	PossibleSharingCounter_NonPC++;
                                                //}
                                                //cout << "SharingAchievement!!!" << "\n";

                                            }
                                            AveIsum = 0.0;
                                            AveSum_nofade = 0.0;
                                            AveSINR = 0.0;
                                            AveSINR_nofade = 0.0;
                                        }


                                        writing_file1 << endl;
                                        //通信路容量の導出
                                        //cmin = log(1 + (pow(10.0,-75.0 / 10.0)) / pow(10, Noise_dBm / 10));
                                        //c = log(1 + (pow(10.0, Tx_P[] / 10.0)) / pow(10, Noise_dBm / 10));


                                        //すべてのセルで電力制御を行った場合でもスループットが最小値より大きく，さらに一から電力制御を繰り返す場合

                                        //スループットの計算
                                        //cf = (BAND_MAX * pow(10.0, 9.0) + BAND_MIN * pow(10.0, 9.0)) / 2.0;
                                        writing_file1 << "wavelength," << wavelength << ",cf," << cf << endl;

                                        //セル端のスループットの計算を行う
                                        for (k= 0; k < SmallNum; k++)
                                        {
                                            rx_edge[k] = Tx_P[k] - Pathloss(InitialSmallCellRadius, wavelength);
                                            writing_file1 << "Tx_P[k]," << Tx_P[k] << ",SmallCellNum," << k << ",rx_edge," << rx_edge[k] << ",InitialSmallCellRadius," << InitialSmallCellRadius << ",SmallCellUseBand," << SmallCell[k].UseBand << endl;
                                            writing_file1 << "rx_edge[k]," << rx_edge[k] << endl;

                                            snr_edge[k] = pow(10.0, rx_edge[k] / 10.0) / pow(10.0, NOISE_DB / 10.0);

                                            writing_file1 << "SmallCellNum," << k << ",snr_edge," << snr_edge[k] << endl;

                                            SmallCell[k].Throughput = (INITIAL_BANDWIDTH) * log2(1 + snr_edge[k]);

                                            if (pc_count != 1)
                                            {
                                                SmallCell[k].PowerControlThroughput = SmallCell[k].Throughput;
                                            }

                                            writing_file1 << "SmallCellNum," << k << ",Throughput," << SmallCell[k].Throughput << ",Throughput[Mbps]," << SmallCell[k].Throughput / pow(10.0, 6.0) << endl;
                                            c_sum += SmallCell[k].Throughput;
                                        }

                                        writing_file1 << "ThroughputSum," << c_sum << "," << c_sum / pow(10.0, 6.0) << endl;
                                        writing_file1 << "AverageThroughput," << c_sum / pow(10.0, 6.0) / SmallNum << endl;

                                        //すべて共用できていたらbreak(抜く必要あり)
                                        if (PossibleSharingCounter == SmallNum)
                                        {
                                            SharingAchievementFlag = 1;
                                            break;
                                        }

                                        for (k = 0; k < SmallNum; k++)
                                        {
                                            cout << SmallCell[k].Throughput / pow(10.0, 6.0) << endl;
                                            if ((SmallCell[k].Throughput / pow(10.0, 6.0)) < 10)
                                            {
                                                cout << "lower than 10Mpbs" << endl;
                                                writing_file1 << "SmallCell[k].Throughput_is_lower_than_10Mbps," << SmallCell[k].Throughput / pow(10.0, 6.0) << endl;
                                                breakflag = 1;
                                                break;
                                            }
                                        }

                                        if (breakflag == 1)
                                        {
                                            break;
                                        }

                                        if (NumofPowerControlDoneFlag == SmallNum)
                                        {
                                            for (i = 0; i < SmallNum; i++)
                                            {
                                                SmallCell[i].PowerControlDoneFlag = 0;
                                            }
                                            NumofPowerControlDoneFlag = 0;
                                        }

                                        //諸々初期化
                                        ComparisonTotalCellAveInter = 0.0;
                                        MaxTotalCellAveInterNum = 0;
                                        PossibleSharingCounter = 0;
                                        c_sum = 0.0;

                                        for (i = 0; i < SmallNum; i++)
                                        {
                                            TotalCellAveInter[i] = 0.0;
                                            rx_edge[i] = 0.0;
                                            snr_edge[i] = 0.0;
                                            SmallCell[i].Throughput = 0.0;

                                        }
                                    }

                                    //諸々初期化
                                    NumofPowerControlDoneFlag = 0;

                                    FirstTotalCellAveInterFlag = 0;


                                    //諸々リセット
                                    for (i = 0; i < SmallNum; i++)
                                    {
                                        for (j = 0; j < SmallNum; j++)
                                        {
                                            CellAveInter[i][j] = InitialCellAveInter[i][j];
                                        }
                                    }
                                */

                                }

                                cout << "whiel no soto" << endl;

                                //諸々初期化
                                PossibleSharingCounter_NonPC = 0;


                                //cout << "SharingAchievement" << "\n";
                                writing_file1 << "\n**************************************************************************************************************************************************************************" << endl;
                                writing_file1 << "SharingAchievement(PowerControlOnly)!!!" << endl;
                                writing_file1 << "**************************************************************************************************************************************************************************\n" << endl;
                                //cout << "\n";

                                for (i = 0; i < SmallNum; i++)
                                {
                                    pout_nonpc[i] = (double)SpectrumSharingImpossible_nonpc[i] / (double)(loop * SensorNum);
                                    pout[i] = (double)SpectrumSharingImpossible_pc[i] / (double)(loop * SensorNum);
                                    pout_nonpc_initial[i] = (double)SpectrumSharingImpossible_nonpc_initial[i] / (double)(loop * SensorNum);

                                    writing_file1 << "pout_nonpc_initial[i]," << pout_nonpc_initial[i] << ",pout_nonpc[i]," << pout_nonpc[i] << ",pout[i]," << pout[i] << ",DoNotSharingCellFlag[i]," << DoNotSharingCellFlag[i] << endl;
                                }

                                for (i = 0; i < SmallNum; i++)
                                {
                                    //アウテージ確率の初期化
                                    pout_nonpc[i] = 0.0;
                                    pout[i] = 0.0;
                                    pout_nonpc_initial[i] = 0.0;

                                    writing_file1 << "pout_nonpc_initial[i]," << pout_nonpc_initial[i] << ",pout_nonpc[i]," << pout_nonpc[i] << ",pout[i]," << pout[i] << ",DoNotSharingCellFlag[i]," << DoNotSharingCellFlag[i] << endl;
                                }


                        }
                        else
                        {
                            //このバンドを使用している基地局は1つ
                            writing_file1 << "numf," << numf[y] << endl;
                            //トータルスループットの初期化
                            c_sum = 0.0;
                            GB_c_sum = 0.0;
                            //スループットの計算を行う
                            for (k = 0; k < SmallNum; k++)
                            {
                                writing_file1 << "SmallCell[k].UseBand," << SmallCell[k].UseBand << ",f[y]," << f[y] << endl;
                                if (SmallCell[k].UseBand != (f[y] * pow(10.0, 9.0)))
                                {
                                    wavelength = 3.0 * pow(10.0, 8.0) / (((BAND_MAX + BAND_MIN) / 2.0) * pow(10.0, 9.0));
                                    writing_file1 << "wavelength," << wavelength << ",SmallCell[k].UseBand," << (((BAND_MAX + BAND_MIN) / 2.0) * pow(10.0, 9.0)) << endl;
                                    Tx_P[k] = SNR + Pathloss(InitialSmallCellRadius, wavelength) + Noise_dBm;
                                    rx_edge[k] = Tx_P[k] - Pathloss(InitialSmallCellRadius, wavelength);
                                    writing_file1 << "Tx_P[k]," << Tx_P[k] << ",SmallCellNum," << k << ",rx_edge," << rx_edge[k] << ",InitialSmallCellRadius," << InitialSmallCellRadius << ",SmallCellUseBand," << SmallCell[k].UseBand << endl;
                                    writing_file1 << "SmallCellNum," << k << ",rx_edge," << rx_edge[k] << ",InitialSmallCellRadius," << InitialSmallCellRadius << ",SmallCellUseBand," << SmallCell[k].UseBand << endl;
                                    //cout << "rx_edge[k]" << rx_edge[k] << endl;

                                    snr_edge[k] = pow(10.0, rx_edge[k] / 10.0) / pow(10.0, NOISE_DB / 10.0);

                                    writing_file1 << "SmallCellNum," << k << ",snr_edge," << snr_edge[k] << endl;
                                    //cout << "snr_edge[k]" << snr_edge[k] << endl;

                                    if (SmallCell[k].AllSpectrumBandFlag == 0)
                                    {
                                        writing_file1 << "SmallCell[k].AllSpectrumBandFlag=," << SmallCell[k].AllSpectrumBandFlag << endl;
                                        SmallCell[k].Throughput = (INITIAL_BANDWIDTH) * log2(1 + snr_edge[k]);
                                        SmallCell[k].Throughput_GB = (INITIAL_BANDWIDTH) * log2(1 + snr_edge[k]);
                                    }
                                    else
                                    {
                                        SmallCell[k].Throughput = (INITIAL_BANDWIDTH / DividedNum) * log2(1 + snr_edge[k]);
                                        SmallCell[k].Throughput_GB = ((INITIAL_BANDWIDTH - (GB * pow(10.0, 6.0) * (DividedNum - 1))) / DividedNum) * log2(1 + snr_edge[k]);
                                    }


                                    writing_file1 << "SmallCellNum," << k << ",Throughput," << SmallCell[k].Throughput / pow(10.0, 6.0) << "Mbps" << endl;
                                    writing_file1 << "SmallCellNum," << k << ",GB_Throughput," << SmallCell[k].Throughput_GB / pow(10.0, 6.0) << "Mbps" << endl;

                                    c_sum += SmallCell[k].Throughput / pow(10.0, 6.0);
                                    GB_c_sum += SmallCell[k].Throughput_GB / pow(10.0, 6.0);
                                }
                            }
                            writing_file1 << "ThroughputSum," << c_sum << ",AverageThroughput," << c_sum / SmallNum << endl;
                            writing_file1 << "GB_ThroughputSum," << GB_c_sum << ",AverageThroughput," << GB_c_sum / SmallNum << endl;
                            GB_c_sum = 0.0;
                            c_sum = 0.0;
                        }

                    }
                }


            }

            for (j = 0; j < SmallNum; j++)
            {
                if (SmallCell[j].UseBand != (f[y] * pow(10.0, 9.0)))
                {
                    if (SmallCell[j].AllSpectrumBandFlag == 0)
                    {
                        SmallCell[j].AlreadyAllocatedBandFlag = 0;
                    }
                }
            }

            //初期化
            c_sum = 0.0;
            numf[0] = numf[1] = numcf = 0;
        }
    }
    else
    {
        //トータルスループットの初期化
        c_sum = 0.0;
        GB_c_sum = 0.0;
        //スループットの計算を行う
        for (k = 0; k < SmallNum; k++)
        {
            wavelength = 3.0 * pow(10.0, 8.0) / (((BAND_MAX + BAND_MIN) / 2.0) * pow(10.0, 9.0));
            writing_file1 << "wavelength," << wavelength << ",SmallCell[k].UseBand," << (((BAND_MAX + BAND_MIN) / 2.0) * pow(10.0, 9.0)) << endl;
            Tx_P[k] = SNR + Pathloss(InitialSmallCellRadius, wavelength) + Noise_dBm;
            rx_edge[k] = Tx_P[k] - Pathloss(InitialSmallCellRadius, wavelength);
            writing_file1 << "Tx_P[k]," << Tx_P[k] << ",SmallCellNum," << k << ",rx_edge," << rx_edge[k] << ",InitialSmallCellRadius," << InitialSmallCellRadius << ",SmallCellUseBand," << SmallCell[k].UseBand << endl;
            writing_file1 << "SmallCellNum," << k << ",rx_edge," << rx_edge[k] << ",InitialSmallCellRadius," << InitialSmallCellRadius << ",SmallCellUseBand," << SmallCell[k].UseBand << endl;
            //cout << "rx_edge[k]" << rx_edge[k] << endl;

            snr_edge[k] = pow(10.0, rx_edge[k] / 10.0) / pow(10.0, NOISE_DB / 10.0);

            writing_file1 << "SmallCellNum," << k << ",snr_edge," << snr_edge[k] << endl;
            //cout << "snr_edge[k]" << snr_edge[k] << endl;

            if (SmallCell[k].AllSpectrumBandFlag == 0)
            {
                writing_file1 << "SmallCell[k].AllSpectrumBandFlag=," << SmallCell[k].AllSpectrumBandFlag << endl;
                SmallCell[k].Throughput = (INITIAL_BANDWIDTH) * log2(1 + snr_edge[k]);
                SmallCell[k].Throughput_GB = (INITIAL_BANDWIDTH) * log2(1 + snr_edge[k]);
            }
            else
            {
                SmallCell[k].Throughput = (INITIAL_BANDWIDTH / DividedNum) * log2(1 + snr_edge[k]);
                SmallCell[k].Throughput_GB = ((INITIAL_BANDWIDTH - (GB * pow(10,6) * (DividedNum-1))) / DividedNum) * log2(1 + snr_edge[k]);
            }


            writing_file1 << "SmallCellNum," << k << ",Throughput," << SmallCell[k].Throughput / pow(10.0, 6.0) << "Mbps" << endl;
            writing_file1 << "SmallCellNum," << k << ",GB_Throughput," << SmallCell[k].Throughput_GB / pow(10.0, 6.0) << "Mbps" << endl;

            c_sum += SmallCell[k].Throughput / pow(10.0, 6.0);
            GB_c_sum += SmallCell[k].Throughput_GB /pow(10.0, 6.0);
        }


        writing_file1 << "ThroughputSum," << c_sum << ",AverageThroughput," << c_sum / SmallNum << endl;
        writing_file1 << "GB_ThroughputSum," << GB_c_sum << ",AverageThroughput," << GB_c_sum / SmallNum << endl;
        c_sum = 0.0;
        GB_c_sum = 0.0;
    }
    return 0;
}
