
#ifndef _INITIAL_CONSTRUCTION_HPP_INCLUDED_
#define _INITIAL_CONSTRUCTION_HPP_INCLUDED_

#pragma once

//#define  oppai  Hatada 0.72;
#define SmallNum 3//�X���[���Z���̐�
#define InitialSmallCellRadius 150.0// �X���[���Z���̔��a

#define DeltaRadius 50.0 //�d�͐��䎞�ɋ��߂锼�a

#define MapSize 1000 //�}�b�v�̍ő勗��
#define GridNum 100  //�O���b�h�̐�

#define SensorNum 10 //�Z���T�[�̐�

#define DeltaTx_P 0.1 //�d�͐��䎞�ɉ�����

//#define Tx_P 30.0			/* ���M�d��[dBm] */
#define Pout 0.01			/* �A�E�e�[�W�m�� */

#define Margin 0			/* �����B�G���A�̃}�[�W�� */

#define PathlossExponent 3.5	/* �p�X���X�����W�� */
#define RefDistance	 10		/* �p�X���X�Q�Ƌ���10[m] */

#define MinTale		 -70.0 /* �V���h�E�C���O��pdf�̉����ŏ��l */
#define MaxTale		 100.0	/* �V���h�E�C���O��pdf�̉����ő�l */
#define ShadowMean	 0.0	/* �V���h�E�C���O�̑ΐ����K���z�̕��� */
#define ShadowSD	 8.0	/* �V���h�E�C���O�̑ΐ����K���z�̕W���΍� (standard deviation : SD) */
#define DeltaD		 0.1    /* ������ */

#define DesiredSINR	 5.0	/* ���]SINR[dB] */

#define MinTale1	 -40.0  /* SIR��pdf�̉����ŏ��l */
#define MaxTale1	 60.0   /* SIR��pdf�̉����ő�l */

#define PI			 3.141592653589793
#define loop		 1000

#define SNR			 30.0		 /* ���]SNR */
#define Boltzmann	 1.38064852  /*�{���c�}���萔 10^(-23)�Ȃ� */
#define Kelvin		 293.15		 /* ���̂̉��x 20�� */
#define FB			 100.0		 /* ���g���ш敝 100MHz 10^(6)�Ȃ� */

#define NOISE_DB -95.0
#define INITIAL_BANDWIDTH 100000000.0
#define BAND_MAX 4.7
#define BAND_MIN 4.6
#define DIVIDED_NUM 2

#define GB 5 /* �K�[�h�o���h 5MHz(MHz��Main��ł���)*/

#define DivideNumMax 10
#define DivideNumMin 2
typedef struct
{
	double x;//map��x���W
	double y;//map��y���W
	double dist[SmallNum];//�X���[���Z���Ԃ̋���
	int InterferenceArrivalAreaFlag[SmallNum];
	int PowerControlDoneFlag;
	int AllSpectrumBandFlag;
	double UseBand;//�g�p���S���g��
	double AlreadyAllocatedBandFlag; //���łɊ��蓖�Ă��o���h�����蓖�Ă��ꍇ�̓t���O����
	double Throughput;
	double Throughput_GB;
	double PowerControlThroughput;
	double BandDivisionThroughput;
}SmallCell;

typedef struct
{
	double x;//�ň��_��x���W
	double y;//�ň��_��y���W
	double difference_x; //�ň��_�ƒ��ڃX���[���Z����x���W�̍���
	double difference_y; //�ň��_�ƒ��ڃX���[���Z����y���W�̍���

}InterferenceWorstPoint;

typedef struct
{
	double x; //�Z���T�[��x���W
	double y; //�Z���T�[��y���W

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