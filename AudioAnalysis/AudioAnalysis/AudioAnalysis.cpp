// AudioAnalysis.cpp : 定义 DLL 应用程序的导出函数。
//

//return:
//0, completed
//1, memory ERROR
//2, NULL input
//3, open file ERROR

#include "stdafx.h"
#define EXPORT_AUDIOANALYSIS_DLL
#include "AudioAnalysis.h"
#include "stdio.h"
#include <cmath>

#define SIMPLE_RATE 44100

// the info field
/*
AVE_BIG_WINDOW_LENGTH = musicAttributes[0]
AVERAGE_INTENSITY_L = musicAttributes[1]
AVERAGE_INTENSITY_R = musicAttributes[2]
AVERAGE_INTENSITY = musicAttributes[3]
MAX_VALUE_L = musicAttributes[4]
MAX_VALUE_R = musicAttributes[5]
MAX_VALUE = musicAttributes[6]
MIN_VALUE_L = musicAttributes[7]
MIN_VALUE_R = musicAttributes[8]
MIN_VALUE = musicAttributes[9]
MAX_L / AVERAGE_L = musicAttributes[10]
MAX_R / AVERAGE_R = musicAttributes[11]
MAX / AVERAGE = musicAttributes[12]
NONZERO_EXCEPTSILENCE_LENGTH = musicAttributes[13]// NOT USED
AVE_BEAT_INTENSITY_L = musicAttributes[14]
AVE_BEAT_INTENSITY_H = musicAttributes[15]
BEATS_L = musicAttributes[16]
BEATS_H = musicAttributes[17]
LOW_PASS_AVE_REDUCE_RATIO = musicAttributes[18]
HIGH_PASS_AVE_REDUCE_RATIO = musicAttributes[19]
*/

//FILE* info;
//errno_t err_info;

/*  -- Function: getMusicAttributes
	-- short int* data: the incoming music short data
	-- int len: length of the incoming data
	-- float* musicAttributes: the info field

	this function calculates the basic attributes of music
*/
AUDIOANALYSIS_API int getMusicAttributes(short int* data, int len, float* musicAttributes)
{
//	err_info = fopen_s(&info, "H:\\info.txt", "w");
	float dataSum_L = 0.0;
	float dataSum_R = 0.0;
	int nonzeroDataCount_L = 0;
	int nonzeroDataCount_R = 0;

	float absMax_L = 0.0;
	float absMax_R = 0.0;
	float absMax = 0.0;

//	musicAttributes[0] = 0.0;// keep this value//AVE_BIG_WINDOW_LENGTH
	musicAttributes[1] = 0.0;
	musicAttributes[2] = 0.0;
	musicAttributes[3] = 0.0;
	musicAttributes[4] = -1000000;
	musicAttributes[5] = -1000000;
	musicAttributes[6] = -1000000;
	musicAttributes[7] = 1000000;
	musicAttributes[8] = 1000000;
	musicAttributes[9] = 1000000;
	musicAttributes[10] = 0.0;
	musicAttributes[11] = 0.0;
	musicAttributes[12] = 0.0;
//	musicAttributes[13] = 0.0;// keep this value// NONZERO_EXCEPTSILENCE_LENGTH = musicAttributes[13]// NOT USED
//	musicAttributes[14] = 0.0;// keep this value// AVE_BEAT_INTENSITY_L = musicAttributes[14]
//	musicAttributes[15] = 0.0;// keep this value// AVE_BEAT_INTENSITY_H = musicAttributes[15]
//  musicAttributes[16] = 0.0;// keep this value// BEATS_L = musicAttributes[16]
//	musicAttributes[17] = 0.0;// keep this value// BEATS_H = musicAttributes[17]
//	musicAttributes[18] = 0.0;// keep this value// LOW_PASS_AVE_REDUCE_RATIO = musicAttributes[18]
//	musicAttributes[19] = 0.0;// keep this value// HIGH_PASS_AVE_REDUCE_RATIO = musicAttributes[19]

	if (data == nullptr || musicAttributes == nullptr || len == 0)
	{
		return 2;// NULL input ERROR
	}

	for (int i = 0; i < len; i = i + 2)
	{
		//if (i >= 0)
		//{
		//	fprintf(info, "%d\n", data[i]);
		//}
		if (data[i] != 0)
		{
			dataSum_L = dataSum_L + abs(data[i]);
			nonzeroDataCount_L ++;
		}
		if (data[i + 1] != 0)
		{
			dataSum_R = dataSum_R + abs(data[i + 1]);
			nonzeroDataCount_R++;
		}

		if (data[i] > musicAttributes[4])// MAX_VALUE_L
		{
			musicAttributes[4] = data[i];// MAX_VALUE_L
		}

		if (data[i] < musicAttributes[7])// MIN_VALUE_L
		{
			musicAttributes[7] = data[i];// MIN_VALUE_L
		}

		if (data[i + 1] > musicAttributes[5])// MAX_VALUE_R
		{
			musicAttributes[5] = data[i + 1];// MAX_VALUE_R
		}

		if (data[i + 1] < musicAttributes[8])// MIN_VALUE_R
		{
			musicAttributes[8] = data[i + 1];// MIN_VALUE_R
		}
	}

	musicAttributes[1] = dataSum_L / nonzeroDataCount_L;// AVERAGE_INTENSITY_L
	musicAttributes[2] = dataSum_R / nonzeroDataCount_R;// AVERAGE_INTENSITY_R
	musicAttributes[3] = (musicAttributes[1] + musicAttributes[2]) / 2.0;// AVERAGE_INTENSITY
	musicAttributes[6] = (musicAttributes[4] > musicAttributes[5]) ? (musicAttributes[4]) : (musicAttributes[5]);// MAX_VALUE_L
	musicAttributes[9] = (musicAttributes[7] > musicAttributes[8]) ? (musicAttributes[7]) : (musicAttributes[8]);// MAX_VALUE_R
	absMax_L = (musicAttributes[4] > abs(musicAttributes[7])) ? (musicAttributes[4]) : abs((musicAttributes[7]));// ABS_MAX_VALUE_L
	absMax_R = (musicAttributes[5] > abs(musicAttributes[8])) ? (musicAttributes[5]) : abs((musicAttributes[8]));// ABS_MAX_VALUE_R
	absMax = (musicAttributes[6] > abs(musicAttributes[9])) ? (musicAttributes[6]) : abs((musicAttributes[9]));// ABS_MAX_VALUE
	musicAttributes[10] = absMax_L / musicAttributes[1];// MAX_L / AVERAGE_L
	musicAttributes[11] = absMax_R / musicAttributes[2];// MAX_R / AVERAGE_R
	musicAttributes[12] = absMax / musicAttributes[3];// MAX / AVERAGE
	//musicAttributes[13] = (float)(nonzeroDataCount_L + nonzeroDataCount_R) / 
	//					  (float)(nonZeroExceptSilenceLen_L + nonZeroExceptSilenceLen_R);// NONZERO_EXCEPTSILENCE_LENGTH

	//if (info)
	//{
	//	fclose(info);
	//}

	return 0;
}

// debug output txt file

//FILE* fvarVar_L;
//errno_t err_varVar_L;
//FILE* fvarVar_R;
//errno_t err_varVar_R;
//
//FILE* fHvarVar_L;
//errno_t err_HvarVar_L;
//FILE* fHvarVar_R;
//errno_t err_HvarVar_R;
//
//FILE* fvar_L;
//errno_t err_var_L;
//FILE* fvar_R;
//errno_t err_var_R;

/*  -- Function: deviation
	-- short int* localExp: the incoming array
	-- int len: length of the incoming array
	-- bool channel: true means the data is from left channel
				     false means the data is from right channel
	-- float* musicAttributes: the info field

	this function calculates the standard deviation(aqrt(variance)) of the given array
*/
float deviation(short int* localExp, int len, bool channel, float* musicAttributes)
{
	float sum = 0;
	short int expectation = 0;
	float weight = 0;

	for (int i = 0; i < len; i++)
	{
		sum = sum + localExp[i];
	}
	expectation = (short int)(sum / (float)len);
	sum = 0;
	for (int i = 0; i < len; i++)
	{
		sum = sum + (localExp[i] - expectation) * (localExp[i] - expectation);
	}
	if (expectation != 0)
	{
		if (channel)// left channel
		{
			weight = musicAttributes[1] / (float)expectation;
		}
		else// right channel
		{
			weight = musicAttributes[2] / (float)expectation;
		}
	}
	if (weight < 1)
	{
		weight = 1;
	}
	
	if (channel)// left channel
	{
		return sqrt(((sum / (float)len) * weight));
	}
	else// right channel
	{
		return sqrt(((sum / (float)len) * weight));
	}
}

/*  -- Function: deviation_deviation
	-- float* localVar: the incoming array
	-- int len: length of the incoming array
	-- bool channel: true means the data is from left channel
					 false means the data is from right channel
	-- float* musicAttributes: the info field

	this function calculates the standard deviation(aqrt(variance)) of the given array
	more or less the same function as float deciation(), the only defference is the parameter localExp/localVar
*/
float deviation_deviation(float* localVar, int len, bool channel, float* musicAttributes)
{
	float sum = 0;
	float expectation = 0;

	for (int i = 0; i < len; i++)
	{
		sum = sum + localVar[i];
	}
	expectation = sum / (float)len;
	sum = 0;
	for (int i = 0; i < len; i++)
	{
		sum = sum + (localVar[i] - expectation) * (localVar[i] - expectation);
	}

	if (channel)// left channel
	{
		return sqrt(sum / (float)len);
	}
	else// right channel
	{
		return sqrt(sum / (float)len);
	}
}

/*  -- Function: homogenization
	-- float* incoming: the incoming array
	-- int len: length of the incoming array
	-- int radius: the radius of homogenization, the bigger this value is, the smoother the array is
	-- bool style: true means the data is from left channel
				   false means the data is from right channel
				   (only used in debug, when the output txt file is needed)

	this function makes the incoming array smoother by using the parameter radius
*/
int homogenization(float* incoming, int len, int radius, bool style)
{
	float* temp = nullptr;

	if (!(temp = new float[len]))
	{
		return 1;// memory allocation ERROR
	}
		
	for (int i = 0; i < len; i++)
	{
		temp[i] = incoming[i];
	}

	for (int i = 0; i < len; i++)
	{
		if (i >= radius && i <= (len - radius) - 1)
		{
			float leftSide = 0.0;
			float rightSide = 0.0;
			for (int j = 1; j <= radius; j++)
			{
				leftSide = leftSide + temp[i - j];
			}
			for (int j = 1; j <= radius; j++)
			{
				rightSide = rightSide + temp[i + j];
			}
			incoming[i] = (leftSide + rightSide + incoming[i]) / (2 * (float)radius + 1);
		}

		//if (style)// left channel
		//{
		//	fprintf(fHvarVar_L, "%f\n", incoming[i]);
		//}
		//if (!style)// right channel
		//{
		//	fprintf(fHvarVar_R, "%f\n", incoming[i]);
		//}
	}
	delete[] temp;
}

/*  -- Function: beatFilter
	-- short int* data: the incoming music short data
	-- int len: length of the incoming data
	-- bool style: true means the data is from left channel
				   false means the data is from right channel
	-- float* musicAttributes: the info field
	-- bool updateMusicAttributes: true means update the info field, and vice versa

		this function recieves the music data stream that is in format of short int[], and uses self-adaptive filter algorithm to
	recognize beats in music. At the same time, the second deviation of music stream can be used to detect if the data stream is
	stable or not. If the data stream is stable means the stream fragment doesn't have much beats, so the rate of convergence of
	self-adaptive filter should be decreased, otherwise the error rate will increase.
*/
AUDIOANALYSIS_API int beatFilter(short int* data, int len, bool style, float* musicAttributes, bool updateMusicAttributes)
{
	// parameters for self adaptive filter
	float adaptive[2] = { 1.0, 1.0 };
	
	float threshValue[2] = { adaptive[0] * musicAttributes[1], adaptive[1] * musicAttributes[2] };
	int maxOfShort = 0;

	int smallWindow_start[2] = { 0, 0 };
	int smallWindow_end[2] = { 0, 0 };
	int bigWindow_start[2] = { 0, 0 };
	int bigWindow_end[2] = { 0, 0 };
	int bigWindowCount[2] = { 0, 0};
	int bigWindowLength[2] = { 0, 0 };
	float AVE_bigWindowLength[2] = { 0.0, 0.0};
	int nonZeroCount[2] = { 0, 0 };
	float nonZeroRation[2] = { 0.0, 0.0 };

	int smallWindowLength = 0;
	float step = 0.0;
	float nonZeroThreshVal = 0.0;

	bool isBeat[2] = { true, true };
	int zeroCount[2] = { 0 , 0};
	int bigWindow = 0;
	float decreaseFactor[2] = { 0, 0 };

	int multiple[2] = { 0, 0 };
	int multipleCheck[2] = { 0, 0 };

	// parameters for prolonged sound detection
	float checkWindowTime = 0.001;
	float varFrameLength = 0.02;
	float varVarFrameLength = 0.1;
	float checkWindowLength[2] = { checkWindowTime * SIMPLE_RATE, checkWindowTime * SIMPLE_RATE };
	int localExpLength = (int)(varFrameLength / checkWindowTime);
	int localVarLength = (int)(varVarFrameLength / varFrameLength);
	int localVarCount[2] = { 0, 0 };
	float* localVar_L = nullptr;
	float* localVar_R = nullptr;
	short int* localExp_L = nullptr;
	short int* localExp_R = nullptr;
	int checkWindowCount[2] = { 0, 0 };
	int checkSum[2] = { 0, 0 };
	int checkCount[2] = { 0, 0 };
	int varVectorLength = (int)((len) / (varFrameLength * SIMPLE_RATE)) + 1;
	int varVarVectorLength = (int)((varVectorLength) / (localVarLength)) + 1;
	int varCount[2] = { 0, 0 };
	int varVarCount[2] = { 0, 0 };
	float* var_L = nullptr;
	float* var_R = nullptr;
	float* varVar_L = nullptr;
	float* varVar_R = nullptr;
	float gradient[2] = { 0.0, 0.0 };
	float varVarThresh_high = 40;
	float varVarThresh_low = 500;
	float parMinThreshVal = 0;
	float minThreshVal = 0.0;

	if (!(localExp_L = new short int[localExpLength]))
		return 1;// memory allocation ERROR
	if (!(localExp_R = new short int[localExpLength]))
		return 1;// memory allocation ERROR
	if (!(localVar_L = new float[localVarLength]))
		return 1;// memory allocation ERROR
	if (!(localVar_R = new float[localVarLength]))
		return 1;// memory allocation ERROR
	if (!(var_L = new float[varVectorLength]))
		return 1;// memory allocation ERROR
	if (!(var_R = new float[varVectorLength]))
		return 1;// memory allocation ERROR
	if (!(varVar_L = new float[varVarVectorLength]))
		return 1;// memory allocation ERROR
	if (!(varVar_R = new float[varVarVectorLength]))
		return 1;// memory allocation ERROR
	
	if (!style)// False: highPass stream
	{
		smallWindowLength = 0.01 * SIMPLE_RATE;
		step = 0.25;
		nonZeroThreshVal = 0.05;
		AVE_bigWindowLength[0] = 2 * SIMPLE_RATE;
		AVE_bigWindowLength[1] = 2 * SIMPLE_RATE;
		parMinThreshVal = 200;
	}
	if (style)// True: lowPass stream
	{
		smallWindowLength = 0.05 * SIMPLE_RATE;
		step = 0.1;
		nonZeroThreshVal = 0.05;
		AVE_bigWindowLength[0] = 1 * SIMPLE_RATE;
		AVE_bigWindowLength[1] = 1 * SIMPLE_RATE;
		parMinThreshVal = 1400;
	}

	if (data == nullptr || musicAttributes == nullptr || len == 0)
	{
		return 2;// NULL input ERROR
	}
	
	//err_var_L = fopen_s(&fvar_L, "H:\\var_L.txt", "w");
	//err_var_R = fopen_s(&fvar_R, "H:\\var_R.txt", "w");
	//err_varVar_L = fopen_s(&fvarVar_L, "H:\\varVar_L.txt", "w");
	//err_varVar_R = fopen_s(&fvarVar_R, "H:\\varVar_R.txt", "w");
	//err_HvarVar_L = fopen_s(&fHvarVar_L, "H:\\HvarVar_L.txt", "w");
	//err_HvarVar_R = fopen_s(&fHvarVar_R, "H:\\HvarVar_R.txt", "w");

	for (int i = 0; i < len - 1; i = i + 2)
	{
		//  ---- calculate the variance in every 0.5s of left channel ----  //
		checkSum[0] = checkSum[0] + abs(data[i]);
		checkCount[0]++;
		if (checkCount[0] >= (int)checkWindowLength[0])
		{
			localExp_L[checkWindowCount[0]] = (short int)((float)checkSum[0] / (float)checkCount[0]);
			checkWindowCount[0]++;
			if (checkWindowCount[0] >= localExpLength)
			{
				var_L[varCount[0]] = deviation(localExp_L, localExpLength, true, musicAttributes);// left channel
				localVar_L[localVarCount[0]] = var_L[varCount[0]];
				localVarCount[0]++;
				if (localVarCount[0] >= localVarLength)
				{
					varVar_L[varVarCount[0]] = deviation_deviation(localVar_L, localVarLength, true, musicAttributes);// left channel
					//fprintf(fvarVar_L, "%f\t", varVar_L[varVarCount[0]]);
					varVarCount[0]++;
					localVarCount[0] = 0;
				}
				checkWindowCount[0] = 0;
//				fprintf(fvar_L, "%f\n", var_L[varCount[0]]);
				varCount[0]++;
			}
			checkSum[0] = 0;
			checkCount[0] = 0;
		}

		//  ---- calculate the variance in every 0.5s of right channel ----  //
		checkSum[1] = checkSum[1] + abs(data[i + 1]);
		checkCount[1]++;
		if (checkCount[1] >= (int)checkWindowLength[1])
		{
			localExp_R[checkWindowCount[1]] = (short int)((float)checkSum[1] / (float)checkCount[1]);
			checkWindowCount[1]++;
			if (checkWindowCount[1] >= localExpLength)
			{
				var_R[varCount[1]] = deviation(localExp_R, localExpLength, false, musicAttributes);// right channel
				localVar_R[localVarCount[1]] = var_R[varCount[1]];
				localVarCount[1]++;
				if (localVarCount[1] >= localVarLength)
				{
					varVar_R[varVarCount[1]] = deviation_deviation(localVar_R, localVarLength, false, musicAttributes);// right channel
//					fprintf(fvarVar_R, "%f\n", varVar_R[varVarCount[1]]);
					varVarCount[1]++;
					localVarCount[1] = 0;
				}
				checkWindowCount[1] = 0;
//				fprintf(fvar_R, "%f\n", var_R[varCount[1]]);
				varCount[1]++;
			}
			checkSum[1] = 0;
			checkCount[1] = 0;
		}
	}

	if (homogenization(varVar_L, varVarVectorLength / 2, 2, true) == 1)
	{
		return 1;// memory allocation ERROR
	}
	if (homogenization(varVar_R, varVarVectorLength / 2, 2, false) == 1)
	{
		return 1;// memory allocation ERROR
	}

	for (int i = 0; i < len - 1; i = i + 2)
	{
		// ----  filter  ---- // 
		if (abs(data[i]) < threshValue[0])
		{
			data[i] = 0;
		}
		if (abs(data[i + 1]) < threshValue[1])
		{
			data[i + 1] = 0;
		}
		//  ---- filter finish ----  //

		//  ---- LEFT CHANNEL ----  //
		varVarCount[0] = (int)((int)(i / (int)(varFrameLength * SIMPLE_RATE * 2)) / localVarLength);
		if (style)// low pass
		{
			gradient[0] = varVar_L[varVarCount[0]] / varVarThresh_low;
		}
		else// high pass
		{
			gradient[0] = varVar_L[varVarCount[0]] / varVarThresh_high;
		}
		//  ---- calculate AVE big window length ----  //
		if (data[i] == 0)
		{
			zeroCount[0]++;
		}
		else
		{
			if (isBeat[0])
			{
				if (bigWindow_end[0] != 0)// not origin
				{
					bigWindow_start[0] = i;
					bigWindowLength[0] = bigWindowLength[0] + (bigWindow_start[0] - bigWindow_end[0]);
					bigWindowCount[0]++;
					AVE_bigWindowLength[0] = (float)bigWindowLength[0] / (float)bigWindowCount[0];

					float temp = (AVE_bigWindowLength[0] + AVE_bigWindowLength[1]) / 2;
					AVE_bigWindowLength[0] = temp;
					AVE_bigWindowLength[1] = temp;
				}
				isBeat[0] = false;
			}
			bigWindow_end[0] = i;
			zeroCount[0] = 0;
		}
		if (zeroCount[0] >= 0.08 * SIMPLE_RATE)
		{
			if (style)// low pass
			{
				//if (varVarCount[0] >= 5)
				if (varVar_L[varVarCount[0]] > varVarThresh_low / 2)
				{
					isBeat[0] = true;
				}
			}
			else// high pass
			{
				if (varVar_L[varVarCount[0]] > varVarThresh_high / 2)
				{
					isBeat[0] = true;
				}
			}
		}
		if (!style && (zeroCount[0] >= AVE_bigWindowLength[0] / 3))// high pass
		{
			multiple[0] = (int)((float)zeroCount[0] / (AVE_bigWindowLength[0] / 1000));
			if (multipleCheck[0] != multiple[0] && multiple[0] >= 400)
			{
				if (varVar_L[varVarCount[0]] > varVarThresh_high / 2)
				{
					if (gradient[0] > 3)
					{
						gradient[0] = 3;
					}
					decreaseFactor[0] = (zeroCount[0] / (AVE_bigWindowLength[0] / 3)) * 0.3;
					adaptive[0] = adaptive[0] - 0.01 * step * pow(2, decreaseFactor[0]) * gradient[0];
				}
				multipleCheck[0] = multiple[0];
			}
		}
		if (style && (zeroCount[0] >= AVE_bigWindowLength[0] / 10))// low pass
		{
			multiple[0] = (int)((float)zeroCount[0] / (AVE_bigWindowLength[0] / 1000));
			if (multipleCheck[0] != multiple[0])
			{
				if (varVar_L[varVarCount[0]] > varVarThresh_low / 2)
				{
					if (gradient[0] > 1.5)
					{
						gradient[0] = 1.5;
					}
					decreaseFactor[0] = (zeroCount[0] / (AVE_bigWindowLength[0] / 10)) * 0.3;
					adaptive[0] = adaptive[0] - 0.01 * step * pow(2, decreaseFactor[0]) * gradient[0];
				}
				multipleCheck[0] = multiple[0];
			}
		}
		//  ---- calculate nonzero ratio in small window ----  //
		smallWindow_start[0] = i;
		if (smallWindow_start[0] - smallWindow_end[0] == smallWindowLength * 2)
		{
			smallWindow_end[0] = smallWindow_start[0] - smallWindowLength * 2;
		}

		if (data[smallWindow_start[0]] != 0)
		{
			nonZeroCount[0]++;
		}
		if ((smallWindow_start[0] - smallWindow_end[0] == smallWindowLength * 2) && data[smallWindow_end[0]] != 0)
		{
			nonZeroCount[0]--;
		}
		nonZeroRation[0] = (float)nonZeroCount[0] / (float)smallWindowLength;

		if (nonZeroRation[0] > nonZeroThreshVal)
		{
			if (!style)// high pass
			{
				adaptive[0] = adaptive[0] + 5 * (nonZeroRation[0] / nonZeroThreshVal) * step;
			}
			else// low pass
			{
				adaptive[0] = adaptive[0] + 5 * (nonZeroRation[0] / nonZeroThreshVal) * step;
			}
			smallWindow_end[0] = i;
			nonZeroCount[0] = 0;
		}
		if (adaptive[0] > musicAttributes[10])
		{
			adaptive[0] = musicAttributes[10];
		}

		if (varVar_L[varVarCount[0]] != 0)
		{
			minThreshVal = parMinThreshVal / varVar_L[varVarCount[0]];
		}
		else
		{
			minThreshVal = 1000;
		}

		if (adaptive[0] < minThreshVal)
		{
			adaptive[0] = minThreshVal;
		}
		//  ---- left channel finish ----  //

		//  ---- RIGHT CHANNEL ----  //
		varVarCount[1] = (int)((int)((i + 1) / (int)(varFrameLength * SIMPLE_RATE * 2)) / localVarLength);
		if (style)// low pass
		{
			gradient[1] = varVar_R[varVarCount[1]] / varVarThresh_low;
		}
		else// high pass
		{
			gradient[1] = varVar_R[varVarCount[1]] / varVarThresh_high;
		}
		//  ---- calculate AVE big window length ----  //
		if (data[i + 1] == 0)
		{
			zeroCount[1]++;
		}
		else
		{
			if (isBeat[1])
			{
				if (bigWindow_end[1] != 0)// not origin
				{
					bigWindow_start[1] = i + 1;
					bigWindowLength[1] = bigWindowLength[1] + (bigWindow_start[1] - bigWindow_end[1]);
					bigWindowCount[1]++;
					AVE_bigWindowLength[1] = (float)bigWindowLength[1] / (float)bigWindowCount[1];

					float temp = (AVE_bigWindowLength[0] + AVE_bigWindowLength[1]) / 2;
					AVE_bigWindowLength[0] = temp;
					AVE_bigWindowLength[1] = temp;
				}
				isBeat[1] = false;
			}
			bigWindow_end[1] = i + 1;
			zeroCount[1] = 0;
		}
		if (zeroCount[1] >= 0.08 * SIMPLE_RATE)
		{
			if (style)// low pass
			{
				if (varVar_R[varVarCount[1]] > varVarThresh_low / 2)
				{
					isBeat[1] = true;
				}
			}
			else// high pass
			{
				if (varVar_R[varVarCount[1]] > varVarThresh_high / 2)
				{
					isBeat[1] = true;
				}
			}
		}
		if (!style && (zeroCount[1] >= AVE_bigWindowLength[1] / 3))// high pass
		{
			multiple[1] = (int)((float)zeroCount[1] / (AVE_bigWindowLength[1] / 1000));
			if (multipleCheck[1] != multiple[1] && multiple[1] >= 400)
			{
				if (varVar_R[varVarCount[1]] > varVarThresh_high / 2)
				{
					if (gradient[1] > 3)
					{
						gradient[1] = 3;
					}
					decreaseFactor[1] = (zeroCount[1] / (AVE_bigWindowLength[1] / 3)) * 0.3;
					adaptive[1] = adaptive[1] - 0.01 * step * pow(2, decreaseFactor[1]) * gradient[1];
				}
				multipleCheck[1] = multiple[1];
			}
		}
		if (style && (zeroCount[1] >= AVE_bigWindowLength[1] / 10))// low pass
		{
			multiple[1] = (int)((float)zeroCount[1] / (AVE_bigWindowLength[1] / 1000));
			if (multipleCheck[1] != multiple[1])
			{
				if (varVar_R[varVarCount[1]] > varVarThresh_low / 2)
				{
					if (gradient[1] > 1.5)
					{
						gradient[1] = 1.5;
					}
					decreaseFactor[1] = (zeroCount[1] / (AVE_bigWindowLength[1] / 10)) * 0.3;
					adaptive[1] = adaptive[1] - 0.01 * step * pow(2, decreaseFactor[1]) * gradient[1];
				}
				multipleCheck[1] = multiple[1];
			}
		}
		//  ---- calculate nonzero ratio in small window ----  //
		smallWindow_start[1] = i + 1;
		if (smallWindow_start[1] - smallWindow_end[1] == smallWindowLength * 2)
		{
			smallWindow_end[1] = smallWindow_start[1] - smallWindowLength * 2;
		}

		if (data[smallWindow_start[1]] != 0)
		{
			nonZeroCount[1]++;
		}
		if ((smallWindow_start[1] - smallWindow_end[1] == smallWindowLength * 2) && data[smallWindow_end[1]] != 0)
		{
			nonZeroCount[1]--;
		}
		nonZeroRation[1] = (float)nonZeroCount[1] / (float)smallWindowLength;

		if (nonZeroRation[1] > nonZeroThreshVal)
		{
			if (!style)// high pass
			{
				adaptive[1] = adaptive[1] + 5 * (nonZeroRation[1] / nonZeroThreshVal) * step;
			}
			else// low pass
			{
				adaptive[1] = adaptive[1] + 5 * (nonZeroRation[1] / nonZeroThreshVal) * step;
			}
			smallWindow_end[1] = i + 1;
			nonZeroCount[1] = 0;
		}
		if (adaptive[1] > musicAttributes[11])
		{
			adaptive[1] = musicAttributes[11];
		}

		if (varVar_R[varVarCount[1]] != 0)
		{
			minThreshVal = parMinThreshVal / varVar_R[varVarCount[1]];
		}
		else
		{
			minThreshVal = 1000;
		}
		
		if (adaptive[1] < minThreshVal)
		{
			adaptive[1] = minThreshVal;
		}

		//  ---- right channel finish ----  //

		threshValue[0] = adaptive[0] * musicAttributes[1];
		threshValue[1] = adaptive[1] * musicAttributes[2];
	}
	musicAttributes[0] = (AVE_bigWindowLength[0] + AVE_bigWindowLength[1]) / (2 * SIMPLE_RATE);// AVE_BIG_WINDOW_LENGTH

	//if (fvar_L)
	//{
	//	fclose(fvar_L);
	//}
	//if (fvar_R)
	//{
	//	fclose(fvar_R);
	//}
	/*if (fvarVar_L)
	{
		fclose(fvarVar_L);
	}*/
	//if (fvarVar_R)
	//{
	//	fclose(fvarVar_R);
	//}
	//if (fHvarVar_L)
	//{
	//	fclose(fHvarVar_L);
	//}
	//if (fHvarVar_R)
	//{
	//	fclose(fHvarVar_R);
	//}

	if (updateMusicAttributes)
	{
		getMusicAttributes(data, len, musicAttributes);
	}

	//localExp_L = nullptr;
	//localExp_R = nullptr;
	//localVar_L = nullptr;
	//localVar_R = nullptr;
	//var_L = nullptr;
	//var_R = nullptr;
	//varVar_L = nullptr;
	//varVar_R = nullptr;

	delete[] localExp_L;
	delete[] localExp_R;
	delete[] localVar_L;
	delete[] localVar_R;
	delete[] var_L;
	delete[] var_R;
	delete[] varVar_L;
	delete[] varVar_R;

	return 0;
}

/*  -- Function: differenceDetection
	-- short int* data: the incoming music short data
	-- int len: length of the incoming data
	-- float ratio: the enlargement factor, the factor controls the range of output
	-- float* musicAttributes: the info field
	-- bool updateMusicAttributes: true means update the info field, and vice versa

		this function recieves the music data stream that is in format of short int[], and calculates the difference of each
	adjacent value.
*/
AUDIOANALYSIS_API int differenceDetection(short int* data, int len, float ratio, float* musicAttributes, bool updateMusicAttributes)
{
	float diff = 0.0;
	float tempL = 0.0, tempR = 0.0;

	int maxOfShort = 0;

	tempL = data[0];
	tempR = data[1];

	if (data == nullptr || musicAttributes == nullptr || len == 0)
	{
		return 2;// NULL input ERROR
	}

	for (int i = 2; i < len - 1; i = i + 2)
	{
		diff = abs(data[i] - tempL);
		tempL = data[i];
		maxOfShort = (int)(diff * ratio);
		if (maxOfShort >= 32767.0)
		{
			data[i] = 32767.0;
		}
		else
		{
			data[i] = (float)maxOfShort;
		}
		if (data[i] < 0.0)
		{
			data[i] = 0.0;
		}

		diff = abs(data[i + 1] - tempR);
		tempR = data[i + 1];
		maxOfShort = (int)(diff * ratio);
		if (maxOfShort >= 32767.0)
		{
			data[i + 1] = 32767.0;
		}
		else
		{
			data[i + 1] = (float)maxOfShort;
		}
		if (data[i + 1] < 0.0)
		{
			data[i + 1] = 0.0;
		}
	}

	if (updateMusicAttributes)
	{
		getMusicAttributes(data, len, musicAttributes);
	}

	return 0;
}

/*  -- Function: makeSilence
	-- short int* data: the incoming music short data
	-- int len: length of the incoming data
	-- float delta_L: the thresh value of left channel
	-- float delta_R: the thresh value of right channel
	-- float* musicAttributes: the info field
	-- bool updateMusicAttributes: true means update the info field, and vice versa

		this function recieves the music data stream that is in format of short int[], and for each channel, if the value is 
	smaller the thresh value, then makes it zero.
*/
AUDIOANALYSIS_API int makeSilence(short int* data, int len, float delta_L, float delta_R, 
								  float* musicAttributes, bool updateMusicAttributes)
{
	if (data == nullptr || musicAttributes == nullptr || len == 0)
	{
		return 2;// NULL input ERROR
	}

	for (int i = 0; i < len - 1; i = i + 2)
	{
		if ((int(data[i]) < delta_L) && ((int(data[i]) > -delta_L)))
		{
			data[i] = 0;
		}
		if ((int(data[i + 1]) < delta_R) && ((int(data[i + 1]) > -delta_R)))
		{
			data[i + 1] = 0;
		}
	}

	if (updateMusicAttributes)
	{
		getMusicAttributes(data, len, musicAttributes);
	}

	return 0;
}

/*  -- Function: lowPassFilter
	-- short int* data: the incoming music short data
	-- int len: length of the incoming data
	-- float alpha: the low pass parameter(normally this value is very small)
	-- float sigma: the enlargement factor
	-- float* musicAttributes: the info field
	-- bool updateMusicAttributes: true means update the info field, and vice versa

		this function recieves the music data stream that is in format of short int[], and uses the formula: 
	Xn = a * (Xn) + (1 - a) * Xn-1, to get the low pass output of original music stream
*/
AUDIOANALYSIS_API int lowPassFilter(short int* data, int len, float alpha, float sigma, float* musicAttributes, bool updateMusicAttributes)
{
	int maxShort = 0;
	float oriAVE = 0;

	data[0] = alpha * data[0];
	data[1] = alpha * data[1];

	if (data == nullptr || musicAttributes == nullptr || len == 0)
	{
		return 2;// NULL input ERROR
	}

	for (int i = 2; i < len; i++)
	{
		maxShort = (int)(alpha * (float)data[i] + (1.0 - alpha) * (float)data[i - 2]);
		if (maxShort >= 32767.0)
		{
			data[i] = 32767.0;
		}
		if (maxShort <= -32767.0)
		{
			data[i] = -32767.0;
		}
		if (maxShort < 32767.0 && maxShort > -32767.0)
		{
			data[i] = (short int)maxShort;
		}
		data[i - 2] = sigma * data[i - 2];
	}

	oriAVE = musicAttributes[3];

	if (updateMusicAttributes)
	{
		getMusicAttributes(data, len, musicAttributes);
	}

	if (oriAVE != 0)
	{
		musicAttributes[18] = musicAttributes[3] / oriAVE;
	}
	else
	{
		musicAttributes[18] = 1;
	}

	return 0;
}

/*  -- Function: highPassFilter
	-- short int* data: the incoming music short data
	-- int len: the low pass parameter(normally this value is very small)
	-- float alpha: the thresh value of left channel
	-- float* musicAttributes: the info field
	-- bool updateMusicAttributes: true means update the info field, and vice versa

		this function recieves the music data stream that is in format of short int[], and uses the formula:
	Xn = a * (Xn) + (1 - a) * Xn-1, to get the low pass output of original music stream, then minus by the original stream
	and get the high pass output of original music stream
*/
AUDIOANALYSIS_API int highPassFilter(short int* data, int len, float alpha, float* musicAttributes, 
									 bool updateMusicAttributes)
{
	short int* temp = nullptr;
	int maxShort = 0;
	float oriAVE = 0;

	if (!(temp = new short int[len]))
		return 1;// memory allocation ERROR

	for (int i = 0; i < len; i++)
	{
		temp[i] = data[i];
	}

	temp[0] = alpha * temp[0];
	temp[1] = alpha * temp[1];

	for (int i = 2; i < len; i++)
	{
//		temp[i] = alpha * temp[i] + (1.0 - alpha) * temp[i - 2];
//		maxShort = data[i] - temp[i];
		maxShort = (int)(((float)temp[i] - alpha * (float)temp[i - 2]));
		if (maxShort >= 32767.0)
		{
			data[i] = 32767.0;
		}
		if (maxShort <= -32767.0)
		{
			data[i] = -32767.0;
		}
		if (maxShort < 32767.0 && maxShort > -32767.0)
		{
			data[i] = maxShort;
		}
	}

	data[0] = temp[0];
	data[1] = temp[1];

	oriAVE = musicAttributes[3];

	if (updateMusicAttributes)
	{
		getMusicAttributes(data, len, musicAttributes);
	}

	if (oriAVE != 0)
	{
		musicAttributes[19] = musicAttributes[3] / oriAVE;
	}
	else
	{
		musicAttributes[19] = 1;
	}

	delete[] temp;
	return 0;
}

/*  -- Function: getBeats
	-- short int* data: the incoming music short data
	-- int len: the low pass parameter(normally this value is very small)
	-- float alpha_h: the alpha parameter in high pass filter
	-- float alpha_l: the alpha parameter in low pass filter
	-- float* musicAttributes: the info field
	-- bool updateMusicAttributes: true means update the info field, and vice versa

		this function recieves the music data stream that is in format of short int[], and centralized processes the data 
	by using the other functions
*/
AUDIOANALYSIS_API int getBeats(short int* data, int len, float alpha_h, float alpha_l, float* musicAttributes,
	bool updateMusicAttributes)
{
	short int* temp = nullptr;

	int maxShort_L = 0;
	int maxShort_R = 0;

	if (data == nullptr || len == 0)
	{
		return 2;// NULL input ERROR
	}

	if (!(temp = new short int[len]))
		return 1;// memory allocation ERROR

	for (int i = 0; i < len; i++)
	{
		temp[i] = data[i];
	}

	lowPassFilter(temp, len, alpha_l, 2, musicAttributes, true);

	beatFilter(temp, len, true, musicAttributes, false);
	countBeats(temp, len, SIMPLE_RATE * 0.08, true, musicAttributes, false);// low pass

	highPassFilter(data, len, alpha_h, musicAttributes, true);

	beatFilter(data, len, false, musicAttributes, false);
	countBeats(data, len, SIMPLE_RATE * 0.08, false, musicAttributes, false);// high pass

	for (int i = 2; i < len - 1; i = i + 2)
	{
		maxShort_L = abs(data[i]) + abs(temp[i]);
		maxShort_R = abs(data[i + 1]) + abs(temp[i + 1]);

		if (maxShort_L > 32767.0)
		{
			data[i] = 32767.0;
		}
		else
		{
			data[i] = maxShort_L;
		}
		if (maxShort_R > 32767.0)
		{
			data[i + 1] = 32767.0;
		}
		else
		{
			data[i + 1] = maxShort_R;
		}
	}

	data[0] = temp[0];
	data[1] = temp[1];

	if (updateMusicAttributes)
	{
		getMusicAttributes(data, len, musicAttributes);
	}

	delete[] temp;

	return 0;
}

/*  -- Function: convertToIntensity
	-- short int* data: the incoming music short data
	-- int len: the low pass parameter(normally this value is very small)
	-- float* musicAttributes: the info field
	-- bool updateMusicAttributes: true means update the info field, and vice versa

		this function recieves the music data stream that is in format of short int[], and calculates the intensity of the
	data(X^2)
*/
AUDIOANALYSIS_API int convertToIntensity(short int* data, int len, float* musicAttributes, bool updateMusicAttributes)
{
	if (data == nullptr || musicAttributes == nullptr || len == 0)
	{
		return 2;// NULL input ERROR
	}

	for (int i = 0; i < len; i++)
	{
		data[i] = abs(data[i]);
	}

	if (updateMusicAttributes)
	{
		getMusicAttributes(data, len, musicAttributes);
	}

	return 0;
}

/*  -- Function: countBeats
	-- short int* data: the incoming music short data
	-- int len: the low pass parameter(normally this value is very small)
	-- int windowSize: the minimum length of the gap between beats
	-- bool style: true means the data is from low pass filter, and vice versa
	-- float* musicAttributes: the info field
	-- bool updateMusicAttributes: true means update the info field, and vice versa

		this function recieves the music data stream that is in format of short int[], and count the beats that is already
	calculated by other function. And also calculates the average intensity of the beats.
*/
AUDIOANALYSIS_API int countBeats(short int* data, int len, int windowSize, bool style, float* musicAttributes, bool updateMusicAttributes)
{
	bool keep_L = true;
	int zeroCount_L = 0;
	bool keep_R = true;
	int zeroCount_R = 0;
	int beatsNum = 0;

	float sumBeatIntensity = 0;

	for (int i = 0; i < len; i = i + 2)
	{
		if (data[i] < 0)
		{
			data[i] = abs(data[i]);
		}
		if (data[i + 1] < 0)
		{
			data[i + 1] = abs(data[i + 1]);
		}

		if (data[i] == 0)
		{
			zeroCount_L++;
		}
		if (data[i] != 0)
		{
			if (zeroCount_L < windowSize && !keep_L)
			{
				data[i] = 0;
				zeroCount_L = 0;
			}
			if (keep_L)
			{
				beatsNum++;
				keep_L = false;
				zeroCount_L = 0;
			}
		}
		if (zeroCount_L >= windowSize)
		{
			zeroCount_L = 0;
			keep_L = true;
		}

		if (data[i + 1] == 0)
		{
			zeroCount_R++;
		}
		if (data[i + 1] != 0)
		{
			if (zeroCount_R < windowSize && !keep_R)
			{
				data[i + 1] = 0;
				zeroCount_R = 0;
			}
			if (keep_R)
			{
				beatsNum++;
				keep_R = false;
				zeroCount_R = 0;
			}
		}
		if (zeroCount_R >= windowSize)
		{
			zeroCount_R = 0;
			keep_R = true;
		}
		sumBeatIntensity = sumBeatIntensity + data[i];
		sumBeatIntensity = sumBeatIntensity + data[i + 1];
	}


	if (updateMusicAttributes)
	{
		getMusicAttributes(data, len, musicAttributes);
	}

	if (style == true)// low pass
	{
		if (beatsNum != 0)
		{
			musicAttributes[14] = (float)sumBeatIntensity / (float)beatsNum;// AVE_BEAT_INTENSITY_L
			musicAttributes[16] = (float)beatsNum * 1.0e5 / (float)len;// BEATS_L
		}
		else
		{
			musicAttributes[14] = 0;// AVE_BEAT_INTENSITY_L
			musicAttributes[16] = 0;// BEATS_L
		}
	}
	if (style == false)// high pass
	{
		if (beatsNum != 0)
		{
			musicAttributes[15] = (float)sumBeatIntensity / (float)beatsNum;// AVE_BEAT_INTENSITY_H
			musicAttributes[17] = (float)beatsNum * 1.0e5 / (float)len;// BEATS_H
		}
		else
		{
			musicAttributes[15] = 0;// AVE_BEAT_INTENSITY_H
			musicAttributes[17] = 0;// BEATS_H
		}
	}


	return 0;
}

/*  -- Function: CharToWchar
	-- const char* c: the incoming char array

		this function recieves the char* array, and transforms the char* array to Wchar* array
*/
wchar_t* CharToWchar(const char* c)
{
	wchar_t *m_wchar;
	int len = MultiByteToWideChar(CP_ACP, 0, c, strlen(c), NULL, 0);
	m_wchar = new wchar_t[len + 1];
	MultiByteToWideChar(CP_ACP, 0, c, strlen(c), m_wchar, len);
	m_wchar[len] = '\0';
	return m_wchar;
}

AUDIOANALYSIS_API int combineMP3Frames(char* path, int frameNum, char* frameData)
{
	FILE *readMP3 = NULL;
	errno_t err;
	wchar_t* MP3TxtPath = nullptr;
	if (!(MP3TxtPath = new wchar_t[strlen(path) + 1]))
		return 1;// memory allocation ERROR
	char temp = 0;

	MP3TxtPath = CharToWchar(path);

	if ((err = _wfopen_s(&readMP3, MP3TxtPath, TEXT("rb"))) != 0)// open return 0
	{
		return 3;// open file ERROR
	}
	for (int i = 0; i < 4096 * frameNum; i++)
	{
		if (!feof(readMP3))
		{
			frameData[i] = fgetc(readMP3);
		}
		else
		{
			break;
		}
	}

	if (readMP3)
	{
		fclose(readMP3);
	}
	
	readMP3 = NULL;
	return 0;
}

/*FILE* MFCC_info;
errno_t err_MFCC_info*/;


AUDIOANALYSIS_API int getMelIntensity(float* data, int len, int* melScale, int melWin, float* melIntensity, int simpleRate,
	float overlap)
{
//	err_MFCC_info = fopen_s(&MFCC_info, "H:\\MFCC_info.txt", "a");
	if (data == nullptr || melScale == nullptr || len == 0 || melIntensity == nullptr)
	{
		return 2;// NULL input ERROR
	}

//	float overlap = 0.5;
	int* melScaleInData = nullptr;
//	float* melIntensity = nullptr;

//	int focus_start = -1;
//	int focus_end = -1;

	if (!(melScaleInData = new int[melWin]))
		return 1;// memory allocation ERROR

	for (int i = 0; i < melWin; i++)
	{
		melScaleInData[i] = (int)(((float)melScale[i] / ((float)simpleRate / 2.0)) * len);
		//if (i >= 1)
		//{
		//	if (melScale[i - 1] < 100 && melScale[i] >= 100)
		//	{
		//		focus_start = i;
		//	}
		//	if (melScale[i - 1] < 1100 && melScale[i] >= 1400)
		//	{
		//		focus_end = i;
		//	}
		//}
	}
//	fprintf(MFCC_info, "start: %d\n", focus_start);
//	fprintf(MFCC_info, "end: %d\n", focus_end);

	int winStart = 0;
	int winEnd = 0;
	int winMid = 0;
	float weight = 0;
	float focusRate = 1.5;

	//fprintf(MFCC_info, "dataLength: %d\n", len);
	//for (int i = 0; i < 20; i++)
	//{
	//	fprintf(MFCC_info, "%d\n", melScaleInData[i]);
	//}
	//fprintf(MFCC_info, "\n\n");

	for (int window = 0; window < melWin; window++)
	{
		melIntensity[window] = 0;// set 0 for melIntensity[]
		if (window == 0)
		{
			winStart = 0;
		}
		if (window == 1)
		{
			winStart = melScaleInData[0] - int(overlap * (float)(melScaleInData[0]));
		}
		if (window > 1)
		{
			winStart = melScaleInData[window - 1] - int(overlap * (melScaleInData[window - 1] - melScaleInData[window - 2]));
		}
		winEnd = melScaleInData[window];
		winMid = winStart + int((winEnd - winStart) / 2);

		for (int i = winStart; i < winEnd; i++)
		{
			if (i < winMid)
			{
				//if (window >= focus_start && window <= focus_end)
				//{
				//	weight = float(i - winStart) / float(winMid - winStart) * focusRate;
				//}
				//else
				//{
				//	weight = float(i - winStart) / float(winMid - winStart);
				//}
				weight = float(i - winStart) / float(winMid - winStart);
			}
			else
			{
				//if (window >= focus_start && window <= focus_end)
				//{
				//	weight = float(winEnd - i) / float(winEnd - winMid) * focusRate;
				//}
				//else
				//{
				//	weight = float(winEnd - i) / float(winEnd - winMid);
				//}
				weight = float(winEnd - i) / float(winEnd - winMid);
			}
			melIntensity[window] = melIntensity[window] + data[i] * weight;
//			fprintf(MFCC_info, "%f * %f\n", data[i], weight);
		}

		//fprintf(MFCC_info, "%d start: %d\n", window, winStart);
		//fprintf(MFCC_info, "%d mid: %d\n", window, winMid);
		//fprintf(MFCC_info, "%d end: %d\n\n", window, winEnd);
	}

//	for (int i = 0; i < melWin; i++)
//	{
//		if (melIntensity[i] != 0)
//		{
//			melIntensity[i] = log(melIntensity[i]);
//		}
//		
////		fprintf(MFCC_info, "%f\n", melIntensity[i]);
//	}
////	fprintf(MFCC_info, "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////\n");
//
//	//  ---- start Discrete Cosine Transform (DCT) ----  //
//
//	for (int i = 0; i < MFCC_dimension; i++)
//	{
//		MFCCs[i] = 0;
//		for (int j = 0; j < melWin; j++)
//		{
//			MFCCs[i] = MFCCs[i] + melIntensity[j] * cos(((float)i * 3.1415926 / (float)melWin) * ((float)j - 0.5));
//		}
//		MFCCs[i] = MFCCs[i] * sqrt(2 / (float)melWin);
//	}

	//if (MFCC_info)
	//{
	//	fclose(MFCC_info);
	//}

	delete[] melScaleInData;
//	delete[] melIntensity;
//	melIntensity[0] = 0;
//	melIntensity[94] = 0;
//	melIntensity[95] = 0;
//	melIntensity[96] = 0;
//	melIntensity[97] = 0;
//	melIntensity[98] = 0;
//	melIntensity[99] = 0;

	return 0;
}

AUDIOANALYSIS_API int getMFCCs(float* melIntensity, int MFCC_dimension, int melWin, float* MFCCs)
{
		for (int i = 0; i < melWin; i++)
		{
			if (melIntensity[i] != 0)
			{
				melIntensity[i] = log(melIntensity[i]);
			}
			
	//		fprintf(MFCC_info, "%f\n", melIntensity[i]);
		}
	//	fprintf(MFCC_info, "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////\n");
	
		//  ---- start Discrete Cosine Transform (DCT) ----  //
	
		for (int i = 0; i < MFCC_dimension; i++)
		{
			MFCCs[i] = 0;
			for (int j = 0; j < melWin; j++)
			{
				MFCCs[i] = MFCCs[i] + melIntensity[j] * cos(((float)i * 3.1415926 / (float)melWin) * ((float)j - 0.5));
			}
			MFCCs[i] = MFCCs[i] * sqrt(2 / (float)melWin);
		}

		return 0;
}

//FILE* Amplitude_info;
//errno_t err_Amplitude_info;

AUDIOANALYSIS_API int getAmplitude(double* real, double* imag, int len, float* amplitude)
{
//	err_Amplitude_info = fopen_s(&Amplitude_info, "H:\\Amplitude_info.txt", "w");
	if (real == nullptr || imag == nullptr || len == 0 || amplitude == nullptr)
	{
		return 2;// NULL input ERROR
	}

	for (int i = 0; i < len; i++)
	{
		/*fprintf(Amplitude_info, "real %f\n", real[i]);
		fprintf(Amplitude_info, "imag %f\n\n", imag[i]);*/
		real[i] = real[i] * real[i];
		imag[i] = imag[i] * imag[i];
		amplitude[i] = sqrt(real[i] + imag[i]);
	}

	/*if (Amplitude_info)
	{
		fclose(Amplitude_info);
	}*/

	return 0;
}

// data format: (1.234e-10 + 5.678e-12j)
int getComplexRealFromStr(char *str, int len)
{
	for (int i = 0; i < len; i++)
	{
		if ((*(str + i) == 43 || *(str + i) == 45) && *(str + i - 1) != 101 && i != 1) // no "e" before "-" or "+"
		{
			return i;
		}
	}
	return -1;
}

