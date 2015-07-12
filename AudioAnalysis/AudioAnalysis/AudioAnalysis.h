// AudioAnalysis.h
//define the main functions of modula "AudioAnalysis" in Python
//@author: Liuyt

#ifdef EXPORT_AUDIOANALYSIS_DLL
#define AUDIOANALYSIS_API __declspec(dllexport)
#else 
#define AUDIOANALYSIS_API __declspec(dllimport)
#endif

extern "C"
{
	AUDIOANALYSIS_API int getMusicAttributes(short int*, int, float*);
	AUDIOANALYSIS_API int beatFilter(short int* , int , bool, float* , bool );
	AUDIOANALYSIS_API int differenceDetection(short int*, int, float, float*, bool);
	AUDIOANALYSIS_API int makeSilence(short int*, int, float, float, float*, bool);
	AUDIOANALYSIS_API int lowPassFilter(short int*, int, float, float, float*, bool);
	AUDIOANALYSIS_API int highPassFilter(short int*, int, float, float*, bool);
	AUDIOANALYSIS_API int getBeats(short int*, int, float, float, float*, bool);
	AUDIOANALYSIS_API int convertToIntensity(short int*, int, float*, bool);
	AUDIOANALYSIS_API int countBeats(short int*, int, int, bool, float*, bool);
	AUDIOANALYSIS_API int combineMP3Frames(char*, int, char*);
	AUDIOANALYSIS_API int getMelIntensity(float*, int, int*, int, float* , int, float);
	AUDIOANALYSIS_API int getMFCCs(float*, int, int, float*);
	AUDIOANALYSIS_API int getAmplitude(double*, double*, int, float*);
	AUDIOANALYSIS_API int getComplexRealFromStr(char *, int);
}