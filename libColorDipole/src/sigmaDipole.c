
//Wrapper for C

#ifdef __cplusplus
extern "C" {
#endif
  double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
#ifdef __cplusplus
};
#endif

double dipole_amplitude(double xBj, double r, double b, int param)
{
	return dipole_amplitude_(&xBj, &r, &b, &param);
}