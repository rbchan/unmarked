#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>


void getDetVec2(int y, double *detVec, double* mp) {
	if(y == 0) {
		//detVec[detVec_ind] *= 1;
		detVec[1] *= 1/(1 + exp(mp[0]));
	} else {
		detVec[0] = 0;
		detVec[1] *= exp(mp[0])/(1 + exp(mp[0]));
	}
}

void getDetVec4(int y, double *detVec, double* mp) {
	double den2 = 1 + exp(mp[0]), den3 = (1 + exp(mp[1]) + exp(mp[2])),
	den4 = (1 + exp(mp[3]) + exp(mp[4]) + exp(mp[5]));

	switch (y) {
	case 0:
	{
		//detVec[detVec_ind] *= 1;
		detVec[1] *= exp(mp[0])/den2;
		detVec[2] *= exp(mp[1])/den3;
		detVec[3] *= exp(mp[3])/den4;
	}
	break;
	case 1:
	{
		detVec[0] = 0;
		detVec[1] *= 1/den2;
		detVec[2] *= exp(mp[2])/den3;
		detVec[3] *= exp(mp[4])/den4;
	}
	break;
	case 2:
	{
		detVec[0] = 0;
		detVec[1] = 0;
		detVec[2] *= 1/den3;
		detVec[3] *= exp(mp[5])/den4;
	}
	break;
	case 3:
	{
		detVec[0] = 0;
		detVec[1] = 0;
		detVec[2] = 0;
		detVec[3] *= 1/den4;
	}
	break;
	}

}

SEXP getSingleDetVec(SEXP y_, SEXP mp_, SEXP K_) {
	int y = INTEGER_VALUE(y_), K = INTEGER_VALUE(K_) + 1;
	SEXP detVec;
	PROTECT(detVec = NEW_NUMERIC(K));
	double *mp = NUMERIC_POINTER(mp_), *detVecPtr = NUMERIC_POINTER(detVec);
	for(int k = 0; k != K; ++k) {
		detVecPtr[k] = 1;
	}
	switch (K) {
	case 2:
		getDetVec2(y, detVecPtr, mp);
		break;
	case 4:
		getDetVec4(y, detVecPtr, mp);
		break;
	}
	UNPROTECT(1);
	return(detVec);
}

SEXP getDetVecs(SEXP y_arr, SEXP mp_arr, SEXP J_i, SEXP tin, SEXP K_) {
	int *mp_dims = INTEGER_POINTER(GET_DIM(mp_arr));
	int nDMP = mp_dims[0], J = mp_dims[1], nY = mp_dims[2], M = mp_dims[3];
	int K = INTEGER_VALUE(K_) + 1;
	int t = INTEGER_VALUE(tin) - 1;
	void (*getDetVecPtr) (int, double*, double*);
	switch (K) {
	case 2:
		getDetVecPtr = getDetVec2;
		break;
	case 4:
		getDetVecPtr = getDetVec4;
		break;
	}
	SEXP detVec;
	PROTECT(detVec = NEW_NUMERIC(K*M));
	double *mp = NUMERIC_POINTER(mp_arr), *detVecPtr = NUMERIC_POINTER(detVec);
	int *J_ip = INTEGER_POINTER(J_i), *y = INTEGER_POINTER(y_arr);
	int detVec_ind = 0, mp_ind, y_ind;
	for (int i = 0; i != M; ++i) {
		for(int k = 0; k != K; ++k) {
			detVecPtr[detVec_ind + k] = 1;  // initialize to 1.
		}
		for(int j = 0; j != J_ip[i]; ++j) {
			y_ind = i + t*M + j*M*nY;
			mp_ind = j*nDMP + t*nDMP*J + i*nDMP*J*nY;
			if(y[y_ind] != 99) {
				getDetVecPtr(y[y_ind], detVecPtr + detVec_ind, mp + mp_ind);
			}
		}
		detVec_ind += K;
	}
	UNPROTECT(1);
	return detVec;
}
