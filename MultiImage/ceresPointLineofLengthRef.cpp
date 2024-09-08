#include "mex.h"
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

#include "ceres/ceres.h"
#include "glog/logging.h"

using namespace std;

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

//  ltl 2017/12/5
struct projectionPointError{
	projectionPointError(double observed_x, double observed_y, double normaliser):
		observed_x(observed_x), observed_y(observed_y), normaliser(normaliser) {}

	template <typename T>
	bool operator()(const T* const H, const T* const point, T* residuals) const {
		// H = [h1,h4,h7,h2,h5,h8,h3,h6,h9], point=[px,py]
		T predicted_x = (point[0] * H[0] + point[1] * H[3] + H[6]) / (point[0] * H[2] + point[1] * H[5] + H[8]);
		T predicted_y = (point[0] * H[1] + point[1] * H[4] + H[7]) / (point[0] * H[2] + point[1] * H[5] + H[8]);

		// The error is the difference between the predicted and observed position.
		residuals[0] = (1 / normaliser) * (predicted_x - T(observed_x));
		residuals[1] = (1 / normaliser) * (predicted_y - T(observed_y));

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x, const double observed_y, const double normaliser) {
		return (new ceres::AutoDiffCostFunction<projectionPointError, 2, 9, 2>(
			new projectionPointError(observed_x, observed_y, normaliser)));
	}

	double observed_x;
	double observed_y;
	double normaliser;
};


struct referencePointError{
	referencePointError(double observed_x, double observed_y, double normaliser):
		observed_x(observed_x), observed_y(observed_y), normaliser(normaliser) {}

	template <typename T>
	bool operator()(const T* const point, T* residuals) const {

		// The error is the difference between the predicted and observed position.
		residuals[0] = (1 / normaliser) * (point[0] - T(observed_x));
		residuals[1] = (1 / normaliser) * (point[1] - T(observed_y));

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x, const double observed_y, const double normaliser) {
		return (new ceres::AutoDiffCostFunction<referencePointError, 2, 2>(
			new referencePointError(observed_x, observed_y, normaliser)));
	}

	double observed_x;
	double observed_y;
	double normaliser;
};


struct projectionLineError {
	projectionLineError(double observed_a, double observed_b, double observed_c, double normaliser):
        observed_a(observed_a), observed_b(observed_b), observed_c(observed_c), normaliser(normaliser) {}
		
	template <typename T>
	bool operator()(const T* const H, const T* const linepoint, T* residuals) const {
		// H = [h1,h4,h7,h2,h5,h8,h3,h6,h9], linepoint=[px1, py1, px2, py2]
		T predicted_x1 = (linepoint[0] * H[0] + linepoint[1] * H[3] + H[6]) / (linepoint[0] * H[2] + linepoint[1] * H[5] + H[8]);
		T predicted_y1 = (linepoint[0] * H[1] + linepoint[1] * H[4] + H[7]) / (linepoint[0] * H[2] + linepoint[1] * H[5] + H[8]);
		T predicted_x2 = (linepoint[2] * H[0] + linepoint[3] * H[3] + H[6]) / (linepoint[2] * H[2] + linepoint[3] * H[5] + H[8]);
		T predicted_y2 = (linepoint[2] * H[1] + linepoint[3] * H[4] + H[7]) / (linepoint[2] * H[2] + linepoint[3] * H[5] + H[8]);

		// The error is the difference between the predicted and observed position.
		residuals[0] = (1 / normaliser) * (predicted_x1*T(observed_a) + predicted_y1*T(observed_b) + T(observed_c));
		residuals[1] = (1 / normaliser) * (predicted_x2*T(observed_a) + predicted_y2*T(observed_b) + T(observed_c));
        residuals[2] = (1 / normaliser) * (T(20.0) - sqrt((linepoint[0]-linepoint[2])*(linepoint[0]-linepoint[2])+(linepoint[1]-linepoint[3])*(linepoint[1]-linepoint[3])));
        
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_a, const double observed_b, const double observed_c, const double normaliser) {
		return (new ceres::AutoDiffCostFunction<projectionLineError, 3, 9, 4>(
			new projectionLineError(observed_a, observed_b, observed_c, normaliser)));
	}

	double observed_a;
	double observed_b;
	double observed_c;
	double normaliser;
};


struct referenceLineError {
	referenceLineError(double observed_a, double observed_b, double observed_c, double normaliser):
        observed_a(observed_a), observed_b(observed_b), observed_c(observed_c), normaliser(normaliser) {}
		
	template <typename T>
	bool operator()(const T* const linepoint, T* residuals) const {

		// The error is the difference between the predicted and observed position.
		residuals[0] = (1 / normaliser) * (linepoint[0]*T(observed_a) + linepoint[1]*T(observed_b) + T(observed_c));
		residuals[1] = (1 / normaliser) * (linepoint[2]*T(observed_a) + linepoint[3]*T(observed_b) + T(observed_c));
        residuals[2] = (1 / normaliser) * (T(20.0) - sqrt((linepoint[0]-linepoint[2])*(linepoint[0]-linepoint[2])+(linepoint[1]-linepoint[3])*(linepoint[1]-linepoint[3])));
        
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_a, const double observed_b, const double observed_c, const double normaliser) {
		return (new ceres::AutoDiffCostFunction<referenceLineError, 3, 4>(
			new referenceLineError(observed_a, observed_b, observed_c, normaliser)));
	}

	double observed_a;
	double observed_b;
	double observed_c;
	double normaliser;
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nrhs!=9) {
		mexErrMsgTxt("The number of input should 8");
	}
	if (nlhs != 1) {
		mexErrMsgTxt("The number of output should 1");
	}

	/* Input/output variables */
	double *origba_params;
	double *arr_xik_idx;
	int num_imgs;
	int ref_img;
	double *point_observations;
	double *point_normaliser;
	double *arr_lines_xik_idx;
	double *line_observations;
	double *line_normaliser;
	
	double *finalba_params;

	/* Intermediate variables.*/
	double *ba_params;
	mxArray *xik_idx[1];
	mxArray *lines_xik_idx[1];
	int ba_paramsm, observationsm_db2;
	int id_jc, id_xik, id_i, id_lines_xik;
	int num_Hparams, num_points, lenarr_xik_idx, lenarr_lines_xik_idx;
	int i, j, c = 0;
	int k = 0;
    int line_id_jc, line_c=0;

	/* Assign pointers to inputs. */
	origba_params = mxGetPr(prhs[0]);
	point_normaliser = mxGetPr(prhs[2]);
	num_points = mxGetM(prhs[2]);
	line_normaliser = mxGetPr(prhs[7]);
	num_imgs = mxGetScalar(prhs[3]);
	ref_img = mxGetScalar(prhs[4]);
	point_observations = mxGetPr(prhs[5]);
	line_observations = mxGetPr(prhs[8]);
	num_Hparams = 9 * (num_imgs - 1);

	/* Get sizes of input matrices (images, transformations, etc.).*/
	ba_paramsm = mxGetM(prhs[0]);
	observationsm_db2 = mxGetM(prhs[5]) / 2;

	ba_params = (double*)malloc(ba_paramsm * sizeof(double));
	memcpy(ba_params, origba_params, ba_paramsm * sizeof(double));

	plhs[0] = mxCreateDoubleMatrix(ba_paramsm, 1, mxREAL);
	finalba_params = mxGetPr(plhs[0]);

	/* Start computations. */
	Problem problem;
	for (i = 0; i<num_imgs; i++)
	{
		xik_idx[0] = mxGetCell(prhs[1], i);
		lines_xik_idx[0] = mxGetCell(prhs[6], i);
		arr_xik_idx = mxGetPr(xik_idx[0]);
		arr_lines_xik_idx = mxGetPr(lines_xik_idx[0]);
		lenarr_xik_idx = mxGetN(xik_idx[0]);
		lenarr_lines_xik_idx = mxGetN(lines_xik_idx[0]);
		if(i==ref_img-1)
		{
			for (j = 0; j<lenarr_xik_idx; j++)
			{
				id_xik = arr_xik_idx[j] - 1;
				id_jc = j + c;

				ceres::CostFunction* ref_point_cost_function =
					referencePointError::Create(point_observations[id_jc],
						point_observations[id_jc + observationsm_db2],
						point_normaliser[id_xik]);

				problem.AddResidualBlock(ref_point_cost_function,
					NULL /* squared loss */,
					ba_params + num_Hparams + (id_xik * 2));
			}
			
			for (j = 0; j < lenarr_lines_xik_idx; j++)
			{
				id_lines_xik = arr_lines_xik_idx[j] - 1;
				line_id_jc = j+line_c;
				ceres::CostFunction* ref_line_cost_function =
					referenceLineError::Create(line_observations[line_id_jc],
						line_observations[line_id_jc + lenarr_lines_xik_idx],
						line_observations[line_id_jc + 2*lenarr_lines_xik_idx],
						line_normaliser[id_lines_xik]);

				problem.AddResidualBlock(ref_line_cost_function, NULL, ba_params + num_Hparams + 2*num_points + (4*id_lines_xik));
			}			
		}
		else
		{
			id_i = k*9;
			for (j = 0; j<lenarr_xik_idx; j++)
			{
				id_xik = arr_xik_idx[j] - 1;
				id_jc = j + c;

				ceres::CostFunction* tar_point_cost_function =
					projectionPointError::Create(point_observations[id_jc],
						point_observations[id_jc + observationsm_db2],
						point_normaliser[id_xik]);

				problem.AddResidualBlock(tar_point_cost_function,
					NULL /* squared loss */,
					ba_params + id_i,
					ba_params + num_Hparams + (id_xik * 2));
			}
			
			for (j = 0; j < lenarr_lines_xik_idx; j++)
			{
				id_lines_xik = arr_lines_xik_idx[j] - 1;
				line_id_jc = j+line_c;
				ceres::CostFunction* tar_line_cost_function =
					projectionLineError::Create(line_observations[line_id_jc],
						line_observations[line_id_jc + lenarr_lines_xik_idx],
						line_observations[line_id_jc + 2*lenarr_lines_xik_idx],
						line_normaliser[id_lines_xik]);

				problem.AddResidualBlock(tar_line_cost_function, NULL, ba_params + id_i, 
					ba_params + num_Hparams + 2*num_points + (4*id_lines_xik));
			}
			k++;	
		}
		c += lenarr_xik_idx;
        line_c += 3*lenarr_lines_xik_idx;
	}

	// Make Ceres automatically detect the bundle structure. Note that the
	// standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
	// for standard bundle adjustment problems.
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	//options.minimizer_progress_to_stdout = true;
	options.max_num_iterations = 40;

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	/*std::cout << summary.FullReport() << "\n";    */

	memcpy(finalba_params, ba_params, ba_paramsm * sizeof(double));
	return;
}