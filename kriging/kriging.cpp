// kriging.cpp : Defines the entry point for the console application.
//

using namespace std;

#include <vector>
#include <iostream>

typedef vector<double> VecDoub;
//typedef vector<VecDoub> VecDoub;

#define MIN(a,b) ((a>b)?b:a)
#define MAX(a,b) ((a>b)?a:b)
#define CAP_GRID_SIZE (1000)
#define CAP_GRID_CURVE (1.0)
#define CAP_GRID_MAX (100)
#define CAP_GRID_MIN (0.00001)
#define TOLERANCE (0.00001)
#define BETA (0.95)
#define R (0.99/BETA-1)
#define NUM_STATES 4
#define WAGE 0.5
#define MAXIT 1000

const double TRANSITION[NUM_STATES][NUM_STATES] = { {0.2,0.4,0.3,0.1},{ 0.2,0.4,0.3,0.1 },{ 0.2,0.4,0.3,0.1 },{ 0.2,0.4,0.3,0.1 } };

void createCapitalGrid(VecDoub& capGrid) {
	/*
	!this routine constructs a grid for capital.The cur_a_v parameter defines what curvature
	!to place on the grid, where higher values place more at the left.
	!
	!INPUT: capGrid
	!OUTPUT : NONE
	!USES :
	!EFFECTS :
		!capGrid : the capital level grid
	*/
	capGrid.resize(CAP_GRID_SIZE);
	for (int i = 0; i < CAP_GRID_SIZE; i++) {
		capGrid[i] = 1.0 / (CAP_GRID_SIZE-1) * i;
	}
	for (int i = 0; i < CAP_GRID_SIZE; i++) {
		capGrid[i] = pow(capGrid[i], CAP_GRID_CURVE);
		capGrid[i] = capGrid[i] * (CAP_GRID_MAX - CAP_GRID_MIN) + CAP_GRID_MIN;
	}
	//capGrid[CAP_GRID_SIZE - 1] = CAP_GRID_MAX;
}
double inv_utility(double v_prime) {
	return pow(v_prime,-0.5);
}

double utility(double consumption) {
	return pow(consumption,1-2)/(1-2);
}

double interpolate(const VecDoub &x, const VecDoub::const_iterator y, double a_prime) {

	int right = 0;
	int left = 0;
	double slope;
	double distance;
	double interp_val;
	int size = x.size();

	if (size == 1) {
		std::cerr << "UtilityFunctions.cpp-interpolate():Should never pass a size 1 vector to interpolate" << std::endl << std::flush;
		exit(-1);
	}

	if (a_prime <= x[0]) {
		right = 1;
		left = 0;
	}
	else if (a_prime >= x[size - 1]) {
		right = size - 1;
		left = right - 1;
	}
	else {
		int i;
		for (i = 0; i < size; i++) {
			if (x[i] >= a_prime) {
				right = i;
				left = i - 1;
				break;
			}
		}
		if (i == size) {
			cout << "interpolate: Error! Couldn't find aprime= " << a_prime << endl;
			for (i = 0; i < size; i++) {
				cout << x[i] << ",";
			}
			exit(-1);
		}
	}
	slope = (y[right] - y[left]) / (x[right] - x[left]);
	distance = a_prime - x[left];
	interp_val = y[left] + distance * slope;

	if (interp_val != interp_val) {
		std::cout << "ERROR! utilityFunctions-interpolate() - return value is NaN. "
			<< y[right] << "-" << y[left] << "/(" << x[right] << "-"
			<< x[left] << "), distance=" << distance << " slope=" << slope
			<< " a_prime=" << a_prime << std::endl << std::flush;
		std::cout << "left: " << left << "  right: " << right << "   size: " << size << std::endl;
		exit(-1);
	}

	return interp_val;
}

void consumerProblemEDG(vector<VecDoub>& valFn, vector<VecDoub>& polFn, VecDoub& capGrid)
{
	/*
	!this subroutine is based on the pseudocode on page 2703 - 2704 of "A generalization of the endogenous grid method" by
	!Francisco Barillas and Jesus Fernandez - Villaverde, JEDC(2006)
	!
	!http ://economics.sas.upenn.edu/~jesusfv/Endogenous_Grid.pdf
	!
	!Mapping of variables(Barillas->code)
	!Gk(t + 1)->a_v
	!V - tilde->value_v
	!
	!INPUT:
		!firstCall : On the first call, need to initialize the value function such that it is increasing in k
	!OUTPUT : NONE
	!EFFECTS :
		!value_v_t - at the end of the subroutine, this array is set to the current estimate of the value function
		!policy_v_t - at the end of the subroutine, this array is set to the current estimate of the policy function
	*/
	vector<VecDoub> temp_v_n(NUM_STATES),temp_k(NUM_STATES), temp_v_np1(NUM_STATES), temp_vk(NUM_STATES),
		c_star(NUM_STATES), y_end(NUM_STATES), G_y(NUM_STATES);

	double tmp, vLeft, vRight;

	//Every time subroutine is called, make sure value_v is increasing in k.If not, convergence
	//is highly unlikely
	for (int state = 0; state < NUM_STATES; state++) {
		temp_v_n[state] = VecDoub(CAP_GRID_SIZE,0);
		temp_k[state] = VecDoub(CAP_GRID_SIZE, 0);
		temp_v_np1[state] = VecDoub(CAP_GRID_SIZE, 0);
		temp_vk[state] = VecDoub(CAP_GRID_SIZE, 0);
		c_star[state] = VecDoub(CAP_GRID_SIZE, 0);
		y_end[state] = VecDoub(CAP_GRID_SIZE, 0);
		G_y[state] = VecDoub(CAP_GRID_SIZE, 0);

		for (int i = 1; i < CAP_GRID_SIZE; i++) {
			if (valFn[state][i] <= valFn[state][i-1]) {
				cout << "error: ConsumerProblemEDG: value function must be increasing in size" << endl;
				exit(-1);
			}
		}
	}

	//preliminary setup
	for (int state = 0; state < NUM_STATES; state++) {
		for (int i = 0; i < CAP_GRID_SIZE; i++) {
			G_y[state][i] = capGrid[i] * (1 + R) + WAGE*state;
		}
	}

	for (int state = 0; state < NUM_STATES; state++) {
		//STEP 1: get derivative of V - tilde
		for (int i = 0; i < CAP_GRID_SIZE; i++) {
			if (i == 0) {
				tmp = capGrid[1] - capGrid[0];
			}else {
				tmp = capGrid[i] - capGrid[i - 1];
			}
			vLeft = interpolate(capGrid, valFn[state].begin(), capGrid[i] - tmp / 2);
			vRight = interpolate(capGrid, valFn[state].begin(), capGrid[i] + tmp / 2);
			temp_vk[state][i] = (vRight - vLeft) / (tmp);

			if (temp_vk[state][i] < 0) {
				cout << "error: negative slope. Should not happen." << endl;
				exit(-1);
			}

			if (temp_vk[state][i] != temp_vk[state][i]) {
				cout << "error: slope=NaN" << endl;
				exit(-1);
			}

#if 0
			if (i > 0) {
				if (temp_vk[state][i] > temp_vk[state][i-1]) {
					cout << "error: slope should be decreasing for point " << i << endl;
					cout << temp_vk[state][i - 1] << ":" << temp_vk[state][i] << ":" << temp_vk[state][i]-temp_vk[state][i-1]<<endl;
					cout << vLeft << ":" << vRight << ":" << tmp << endl;
					exit(-1);
				}
			}
#endif
		}
		//STEP 2: calculate C*, Yend, and Vn
		for (int i = 0; i < CAP_GRID_SIZE; i++) {
			c_star[state][i] = MIN(inv_utility(temp_vk[state][i]),capGrid[i]/(1+R)+WAGE*state+CAP_GRID_MIN);
			y_end[state][i] = c_star[state][i] + capGrid[i];
			temp_v_n[state][i] = utility(c_star[state][i]) + valFn[state][i];
#if 1
			if (i > 0) {
				if (temp_v_n[state][i] < temp_v_n[state][i - 1]) {
					cout << "error: temp_v_n must be increasing in k " << state << ":" << i << ":" << temp_v_n[state][i - 1] << ":" << temp_v_n[state][i] << endl;
					cout << c_star[state][i - 1] << ":" << c_star[state][i] << endl;
					cout << temp_vk[state][i - 1] << ":" << temp_vk[state][i] << endl;
					exit(-1);
				}
			}
#endif
		}
	}

	//STEP 3: Calculate V_n + 1.
	for (int state = 0; state < NUM_STATES; state++) {
		for (int i = 0; i < CAP_GRID_SIZE; i++) {
			temp_v_np1[state][i] = interpolate(y_end[state], temp_v_n[state].begin(), G_y[state][i]);

			/*
			!an error check.Basically, we need to make sure that the value function is increasing.If
			!this is not true, then often we don't get convergence, so may as well kill the program immediately
			*/
			if (i > 0) {
				if (temp_v_np1[state][i] < temp_v_np1[state][i - 1]) {
					cout << "error: temp_v_np1 must be increasing in k " << state << ":" << i << ":" << temp_v_np1[state][i-1] << ":" << temp_v_np1[state][i] << endl;
					exit(-1);
				}
			}
		}
	}

	//STEP 4: V_tilde
	for (int state = 0; state < NUM_STATES; state++) {
		for (int i = 0; i < CAP_GRID_SIZE; i++) {
			tmp = 0;
			for (int i2 = 0; i2 < NUM_STATES; i2++) {
				tmp += TRANSITION[state][i2] * temp_v_np1[i2][i];
			}
			valFn[state][i] = BETA*tmp;
			temp_k[state][i] = (y_end[state][i] - (state*WAGE)) / (1 + R);
		}
		polFn[state] = temp_k[state];
	}

	/*
	!penalty function for going below 0. Have found this necessary to get stronger convergence of the value function.
	!Without this, ripples form and prevent convergence in many sceanrios.
	*/
	/*
			DO i = 1, dim_a_v
			tmp = value_v_t(j, i) - &
			max(0.0_dp, -policy_v_t(j, i))**2.0_dp
			value_v_t(j, i) = tmp
			END DO
	*/
}

void solveConsumerProblem(vector<VecDoub>& valueFn, vector<VecDoub>& policyFn, VecDoub& capGrid)
{

	if (valueFn.size() != policyFn.size()) {
		cout << "Value function and policy function should be the same size." << endl;
		exit(-1);
	}

	if (valueFn[0].size() != policyFn[0].size()) {
		cout << "Value function in state 0 and policy function in state 0 should be the same size." << endl;
		exit(-1);
	}

	if (valueFn[0].size() != CAP_GRID_SIZE) {
		cout << "Expect value/policy function grid size of " << CAP_GRID_SIZE << endl;
		exit(-1);
	}

	for (int state = 0; state < NUM_STATES; state++) {
		for (int i = 0; i < CAP_GRID_SIZE; i++) {
			valueFn[state][i] = utility(capGrid[i] + WAGE*state +0.001);
			policyFn[state][i] = 0;
		}
	}

	vector<vector<VecDoub>> history(NUM_STATES);
	for (int i = 0; i < NUM_STATES; i++) {
		history[i].resize(CAP_GRID_SIZE);
	}
	/*
		!this is the fundamental loop solving the consumer problem.The loop is ended when either the supnorm or average
		!(depending on the invocation options) is less than the tolerance.Each iteraction solves one iteration of the
		!consumer problem
	*/
	for (int maximizationCount = 0; maximizationCount < MAXIT; maximizationCount++) {
		double diff = -100, minDiff = 100;

		//save the previous policy and value functions so that we can calculate the difference
		vector<VecDoub> policy_last(NUM_STATES), value_last(NUM_STATES);
		for (int state = 0; state < NUM_STATES; state++) {
			policy_last[state] = VecDoub(policyFn[state]);
			value_last[state] = VecDoub(valueFn[state]);
		}

		consumerProblemEDG(valueFn, policyFn, capGrid);

		//calculate the difference
		for (int state = 0; state < NUM_STATES; state++) {
			for (int i = 0; i < CAP_GRID_SIZE; i++) {
				double temp = abs(valueFn[state][i] - value_last[state][i]);
				diff = MAX(diff, temp);
				minDiff = MIN(minDiff, temp);
				history[state][i].push_back(valueFn[state][i]);
			}
		}

		/* Adjust formatting before print*/
		ios oldState(nullptr);
		oldState.copyfmt(cout);
		cout.precision(10);
		cout << "Iteration " << maximizationCount << " Diff = " << diff << endl;
		cout.copyfmt(oldState);

		if (diff < TOLERANCE) {
			break;
		}
	}
#if 0
	exit(-1);
	for (int i = 0; i < NUM_STATES; i++) {
		cout << i<<",";
		for (int j = 0; j < CAP_GRID_SIZE; j++) {
			cout << j << ",";
			for (int k = 0; k < history[i][j].size(); k++) {
				cout << history[i][j][k] << ",";
			}
			cout << endl;
		}
	}
#endif
}

int main()
{
	VecDoub capGrid;
	vector<VecDoub> valueFunction(NUM_STATES), policyFunction(NUM_STATES);
	
	for (int i = 0; i < NUM_STATES; i++) {
		valueFunction[i] = VecDoub(CAP_GRID_SIZE);
		policyFunction[i] = VecDoub(CAP_GRID_SIZE);
	}
	createCapitalGrid(capGrid);
	solveConsumerProblem(valueFunction, policyFunction, capGrid);

#if 1	
	for (int i = 0; i < CAP_GRID_SIZE; i++) {
		cout << capGrid[i] << ",";
		for (int j = 0; j < NUM_STATES; j++ ) {
			cout << valueFunction[j][i];
			if (j < NUM_STATES - 1) {
				cout << ",";
			}
		}
		cout << endl;
	}

	for (int i = 0; i < CAP_GRID_SIZE; i++) {
		cout << capGrid[i] << ",";
		for (int j = 0; j < NUM_STATES; j++ ) {
			cout << policyFunction[j][i];
			if (j < NUM_STATES - 1) {
				cout << ",";
			}
		}
		cout << endl;
	}
#endif	

    return 0;
}

