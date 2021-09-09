#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <sciplot/sciplot.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;
using namespace sciplot;


typedef vector<Matrix<double,17,1>> vector_rec;
typedef Matrix<double,8,1> Vector8d;

Vector4d CL_term(0.0,0.0,0.0,0.0);
Matrix<double,4,4> YFYF = MatrixXd::Zero(4,4);
Vector4d           YFtheta_star(0.0,0.0,0.0,0.0);

typedef struct SysPara
{
    Vector4d theta_star;
    double   gamma;
    double   gamma_epsilon;
    double   K;
    double   tau;
} syspara;


typedef VectorXd (* DynamiModel) (const double &t, const VectorXd &x, const syspara &para);

class RungeKutta
{
    private:
        double m_time_step;
        DynamiModel m_dynamic_model;

    public:
        RungeKutta(DynamiModel dynamic_model, double time_step);
        void Euler(double &t, VectorXd &x, const syspara &para);
        void Rk45(double &t, VectorXd &x, const syspara &para);
};

RungeKutta::RungeKutta(DynamiModel dynamic_model, double time_step)
{
    m_dynamic_model = dynamic_model;
    m_time_step     = time_step;
}

void RungeKutta::Euler(double &t, VectorXd &x, const syspara &para)
{
    x = x + m_time_step * m_dynamic_model(t,x,para);
    t = t + m_time_step;
}

void RungeKutta::Rk45(double &t, VectorXd &x, const syspara &para)
{
    VectorXd k1 = m_dynamic_model(t, x, para);
    
    t = t + m_time_step / 2;
    VectorXd x1 = x + k1 * m_time_step / 2;
    VectorXd k2 = m_dynamic_model(t, x1, para);
    
    VectorXd x2 = x + k2 * m_time_step / 2;
    VectorXd k3 = m_dynamic_model(t, x2, para);
    
    t = t + m_time_step / 2;
    VectorXd x3 = x + m_time_step * k3;
    VectorXd k4 = m_dynamic_model(t, x3, para);
    
    x = x + m_time_step / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
}
//-----------------------------------------------------------
//|   x   | theta |   xi  |   YF   |
//-----------------------------------------------------------
//| [0:1] | [2:5] | [6:7] | [8:15] |
//-----------------------------------------------------------
VectorXd NonlinearDynamics(const double &t, const VectorXd &states, const syspara &para)
{
    Vector2d x(states(0),states(1));        // system states
    double x1 = states(0);
    double x2 = states(1);
    Vector4d theta(states(2),states(3),states(4),states(5));
    Vector2d xi(states(6),states(7));
    Vector8d YF_vec(states(8),states(9),states(10),states(11),states(12),states(13),states(14),states(15));
    Vector2d x_d(1.0,0.0);

    Vector4d theta_star = para.theta_star;
    double K = para.K;
    Vector8d Y_vec(x1*x1, sin(x2), 0.0, 0.0,0.0, x2*sin(t), x1, x1*x2);
    Matrix<double,2,4> Y  = Y_vec.reshaped<RowMajor>(2,4);
    Matrix<double,2,4> YF = YF_vec.reshaped<RowMajor>(2,4);
    
    //******************************************************************************************
    Vector2d e = x_d - x;
    Vector2d u = -Y * theta + K * e;
    //******************************************************************************************
    Vector2d YFtheta = xi - e/para.tau;

    if (e.norm() > 0.1)
    {
        YFYF    = YFYF + YF.transpose()*YF;
        YFtheta_star = YFtheta_star + YF.transpose() * YFtheta;
    }
    // YFYF    = 0.8*YFYF + 0.01*(YF.transpose()*YF);
    // YFtheta_star = 0.8*YFtheta_star + 0.01*(YF.transpose() * YFtheta);
    CL_term = YFYF * theta - YFtheta_star;
    //******************************************************************************************
    Vector2d dot_x = Y * theta_star + u;
    Vector4d dot_theta = - para.gamma * Y.transpose() * e - para.gamma_epsilon * CL_term;
    Vector2d dot_xi = (Y * theta + (1/para.tau - K) * e - xi)/para.tau;
    Vector8d dot_YF = (Y_vec - YF_vec)/para.tau;

    VectorXd dot_states(states.size());
    dot_states <<dot_x, dot_theta, dot_xi, dot_YF;
    return dot_states;
}

void DataVisiualization(const MatrixXd &plotdata)
{
    Plot plot_x_2D;
    plot_x_2D.xlabel("time/[s]");
    plot_x_2D.ylabel("x(t)");
    plot_x_2D.drawCurve(plotdata(0,all), plotdata(1,all)).label("x_1").lineWidth(2);
    plot_x_2D.drawCurve(plotdata(0,all), plotdata(2,all)).label("x_2").lineWidth(2);
    plot_x_2D.grid().show();

    Plot plot_theta_2D;
    plot_theta_2D.xlabel("time/[s]");
    plot_theta_2D.ylabel("theta");
    plot_theta_2D.drawCurve(plotdata(0,all), plotdata(3,all)).label("theta_1").lineWidth(2);
    plot_theta_2D.drawCurve(plotdata(0,all), plotdata(4,all)).label("theta_2").lineWidth(2);
    plot_theta_2D.drawCurve(plotdata(0,all), plotdata(5,all)).label("theta_3").lineWidth(2);
    plot_theta_2D.drawCurve(plotdata(0,all), plotdata(6,all)).label("theta_4").lineWidth(2);
    plot_theta_2D.grid().show();


    // Use the previous plots as sub-figures in a larger 2x2 figure.
    Figure fig1 = {{plot_x_2D},{plot_theta_2D}};
    fig1.size(600,600);
    fig1.title("Missile-Type UAV Position & Attitude");
    fig1.palette("set1");
    fig1.show();   
}

int main(int argc, char * argv[])
{
    // output simulation results to specified files
    const string save_path = "../result.csv"; 
    ofstream save_file(save_path);
    save_file << "time,x1,x2,theta1,theta2,theta3,theta4,xi1,xi2,YF1,YF2,YF3,YF4,YF5,YF6,YF7,YF8"<<endl;                  

    syspara config;
    config.theta_star = Vector4d(5.0,10.0,15.0,20.0);
    config.gamma = 100;
    config.gamma_epsilon = 2;
    config.K =  5;
    config.tau = 0.01;
    //******************************************************************************************
    Vector2d x0(0.0,0.0);                           // initial states
    Vector4d theta0(0.0,0.0,0.0,0.0);               // initial estimates
    Vector2d xi0(0.0,0.0);                          // initial xi
    Vector8d YF0(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);  // initial yf
    
    int states_dim = x0.size() + theta0.size() + xi0.size() + YF0.size(); // states dimensions
    
    //******************************************************************************************
    
    constexpr double init_time     = 0.0;        // initial simulation time
    constexpr double final_time    = 1.0;        // final simulation time
    constexpr double time_step     = 0.01;       // simulation time step
    constexpr double data_rec_step = 0.01;        // data recording time step
    constexpr int num_rec_step     = (int)(data_rec_step/time_step);// data recording numbers
    //******************************************************************************************

    RungeKutta cac = RungeKutta(NonlinearDynamics,time_step);  // simulation instance
    double t = init_time;                               // simulation time
    VectorXd states_value(states_dim);                  // initial states
    states_value << x0, theta0, xi0, YF0;

    const bool isplot   = true;                  // gnuplot windows
    const bool isoutput = true;                  // console output

    int iter_num = 0;                            // iteration numbers

    Matrix<double,17,1> states;                 // total system states
    
    cout.precision(6);                          // config the terminal output format
    cout.flags(ios::fixed);
    cout.setf(ios::right);

    vector_rec rec_data;                        // Simulations results vector<double>
    // start loop to simulate
    while (t < final_time)
    {   
        // augment time and states
        states << t, states_value;
        if((iter_num % num_rec_step)==0)
        {
            rec_data.push_back(states);
            for (int i = 0; i < states.size(); i++)
            {
                save_file << states(i) << ',';
                if(isoutput == true)
                {
                    cout << states(i) << '\t';
                }
            }
            save_file << endl;
            if(isoutput == true)
            {
                cout << endl;
            }
        }
        cac.Euler(t,states_value,config);
        iter_num += 1;
    }
    save_file.close();
    // plot
    if(isplot == true)
    {
        MatrixXd plotdata = Map<MatrixXd> (rec_data[0].data(),17,rec_data.size());
        DataVisiualization(plotdata);
    }
    return 0;
}