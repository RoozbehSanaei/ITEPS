
//Camera
double payload_1 = camera_weight;


//Motor
double elect_powr_loss = pow(motor_rate_current,2)*Int_resist;
double motor_elect_in  = motor_rate_current*motor_voltage;
double prop_motor_rpm_rated = motor_rpm_per_volt*motor_voltage;
double prop_motor_rps_rated = prop_motor_rpm_rated/60;
double prop_motor_omega_rated = prop_motor_rps_rated*2*pi;
double prop_motor_torque_rated = motor_torque_per_volt*motor_rate_current;
double motor_cont_mpower = prop_motor_omega_rated*prop_motor_torque_rated;
double eta_motor_elect = (motor_elect_in-elect_powr_loss)/motor_elect_in;
double eta_motor_mech = motor_cont_mpower/(motor_elect_in-elect_powr_loss);
double eta_motor = motor_cont_mpower/motor_elect_in;
double eta_propeller = 0.8;
int engine_ON = n_engines;
double motor_rate_current_total = motor_rate_current*engine_ON;
double total_current_drawn = motor_rate_current_total+airborne_systems_current;

//Battery
double cell_cap_rated = cell_cap*cell_cap_factor;
double num_cell = batt_series*batt_parallel;
double battpack_cap = cell_cap*batt_parallel;

double battpack_volt;
if ((num_cell==1) || (batt_parallel==0)) 
battpack_volt = cell_volt;
else 
battpack_volt = batt_series * cell_volt;
double battpack_weight = num_cell*cell_weight;
double Rt = 1/battpack_cap;
double Rt_rated = 1/cell_cap_rated;


//Endurance
double end_hr = (Rt/pow(total_current_drawn,cell_constant))*pow(battpack_cap/Rt,cell_constant);
double end_min = end_hr*60;
double end_hr_rated = (Rt_rated/pow(total_current_drawn,cell_constant))*pow(battpack_cap/Rt_rated,cell_constant);
double end_min_rated = end_hr_rated*60;



double altitude_m  = 0;
double altitude_ft = 0;
double axial_speed = 0;


double h = 1524;
double T_p = 0;
double T_0 = 288 + T_p - (6.5*h/1000);
double P = 101325*pow((1-0.0065*h/288),5.2561);
double rho = P/(287*T_0);
double sigma = rho/1.22500;
double nu = 0.00001646;
///double mach = U_inf/sqrt(1.4*287*T_0);'''

double ru_ft;

if (h<10000)
ru_ft = (1-0.00002615*altitude_ft)*0.076474;
else if (h>40000) 
ru_ft = (-0.0000091*altitude_ft+0.6211)*0.076474;
else 	
ru_ft = (-0.00001681*altitude_ft+0.9066)*0.076474;


double M_tip_lim = 0.7;///input
double b_b_ratio = 10;
double b_b_eta = 0.98;
double ratio_gear = 1;
double mach_limit = 0.8;






int n_blades = 2;
double dia_m = dia_inch*0.0254;


double d_prop_ft = dia_m * 3.28;
double RPM = ratio_gear * prop_motor_rpm_rated;
double n = RPM / 60;
double V_tip_m = (dia_m/2)*(n*2*pi);
double V_tip_ft = V_tip_m * 3.28;
double mach_tip = V_tip_ft/(1036-0.0034*(altitude_ft-20000));


double n_batt_pack = 1;
double w_batt_total = battpack_weight * n_batt_pack;
double w_fuel = w_batt_total;
double w_util = 0.1 * w_fuel;



double empty_frame = 0.2* MDTW;
double payload_2 = 0.1*MDTW;

double MTOW = empty_frame + battpack_weight +(n_engines* W_motor )+ payload_1 + payload_2;


double T_static_thrust_SL = 1.225*pi*pow((0.0254*dia_inch),2)/4*(pow(RPM*0.0254*prop_pitch/60,2)-(RPM*0.0254*prop_pitch/60)*axial_speed)*pow(dia_inch/(prop_pitch*3.29546),1.5);



double thrust_required = MTOW * 9.81;
double thrust_generated = n_engines * T_static_thrust_SL;


//electrical propolsion model

double lift_margin = (MDTW - MTOW) / MDTW;
double thrust_margin = (thrust_generated-thrust_required)/thrust_required;
double propeller_tip_margin = (mach_limit - mach_tip)/mach_limit;



double TM = min(lift_margin, thrust_margin, propeller_tip_margin);


double current_demand_prop_sys = motor_rate_current * batt_parallel * n_engines;
double batt_sys_discharge = current_demand_prop_sys / battpack_cap;
double batt_sys_volt = battpack_volt;
double eta_overall = eta_propeller * eta_motor * 0.98;

bool constraints[4] = {
batt_sys_discharge<max_discharge,
batt_sys_volt>motor_voltage,
(eta_overall>0) && (eta_overall<1),
lift_margin>0,
};
