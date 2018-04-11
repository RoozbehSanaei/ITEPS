if (b_h < R+0.01) {
	result[0] = 0;
	return;
}


// Electric Motor Model


double motor_cont_mpower = prop_motor_torque_rated*prop_motor_rpm_rated/60*2*pi;
double motor_elect_in  = motor_rate_current*motor_voltage;
double eta_motor = motor_cont_mpower/motor_elect_in;
double num_cell = batt_series*batt_parallel;
double battpack_cap = cell_cap*batt_parallel;


double battpack_volt;
if ((num_cell==1) || (batt_parallel==0)) 
	battpack_volt = cell_volt;
else 
battpack_volt = batt_series * cell_volt;


double Battpack_weight = num_cell*cell_weight;
double Rt = 1/battpack_cap;

double U_inf;

if (condition == 1.0)
	U_inf = 13.6; 
else 
	U_inf = 14.3; 



double h = 1524;
double T_p = 0;
double T_0 = 288 + T_p - (6.5*h/1000);
double P = 101325*pow((1-0.0065*h/288),5.2561);
double rho = P/(287*T_0);
double sigma = rho/1.22500;
double nu = 0.00001646;
double mach = U_inf/sqrt(1.4*287*T_0);
 

double M_tip_lim = 0.7;
double b_b_ratio = 10;
double b_b_eta = 0.98;
double ratio_gear = 1;
double U_dive = 20;
double altitude_ft= h*3.28;
double u_inf_ft = U_inf*3.28;

double ru_ft;
if (h<10000)
	ru_ft = (1-0.00002615*altitude_ft)*0.076474;
else if (h>40000) 
	ru_ft = (-0.0000091*altitude_ft+0.6211)*0.076474;
else 	
	ru_ft = (-0.00001681*altitude_ft+0.9066)*0.076474;



int n_blades = 2;

double dia_m;
if (n_engines == 1)
		dia_m = 0.354; 
if (n_engines == 2)
		dia_m = 0.305; 
if (n_engines == 4)		
		dia_m = 0.203; 
if (n_engines == 6)		
		dia_m = 0.127; 
if (n_engines == 8)		
		dia_m = 0.081; 



// Propeller Model (Single)

double d_prop_ft = dia_m * 3.28;
double RPM = ratio_gear * prop_motor_rpm_rated;
double n = RPM / 60;
double J_advanced_ratio = u_inf_ft / (d_prop_ft * n);
double eta_propeller = (-0.7084*pow(J_advanced_ratio,4)) + (0.5461*pow(J_advanced_ratio,3)) - (0.8019*pow(J_advanced_ratio,2)) + (1.678*J_advanced_ratio) - 0.0117;
double Power_absorbed_kW = motor_cont_mpower / 1000 * eta_propeller;
double Power_absorbed_hp = Power_absorbed_kW * 1000 / 745.7;
double V_tip_m = (dia_m/2)*(n*2*pi);
double V_tip_ft = V_tip_m * 3.28;
double M_tip = V_tip_ft/(1036-0.0034*(altitude_ft-20000));
double C_P = 550*Power_absorbed_hp*32.2/(ru_ft*pow(n,3)*pow(d_prop_ft,5));
double h_P = (0.03*(3-n_blades)+1)*eta_propeller; //dont exceed 6;
double C_T = h_P*C_P/J_advanced_ratio;
double Thrust_Power = C_T/C_P;
double T_static = Thrust_Power *550*Power_absorbed_hp/(n*d_prop_ft);
double T_static_N = T_static*4.45;
double T_static_corrected = (0.03*(n_blades-3)+1)*T_static;
double T_static_corrected_N = T_static_corrected*4.45;
double T_cruise =(Thrust_Power)*550*Power_absorbed_hp/(n*d_prop_ft);
double T_cruise_N = T_cruise*4.45l;
double propulsive_power_kw = T_cruise_N*U_inf / 1000;


// Nacelle pylon drag model
double R_wf = 1;
double L_p = 1.2;
double pylon_sweep = 15;
double pylon_cf = 1.2;
double pylon_hf = 1.2;
double L_nacelle =  L_motor * 1.1;
double D_nacelle = D_motor * 1.1;
double sweep_parameter = cos(pylon_sweep/180*pi);
double nacelle_length_to_diameter = L_nacelle/D_nacelle;
double S_wet_nacelle = L_nacelle*D_nacelle*pi;
double Re_nacelle =  rho * U_inf * L_nacelle / nu;
double C_ff_nacelle = 0.0074;
double C_ff_pylon = 0.0074;
double C_D_nacelle = 0.0684;
double c_pylon = pylon_cf * L_nacelle;
double h_pylon = pylon_hf * D_nacelle;
double S_pylon  = c_pylon*h_pylon;
double S_wet_pylon = 2*S_pylon*1.05;
double ratio_of_areas = S_wet_nacelle/S_pylon;

int no_of_engines = 16;
double buck_boost_efficiency = 0.98;

double Cff_nacelle;

if (mach<0.3) {
		Cff_nacelle = 0.455/(pow(log10(Re_nacelle),2.58)*pow(1+0.144*pow(mach,2),0.58));
	}
else {
		Cff_nacelle = 0.455/(pow(log10(Re_nacelle),2.58)*pow(1+0.144*pow(mach,2),0.58)*(1+(0.35/nacelle_length_to_diameter)));
}


double Cff_pylon;

if (mach<0.3) {
		Cff_pylon = 0.455/(pow(log10(Re_nacelle),2.58)*pow(1+(0.144*pow(mach,2)),0.58));
	}
else {
		Cff_pylon = 0.455/(pow(log10(Re_nacelle),2.58)*pow(1+(0.144*pow(mach,2)),0.58))*(1-(0.09*pow(mach,2)));
}


double CD_nacelle = R_wf*Cff_nacelle*(1+(60/pow(nacelle_length_to_diameter,3))+(0.0025*nacelle_length_to_diameter))*ratio_of_areas;
double max_thickness_to_chord = h_pylon/c_pylon;

double R_LS; int ERR = 0;
if (mach<0.3) 
		R_LS = 1.08;
	else if (mach<0.6)
		R_LS = 1.15;
	else if (mach<0.8)
		R_LS = 1.28;
	else if (mach<0.9)
		R_LS = 1.35;
	else
		ERR = 1;



double CD_pylon = 1*R_LS*Cff_pylon*(1+(L_p*max_thickness_to_chord)+(100*(pow(max_thickness_to_chord,4))))*(S_wet_pylon/S_pylon);
double D_pylon = 0.5 * rho * pow(U_inf,2) * Cff_pylon * S_pylon;
double Dnacelle = 0.5 * rho * pow(U_inf,2) * CD_nacelle * S_wet_nacelle;
double D_pylon_nacelle_s = D_pylon+Dnacelle;
double D_pylon_nacelle_c = D_pylon_nacelle_s*n_engines;
double S_wet_nacelle_pylon = n_engines * (S_wet_nacelle + S_wet_pylon);

// Geometric Sizing
double n_margin = 0.15;
double b_f = 2* b_h;
double t_c = 0.13;
double Ls = 2*c_root;
double n_ult = 3;
double eta_plan = 0.99;
double f_ind = 1 - 1.556*pow(R/b_f,2);
double b_w = b_h - R/2;

double Lambda_1_4 = 0;
double Lambda_1_2 = 0;
double LE = b_h/sin(pi/2-Lambda_1_4);
double S_w = (c_root+c_tip)*(b_h-R);
double AR = pow(b_f,2)/S_w;
double Lambda = c_tip/c_root;

if (AR<=0) {
	result[0] = 0;
	return;
}


double CL_alpha_2D_RAD = 2*pi*AR/(2+sqrt(pow(AR,2)+4));
double k_w = CL_alpha_2D_RAD/(2*pi);
double beta = sqrt(1-pow(mach,2));
double CL_alpha_3D_RAD = (2*pi*AR)/(2 + sqrt((pow(AR,2)*pow(beta,2)/pow(k_w,2)*(1 + tan(pow(Lambda_1_2,2)/pow(beta,2))) + 4)));
double CL_alpha_3D_DEG =  CL_alpha_3D_RAD * pi / 180;
double CL_0 = 0.75;

double CL_F;
if (mach<0.3) {
CL_F=CL_0+CL_alpha_3D_DEG*alpha;
}
else {
CL_F=(CL_0+CL_alpha_3D_DEG*alpha)/beta;
}


double c_avg = (c_root + c_tip)/2;
double MAC = (2 * pow(c_avg,2)  * b_f / 2) / S_w;
double S_HT_f = 0.1;
double S_HT = S_HT_f * S_w;
double c_root_HT = 0.6 * c_root;
double c_tip_HT = 0.8 * c_root_HT;
double b_HT_1_2 = S_HT / ((c_tip_HT + c_root_HT) / 2);
double AR_HT = 4 * pow(b_HT_1_2,2) / S_HT;
double CL_alpha_2D_RAD_HT = 2 * pi * AR_HT / (2 + sqrt(pow(AR_HT,2) + 4));
double k_HT = CL_alpha_2D_RAD_HT / (2 * pi);
double CL_alpha_3D_RAD_HT = (2 * pi * AR_HT) / (2 + sqrt((pow(AR_HT * beta / k_HT,2) * (1 + tan(pow(Lambda_1_2 / beta,2))) + 4)));
double CL_alpha_3D_DEG_HT = (2 * pi * AR_HT) / (2 + sqrt((pow(AR_HT * beta / k_HT,2) * (1 + tan(pow(Lambda_1_2 / beta,2))) + 4))) * pi / 180;
double CL_0_HT = 0.75;
double LE_HT = b_HT_1_2 / sin(pi / 2 - Lambda_1_2);
double n_batt_pack = 1; // Battery Pack Number
double w_batt_total = Battpack_weight * n_batt_pack;
double w_fuel = w_batt_total;
double w_util = 0.1 * w_fuel;
double ZFW_design = MDTW - w_fuel;

double w_wing;
if (MDTW > 10){
	w_wing = ((4.22 * S_w * 10.76) + ((0.000001642 * n_ult * pow(b_f,3) * sqrt(MDTW * ZFW_design) * (1 + (2 * Lambda))) / ((t_c * cos(Lambda_1_4) * cos(Lambda_1_4) * S_w * 10.76) * (1 + Lambda)))) * 0.454;
}
else {
	w_wing = ((0.000212 * (MDTW * 2.2)) * n_ult * ((b_h * 2) + 33)) * 0.454;
}



double w_propulsion_sys = n_engines * W_motor;
double w_eng = w_propulsion_sys;

double w_body = ZFW_design - w_wing - w_eng;
double lambda_VT = 0.3;
double S_VT_f = 0.1;
double S_VT = S_VT_f * S_w;
double c_root_VT = 0.9 * c_root;
double c_tip_VT = lambda_VT * c_root_VT;
double b_VT_2 = 2 * S_VT / (c_root_VT + c_tip_VT);

double P_F = 0;
double B_F = 2 * R;
double H_F = B_F;
double n_limit = n_ult - n_margin;
double Ip = 0.0015 * (0.021 * P_F) * (3.28 * B_F);
double Ib = (0.000191 * n_limit * w_body * 2.2 * L_f * 3.28) / pow(H_F * 3.28, 2);
double S_fuse = pi * L_f * B_F;
double a2 = 0.2;
double a1 = L_f - c_root - Ls - c_root_HT - a2;
double h0 = (0.25 * MAC) + a1;
double eta_tail = 0.9;
double V_stab_coeff = (S_HT * Ls) / (S_w * c_avg);
double as_aw = CL_alpha_3D_RAD_HT / CL_alpha_3D_RAD;
double de_alpha = (2 * CL_alpha_3D_RAD / (pi * AR));
double x_np = h0 + eta_tail * V_stab_coeff * as_aw * (1 - de_alpha);

double x_wing = a1 + MAC / 2;
double x_HT = a1 + c_root + Ls + c_root_HT / 2;
double x_VT = a1 + c_root + Ls + c_root_VT / 2;
double x_fuse = L_f / 2;

double x_eng = a1 + c_root + Ls;

double engine_position,x_engine;

if (n_engines == 1) {
	engine_position = 1;	//front
	x_engine = a1 / 2;
}
else if (n_engines > 1) {
	engine_position = 2; //wing
	x_engine = a1 + 0.25 * MAC;
};





double x_util = 0.2 * L_f;
double x_payload = 1.3 * a1;
double x_fuel = 0.7 * L_f;
double S_wet_wing = 2 * (1.9767 + 0.553 * t_c) * S_w;
double S_wet = 2 * S_HT + 2 * S_VT + S_wet_wing + (2 * pi * R) / 4 + 2 * pi * R * L_f + S_wet_nacelle_pylon;
double S_wet_S = S_wet / S_w;
double RE = ((S_wet / b_h) * U_inf) / nu;
double cfe = 0.00258 + (0.00102 * exp(-0.00000000628 * RE)) + (0.00295 * exp(-0.0000000201 * RE));
double CD_0 = cfe * S_wet / S_w;
double k_vis = 0.38 + 0.0000005 * pow(Lambda_1_4 * pi / 180,2) * CD_0;

double e;
if (Lambda_1_4==0) {
	e=1.78*(1-0.045*pow(AR,0.68))-0.64;
	}
else {
	e=1/(pi*k_vis*AR+1/(eta_plan*f_ind));
}

double K = 1/(pi*e*AR);

double CD_i = pow(CL_F,2) / (pi * AR * e);
double CD_tot = CD_0 + CD_i;
double D_BL = 0.5 * rho * pow(U_inf, 2) * S_w * CD_tot;

double I_fuse;
if (Ip > Ib)
	I_fuse = Ib;
else
	I_fuse = (pow(Ip, 2) + pow(Ib, 2)) / (2 * Ib);



double w_fuse,w_payload,w_HT,w_VT;
if (MDTW > 10) {
		w_fuse = ((1.051 + (0.102 * I_fuse)) * (S_fuse * 10.76)) * 0.454;
		w_payload = (1.4950 * r + 0.805) * w_fuse;
		w_HT = (5.25 * 0.25 * S_HT * 10.76) + ((0.0000008 * n_ult * pow((b_HT_1_2 * 2) , 3) * (MDTW * 2.2) * sqrt(0.25 * S_HT * 10.76)) / (t_c * cos(Lambda_1_4) * cos(Lambda_1_4) * (c_root_HT * 3.28) * (pow(S_HT,1.5)))) * 0.454;
		w_VT = ((2.62 * S_VT * 10.76) + ((0.000015 * n_ult * pow((b_VT_2 * 2), 3) * (8 + (0.44 * (MDTW * 2.2) / (S_w * 10.76)))) / (L_f / b_f * cos(Lambda_1_4) * cos(Lambda_1_4)))) * 1.15 * 0.454;
	} else {
		w_fuse = 0.23 * pow(S_fuse, 1.2) * sqrt((U_dive * L_f) / (2 * 2 * R));
		w_payload = (1.53 * r + 0.92) * w_fuse;
		w_HT = ((0.000212 * (MDTW * 2.2)) * n_ult * ((b_HT_1_2 * 2) + 33)) * 0.454;
		w_VT = 1.25 * ((0.000212 * (MDTW * 2.2)) * n_ult * ((b_HT_1_2 * 2) + 33)) * 0.454 * 1.15;
	}

//Moment Arm and CG
double ZFW_actual = w_wing + w_fuse + w_HT + w_VT + w_eng + w_util + w_payload;
double M_wing = w_wing * x_wing;
double M_HT = w_HT * x_HT;
double M_VT = w_VT * x_VT;
double M_eng = w_eng * x_eng;
double M_fuse = w_fuse * x_fuse;
double M_util = w_util * x_util;
double M_payload = w_payload * x_payload;
double M_fuel = w_fuel * x_fuel;
double M_total = M_wing + M_HT + M_VT + M_eng + M_payload + M_fuel + M_fuse + M_util;

double x_cg = M_total / MDTW;
double V_h = S_HT * (x_HT - x_cg) / (S_w * c_avg);
double V_v = S_VT * (x_VT - x_cg) / (S_w * b_f);

double L_sys = 0.5 * rho * pow(U_inf, 2) * S_w * CL_F;
double MTOW = ZFW_actual + w_fuel;
double L_check = L_sys - MTOW * 9.81;
double lift_to_drag = L_sys / D_BL;
double W_b = MTOW - (1 - 0.9) * w_fuel;
double W_a = MTOW - (1 - 0.5) * w_fuel;
double PR_kw = ((0.5 * rho * pow(U_inf,3) * S_w * CD_0) + ((2 * pow(MTOW * 9.81 , 2)) * K) / (rho * U_inf * S_w)) / 1000;

double PR_sys;

if (condition==1) 
		PR_sys = 0.01; 
if (condition==2)
		PR_sys = 0.02;


	double PR_kw_act = PR_kw + PR_sys;
	double mach_limit = 0.7;

	double thurst_required = 1000 * PR_kw / U_inf;
	double propulsive_thrust = 69.7;
	double mach_tip = (-0.0021 * pow(PR_kw, 2) + (0.0948 * PR_kw) + 0.2373) * (((rho / 1.226) - 0.117) / 0.883);

	//electrical propulsion model
	double current_demand_prop_sys = motor_rate_current * batt_parallel * n_engines;
	double batt_sys_discharge = current_demand_prop_sys / battpack_cap;
	double batt_sys_volt = battpack_volt * b_b_ratio;
	double PR_per_split = PR_kw_act / n_batt_pack;
	double P_total = (n_engines * motor_cont_mpower) / 1000;
	double total_propulsive_power = n_engines * propulsive_power_kw;
	double propulsive_thrust_total = T_cruise_N * n_engines;
	double eta_overall = eta_propeller * eta_motor * 0.98;
	double max_no_engines = 0.7 * 2 * b_h / dia_m;


	double lift_margin = (MDTW - MTOW) / MDTW;
	double engine_power_margin = (total_propulsive_power - PR_kw) / PR_kw;
	double thrust_drag_margin = ((propulsive_thrust - D_BL) / D_BL);
	double wing_span_margin = 5 - b_f;
	double propeller_tip_mach_margin = mach_limit - mach_tip;
	double wing_load = MTOW / S_w;
	double static_margin = (x_np - x_cg) / c_avg;
	double lcm = (L_sys / 9.81 - MTOW) / MTOW;
	double TM = min(engine_power_margin,lift_margin);


	double endurance = pow(Rt, (1 - cell_constant)) * pow((eta_overall * cell_volt * battpack_cap) / (PR_per_split * 1000), cell_constant);
	double range = U_inf * 3600 * endurance / 1000;

	bool constraints[10] = {L_check>0,
	batt_sys_discharge<max_discharge,
	batt_sys_volt>motor_voltage,
	total_propulsive_power<P_total,
	propulsive_thrust_total>D_BL,
	(eta_overall>0) && (eta_overall<1),
	lift_margin>0,
	engine_power_margin>0,
	propeller_tip_mach_margin>0,
	n_engines<max_no_engines,
	};
