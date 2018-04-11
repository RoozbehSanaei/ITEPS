import xlrd
import sys;

X = sys.argv[1];
Y = sys.argv[2];
Z = sys.argv[3];

data = xlrd.open_workbook('data.xlsx')
component_names = [c for c in data.sheet_names() if (c!='Inputs' and c!='Outputs')];


def save_to_txt_file(M,filename):
	thefile = open(filename, 'w')
	for item in M:
		thefile.write("%s\n" % item)
	thefile.close();

def load_from_txt_file(filename):
	thefile = open(filename, "r")
	lines = thefile.read().splitlines()
	thefile.close();
	return lines;

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False


def convert(S):
	if is_number(S):
		return str(S)
	else:
		return '"'+str(S)+'"'

def find_all(query,list):
	indices = [i for i, x in enumerate(L) if x == query];
	if len(indices) == 0 :
		return (-1)
	elif len(indices) == 1 :
		return (indices[0])
	else:
		return indices



var_names = dict();
array_defs = [];
var_defs = [];
for sheet_name in data.sheet_names():
	if sheet_name == 'Outputs':
		break;
	var_names[sheet_name] = [];
	inputs_sheet = data.sheet_by_name(sheet_name)
	inputs_vars = [[str(inputs_sheet.cell_value(r, c)) for r in range(0,inputs_sheet.nrows)] for c in range(inputs_sheet.ncols)]
	for input_var_list in inputs_vars:
		if is_number(input_var_list[1]):
			T = 'double';
		else:
			T = 'string';
		
		j = 1;
		while (input_var_list[-j]==''):
			j = j + 1;
	
		L = 'const ' + T + ' '+ input_var_list[0] + '_list[' + str(len(input_var_list)-j)+ '] = {' + convert(input_var_list[1]);
		var_names[sheet_name].append([input_var_list[0],T]);

		for i in range(2,len(input_var_list)-j + 1):
			if (input_var_list[i]!=''):
				L = L + ','+convert(input_var_list[i]);
		L = L + "};"
		array_defs.append(L);


k = 0;
for var_name in var_names['Inputs']:
	L = var_name[1] + ' '+ var_name[0] + ' = ' +  var_name[0] +'_list[inputs[' +  str(k) + ']];';
	var_defs.append(L);
	k = k + 1;

for component_name in component_names:
	L = 'int '+component_name + '_index = inputs[' + str(k) + '];';
	var_defs.append(L);
	k = k + 1;

var_defs.append("");

for component_name in component_names:
	component_vars = var_names[component_name];
	for component_var in component_vars:
		L = component_var[1] + ' ' + component_var[0] +'= '+component_var[0] +'_list['+component_name+'_index'+'];'
		var_defs.append(L);





results = [];
outputs_sheet = data.sheet_by_name('Outputs');
outputs_vars = [[outputs_sheet.cell_value(r, c) for r in range(0,outputs_sheet.nrows)] for c in range(outputs_sheet.ncols)]
k = 1;

for i in range(1,len(outputs_vars[0])):
	L = 'result['+str(k)+'] = ' + outputs_vars[0][i] +';'
	results.append(L);
	k = k + 1;

k2 = 0;
for i in range(0,len(var_names['Inputs'])):
	L = 'result['+str(k)+'] = inputs['+str(k2)+'];';
	results.append(L);
	k = k + 1;
	k2 = k2 + 1;

for component_name in component_names:
	L = 'result['+str(k)+'] = inputs['+str(k2)+'];';
	results.append(L);
	k = k + 1;
	k2 = k2 + 1;




randomization = [];


k = 0;
for var_name in var_names['Inputs']:
	L = 'inputs['+ str(k) + '] = random_index(' + var_name[0]+'_list);'
	randomization.append(L);
	k = k + 1;


for component_name in component_names:
	component_vars = var_names[component_name];
	L = 'inputs['+ str(k) + '] = random_index('+component_vars[0][0] +'_list);'
	randomization.append(L);
	k = k + 1;


randomization.append("");



outputs = [];


k = 0;
for var_name in var_names['Inputs']:
	L = 'results[i]['+ str(len(outputs_vars[0])+k) + '] = '+var_name[0]+'_list[inputs['+str(k)+']];';
	outputs.append(L);
	k = k + 1;


for component_name in component_names:
	component_vars = var_names[component_name];
	L = 'results[i]['+ str(len(outputs_vars[0])+k) + '] = '+component_vars[0][0] +'_list[inputs['+str(k)+']];';
	outputs.append(L);
	k = k + 1;

outputs.append("");



eqs = load_from_txt_file('Eqs.cpp');

consts = ["const int num_inputs = "+str(len(var_names['Inputs'])+len(component_names))+";",
"const int num_outputs = "+str(len(outputs_vars[0]))+";","",
"const int x_index = "+str(X)+";",
"const int y_index = "+str(Y)+";",
"const int z_index = "+str(Z)+";",

];

lines = load_from_txt_file('model_template.h');
lines[lines.index('//data'):lines.index('//data')]=consts+array_defs;
lines[lines.index('//inputs'):lines.index('//inputs')]=var_defs;
lines[lines.index('//equations'):lines.index('//equations')]=eqs;
lines[lines.index('//results'):lines.index('//results')]=results;
lines[lines.index('//randoms'):lines.index('//randoms')]=randomization;
lines[lines.index('//outputs'):lines.index('//outputs')]=outputs;

save_to_txt_file(lines,'model.h');


