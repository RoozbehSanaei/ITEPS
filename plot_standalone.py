import xlrd
import numpy
from model import model
import matplotlib.pyplot as plt
import matplotlib


book = xlrd.open_workbook('data.xlsx')

outputs_sheet = book.sheet_by_name('Outputs')
output_variables = [[outputs_sheet.cell_value(r, c) for r in range(0,outputs_sheet.nrows)] for c in range(outputs_sheet.ncols)]

number_of_threads = 20;

from multiprocessing.pool import ThreadPool
pool = ThreadPool(processes=number_of_threads+2)


async_result = [0]*number_of_threads;
return_val = [0]*number_of_threads;


for i in range(number_of_threads):
	print("thread started : ", i)
	async_result[i] = pool.apply_async(model, ()) # tuple of args for foo

# do some other stuff in the main process

for i in range(number_of_threads):
	return_val[i] = async_result[i].get()
	print("thread ended : ", i)
 



return_vals = [return_val[i] for i in range(number_of_threads)];



results=numpy.concatenate(return_vals, axis=0);

res1 = [res[1] for res in results];
res2 = [res[2] for res in results];
res3 = [res[3] for res in results];

ind1 = 1
ind2 = 2
ind3 = 3

X_points = []
Y_points = []
Z_points = []
I_points = []

X_points1 = []
Y_points1 = []
Z_points1 = []
I_points1 = []

res_x =  0.02;
res_y =  0.02;
res_z =  0.02;


a_min = numpy.amin(results,axis=0)
a_max = numpy.amax(results,axis=0)
print(a_min.shape)
xmin = a_min[ind1];
xmax = a_max[ind1];
ymin = a_min[ind2];
ymax = a_max[ind2];


n_y = int((ymax-ymin)/res_y)
n_x = int((xmax-xmin)/res_x)



Bins = [[[-10000,-10000,-10000] for r in range(n_y + 2)] for c in range(n_x + 2)]
Bins1 = [[-10000,-10000,-10000] for r in range(n_x + 2)]


for i in range(len(results)):
	X_value = results[i][ind1]
	Y_value = results[i][ind2]
	Z_value = results[i][ind3]
	I_value = i


	p_x = int((X_value-xmin)/res_x);
	p_y = int((Y_value-ymin)/res_y);
	
	if Z_value>Bins[p_x][p_y][2]:
		Bins[p_x][p_y] = (X_value,Y_value,Z_value,I_value)
	
	if Y_value>Bins1[p_x][1]:
		Bins1[p_x] = (X_value,Y_value,Z_value,I_value)

for Bt in Bins :
	for B in Bt :
		if (B[2]>-10000) :
			X_points.append(round(B[0]/res_x)*res_x)
			Y_points.append(round(B[1]/res_y)*res_y)
			Z_points.append(round(B[2]/res_z)*res_z)
			I_points.append(B[3])


for B in Bins1 :
	if (B[1]>-10000) :
		X_points1.append(round(B[0]/res_x)*res_x)
		Y_points1.append(round(B[1]/res_y)*res_y)
		Z_points1.append(round(B[2]/res_z)*res_z)
		I_points1.append(B[3])


		
fig, ax = plt.subplots()

p = ax.scatter(X_points, Y_points,c =Z_points, s=10, edgecolor='',cmap=matplotlib.cm.jet,picker = 1)
#ax.scatter(X_points1, Y_points1, color='black',s=20, edgecolor='',cmap=matplotlib.cm.jet,picker = 1)

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
#ax.set_xlabel(output_variables[0][ind1])
#ax.set_ylabel(output_variables[0][ind2])



cbar = fig.colorbar(p)
#cbar.set_label(output_variables[0][ind3], rotation=270)



plt.show()
