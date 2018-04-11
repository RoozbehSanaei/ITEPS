# This code is based on Anaconda 4.1.1
import os
from shutil import copyfile
import tkFileDialog
import tkFont
import Tkinter
from Tkinter import *
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
import subprocess
import pickle
import xlrd
import ttk
import matplotlib
import matplotlib.pyplot as plt
import numpy
from multiprocessing.pool import ThreadPool
import math

global itemlist,E1min_text,E1max_text,E2min_text,E2max_text,output_variables,numpy_input,numpy_data,number_of_threads,results




def compute_model():
	global number_of_threads,results
	os.system("build27.bat "+str(G1.current()+1)+" "+str(G2.current()+1)+" "+str(G3.current()+1)) # run command
	from model import model
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
		print("thread ended : ", i)
		return_val[i] = async_result[i].get() 

	return_vals = [return_val[i] for i in range(number_of_threads)];
	results=numpy.concatenate(return_vals, axis=0);

	G1["state"] = "readonly"
	G2["state"] = "readonly"
	G3["state"] = "readonly"

	B4["state"] = "normal"
	B5["state"] = "disabled"

	label["text"] = "process is completed!"
	B4["state"] = "normal"
	B5["state"] = "disabled"

	print("computation is done!");




	

	

	
def save_results_file():
	global results,output_variables
	ftypes = [('pickle Files', '*.pickle')]
	fl = tkFileDialog.asksaveasfile(defaultextension='.pickle',mode='w', filetypes = ftypes)
	import pickle
	with open(fl.name, 'wb') as fp:
		pickle.dump([output_variables,results], fp)    
	E_File.delete(0,END)
	E_File.insert(0,fl.name)
	B5["state"] = "normal"
	

   
def load_results_file():
	global results,output_variables
	ftypes = [('pickle Files', '*.pickle')]
	fl = tkFileDialog.askopenfilename(defaultextension='.pickle', filetypes = ftypes)
	with open (fl, 'rb') as fp:
		[output_variables,results] = pickle.load(fp)
	E_File.delete(0,END)
	E_File.insert(0,fl)
	E_File.state = "disabled"
	B4["state"] = "normal"

	C1["values"] = output_variables[0]
	C2["values"] = output_variables[0]
	C3["values"] = output_variables[0]
	C1.current(0)
	C2.current(0)
	C3.current(0)
	C1["state"] = "readonly"
	C2["state"] = "readonly"
	C3["state"] = "readonly"
       


def set_model_file():
	if os.path.isfile(dir_path+'/Eqs.cpp'):
		os.remove(dir_path+'/Eqs.cpp')  
	tf = tkFileDialog.Open()
	fl = tf.show()
	E1.delete(0,END)
	E1.insert(0,fl)
	src_model=E1.get()      
	copyfile(src_model, dir_path+'/Eqs.cpp')

def set_data_file():
	if os.path.isfile(dir_path+'\data.xlsx'):
		os.remove(dir_path+'\data.xlsx')    
	global output_variables,numpy_input,numpy_data,number_of_threads
	tf = tkFileDialog.Open()
	fl = tf.show()
	E2.delete(0,END)
	E2.insert(0,fl)
	src_excel_file=E2.get()
	
	copyfile(src_excel_file, dir_path+'/data.xlsx')
	book = xlrd.open_workbook('data.xlsx')
	outputs_sheet = book.sheet_by_name('Outputs')
	output_variables = [[outputs_sheet.cell_value(r, c) for r in range(1,outputs_sheet.nrows)] for c in range(outputs_sheet.ncols)]
	C1["values"] = output_variables[0]
	C2["values"] = output_variables[0]
	C3["values"] = output_variables[0]
	C1.current(0)
	C2.current(0)
	C3.current(0)
	C1["state"] = "readonly"
	C2["state"] = "readonly"
	C3["state"] = "readonly"
	G1["values"] = output_variables[0]
	G2["values"] = output_variables[0]
	G3["values"] = output_variables[0]
	G1.current(0)
	G2.current(0)
	G3.current(0)

   
	
def plot_figure():
	global ind,ax,X_points,Y_points,Z_points,itemlist,E1min_text,E1max_text,E2min_text,E2max_text,output_variables,xmin,xmax,ymin,ymax,X_points1,Y_points1,Z_points1
	X_points = []
	Y_points = []
	Z_points = []
	
	X_points1 = []
	Y_points1 = []
	Z_points1 = []

	res_x =  0.05;
	res_y =  0.05;
	a_min = numpy.amin(results,axis=0)
	a_max = numpy.amax(results,axis=0)
	print a_min.shape
	xmin = a_min[C1.current()+1];
	xmax = a_max[C1.current()+1];
	ymin = a_min[C2.current()+1];
	ymax = a_max[C2.current()+1];
	
	
	n_y = int((ymax-ymin)/res_y)
	n_x = int((xmax-xmin)/res_x)
	


	Bins = [[[-10000,-10000,-10000] for r in range(n_y + 2)] for c in range(n_x + 2)]
	Bins1 = [[-10000,-10000,-10000] for r in range(n_x + 2)]

	
	for i in range(len(results)):
		X_value = results[i][C1.current()+1]
		Y_value = results[i][C2.current()+1]
		Z_value = results[i][C3.current()+1]

		X_value1 = results[i][C1.current()+1]
		Y_value1 = results[i][C2.current()+1]
		Z_value1 = results[i][C3.current()+1]

		p_x = int((X_value-xmin)/res_x);
		p_y = int((Y_value-ymin)/res_y);
		
		if Z_value>Bins[p_x][p_y][2]:
			Bins[p_x][p_y] = (X_value,Y_value,Z_value)
		
		if Y_value>Bins1[p_x][1]:
			Bins1[p_x] = (X_value,Y_value,Z_value)
	
	for Bt in Bins :
		for B in Bt :
			if (B[2]>-10000) :
				X_points.append(B[0])
				Y_points.append(B[1])
				Z_points.append(B[2])

	for B in Bins1 :
		if (B[1]>-10000) :
			X_points1.append(B[0])
			Y_points1.append(B[1])
			Z_points1.append(B[2])


			
	fig, ax = plt.subplots()

	p = ax.scatter(X_points, Y_points,c =Z_points, s=20, edgecolor='',cmap=matplotlib.cm.jet,picker = 1)
	ax.scatter(X_points1, Y_points1, color='black',s=20, edgecolor='',cmap=matplotlib.cm.jet,picker = 1)
    
	ax.set_xlim(xmin,  xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_xlabel(C1.get())
	ax.set_ylabel(C2.get())



	cbar = fig.colorbar(p)
	cbar.set_label(C3.get(), rotation=270)

	fig.canvas.mpl_connect('pick_event', onpick)
	plt.show()



def onpick(event):
	global ind,ax,X_points,Y_points,Z_points,X_points1,Y_points1,Z_points1
	f =  float(Stel.get())
	ind = event.ind
	print('onpick3 scatter:', ind)
	X0 = X_points[ind[0]]
	Y0 = Y_points[ind[0]]
	Z0 = Z_points[ind[0]]
	Xs = []
	Ys = []
	Xs.append(X0)
	Ys.append(Y0)
	plt.cla()	
	ax.scatter(X_points, Y_points,c =Z_points, s=5, edgecolor='',cmap=matplotlib.cm.jet,picker = 1)
	ax.plot(Xs, Ys, marker='o', markersize=3, color="black",linewidth=2.0)
	for i in range(100):
		neighbors = [i for i in range(len(X_points)) if (math.sqrt((X_points[i]-X0)**2+(Y_points[i]-Y0)**2) < f)]
		neighbors_Z = [Z_points[neighbor] for neighbor in neighbors]
		highest_neighbor = neighbors[neighbors_Z.index(max(neighbors_Z))]

		eligible_points = [i for i in range(len(X_points)) if (Z_points[i]-Z0)>f]
		if (len(eligible_points)>0):
			eligible_points_d = [math.sqrt((X_points[point]-X0)**2+(Y_points[point]-Y0)**2) for point in eligible_points]
			closest_point = eligible_points[eligible_points_d.index(min(eligible_points_d))]
			X0 = X_points[closest_point]
			Y0 = Y_points[closest_point]
			Z0 = Z_points[closest_point]
			Xs.append(X0)
			Ys.append(Y0)
	
	plt.cla()	
	ax.scatter(X_points, Y_points,c =Z_points, s=5, edgecolor='',cmap=matplotlib.cm.jet,picker = 1)
	ax.scatter(X_points1, Y_points1, color='black',s=5, edgecolor='',cmap=matplotlib.cm.jet,picker = 1)
	ax.plot(Xs, Ys, marker='o', markersize=3, color="black")

	ax.set_xlim(xmin,  xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_xlabel(C1.get())
	ax.set_ylabel(C2.get())

	plt.show()

def my_round(x):
	return str(round(x*1000)/1000.00)
	

def combobox_selection(event):
	global E1min_text,E1max_text,E2min_text,E2max_text,E3min_text,E3max_text
	print C1.current()
	print C2.current()

	
root = Tk()
root.title('Trade-off analysis GUI')
root.geometry('900x450')
root.resizable(width=False, height=False)

selected_font = tkFont.Font(family='DejaVu Sans', size=12)

label = Label(text=" ",  font=selected_font)
label.place(bordermode=OUTSIDE, x=500,y=100,height=30, width=200)

E_File = Entry(text="File : Not Loaded",  font=selected_font)
E_File.place(bordermode=OUTSIDE, x=160,y=260,height=30, width=600)


B1 = Button(text='model file :', font=selected_font,command = set_model_file)
E1 = Entry(font=selected_font)
B1.place(bordermode=OUTSIDE, x=10,y=10,height=30, width=150)
E1.place(bordermode=OUTSIDE, x=170,y=10,height=30, width=600)

B2 = Button(text='data file : ', font=selected_font,command = set_data_file)
E2 = Entry(font=selected_font)
B2.place(bordermode=OUTSIDE, x=10,y=60,height=30, width=150)
E2.place(bordermode=OUTSIDE, x=170,y=60,height=30, width=600)



B3 = Button(text='compute', font=selected_font,command = compute_model)
B3.place(bordermode=OUTSIDE, x=350,y=200,height=30, width=100)



B4 = Button(text='save...', font=selected_font,command = save_results_file,state = "disabled")
B4.place(bordermode=OUTSIDE, x=20,y=240,height=30, width=100)

B5 = Button(text='load...', font=selected_font,command = load_results_file)
B5.place(bordermode=OUTSIDE, x=20,y=270,height=30, width=100)





C1 = ttk.Combobox(font=selected_font,state = "disabled")
C1.place(bordermode=OUTSIDE, x=20,y=350,height=30, width=200)
C1.bind("<<ComboboxSelected>>", combobox_selection)
C2 = ttk.Combobox(font=selected_font,state = "disabled")
C3 = ttk.Combobox(font=selected_font,state = "disabled")
C2.place(bordermode=OUTSIDE, x=300,y=350,height=30, width=200)
C2.bind("<<ComboboxSelected>>", combobox_selection)
C3.place(bordermode=OUTSIDE, x=570,y=350,height=30, width=200)
C3.bind("<<ComboboxSelected>>", combobox_selection)



G1 = ttk.Combobox(font=selected_font,state = "enabled")
G1.place(bordermode=OUTSIDE, x=20,y=130,height=30, width=200)
G1.bind("<<ComboboxSelected>>", combobox_selection)
G2 = ttk.Combobox(font=selected_font,state = "enabled")
G3 = ttk.Combobox(font=selected_font,state = "enabled")
G2.place(bordermode=OUTSIDE, x=300,y=130,height=30, width=200)
G2.bind("<<ComboboxSelected>>", combobox_selection)
G3.place(bordermode=OUTSIDE, x=570,y=130,height=30, width=200)
G3.bind("<<ComboboxSelected>>", combobox_selection)





B6 = Button(text='plot', font=selected_font,command = plot_figure)
B6.place(bordermode=OUTSIDE, x=350,y=400,height=30, width=100)


S_tel = StringVar()
S_tel.set("0.1")
Stel = Spinbox(textvariable = S_tel,from_=0, to=2,increment=0.01)
Stel.place(bordermode=OUTSIDE, x=800,y=350,height=30, width=50)


root.mainloop()
#compute_model()

p = subprocess.Popen([sys.executable, 'plots_threaded.py'], 
									stdout=subprocess.PIPE, 
									stderr=subprocess.STDOUT)
