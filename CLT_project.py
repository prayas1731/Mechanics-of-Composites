

# Import necessary libraries
import numpy as np
import sympy as sym
from tkinter import *
import math
import matplotlib.pyplot as plt
from IPython.display import display
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)




def Ply_properties(E1f, E2f, G12f, nu12f, Em, num, Vf):
	# Calculate ply properties
	# E1f,E2f,G12f,nu12f,Em,num,Vf = float
    E1_ply = Vf*E1f + (1-Vf)*Em
    E2_ply = E2f*Em/(E2f*(1-Vf) + Em*Vf)
    nu12_ply = nu12f*Vf + num*(1-Vf)
    Gm = Em/(2*(1+num))
    G12_ply = Gm*G12f/(G12f*(1-Vf) + Vf*Gm)
    
    return E1_ply, E2_ply, nu12_ply, G12_ply

def get_1stinput():
	# get the input from the GUI
    global E1, E2,nu12,G12,Em,num,Vf
    E1 =  float(E1_in.get())
    E2 = float(E2_in.get())
    nu12 = float(nu12_in.get())
    G12 = float(G12_in.get())
    Em = float(Em_in.get())
    num = float(num_in.get())
    Vf = float(Vf_in.get()) 


# First frame of GUI
root = Tk() 

titlelbl = Label(root, text="Mechanics Of Composite Term Project",anchor="center",fg="white",bg="red")
titlelbl.grid(column=0,row=0,columnspan=3,sticky='WE')

subtitlelbl = Label(root, text="By Prayas Sambhare", anchor="center",fg="white",bg="black")
subtitlelbl.grid(column=0,row=1,columnspan=3,sticky='WE')

root.geometry("550x600")  
root.title("Composites GUI") 

E1lbl = Label(root, text = "E1 of fibre (GPa) : ",fg="white",bg="SlateGray4")  

E2lbl = Label(root, text = "E2 of fibre (GPa) : ",fg="white",bg="SlateGray4") 

nu12lbl = Label(root, text = "nu12 of fibre : ",fg="white",bg="SlateGray4")

G12lbl = Label(root, text = "G12 of fibre (Gpa) : ",fg="white",bg="SlateGray4") 

Emlbl = Label(root, text = "E of matrix (Gpa) : ",fg="white",bg="SlateGray4") 

numlbl = Label(root, text = "nu of matrix : ",fg="white",bg="SlateGray4") 

Vflbl = Label(root, text = "Vf in decimals: ",fg="white",bg="SlateGray4")



 
E1lbl.grid(column=0,row=2,columnspan=1,sticky='EW',padx = 20, pady = 10)
E2lbl.grid(column=0,row=3,columnspan=1,sticky='EW',padx = 20, pady = 10)
nu12lbl.grid(column=0,row=4,columnspan=1,sticky='EW',padx = 20, pady = 10)
G12lbl.grid(column=0,row=5,columnspan=1,sticky='EW',padx = 20, pady = 10)
Emlbl.grid(column=0,row=6,columnspan=1,sticky='EW',padx = 20, pady = 10)
numlbl.grid(column=0,row=7,columnspan=1,sticky='EW',padx = 20, pady = 10)
Vflbl.grid(column=0,row=8,columnspan=1,sticky='EW',padx = 20, pady = 10)


E1_in = Entry(root)  
E2_in = Entry(root)  
nu12_in = Entry(root)
G12_in = Entry(root)
Em_in = Entry(root)
num_in = Entry(root)
Vf_in = Entry(root)

E1oen = Entry(root)
E2oen = Entry(root)
nu12oen = Entry(root)
G12oen = Entry(root)

 
E1_in.grid(row = 2, column = 1, sticky='Ew',padx = 20, pady = 10)  
E2_in.grid(row = 3, column = 1, sticky='Ew',padx = 20, pady = 10)  
nu12_in.grid(row = 4, column = 1, sticky='Ew',padx = 20, pady = 10)
G12_in.grid(row = 5, column = 1, sticky='Ew',padx = 20, pady = 10)
Em_in.grid(row = 6, column = 1, sticky='Ew',padx = 20, pady = 10)
num_in.grid(row = 7, column = 1, sticky='Ew',padx = 20, pady = 10)
Vf_in.grid(row = 8, column = 1, sticky='Ew',padx = 20, pady = 10)




Button(root, text="S A V E  D A T A", command = get_1stinput).grid(row = 17, column = 0, pady = 10,padx=100)

Button(root, text="N E X T", command= root.destroy).grid(row = 19, column = 0, pady = 10,padx=100)
 
root.mainloop()
# First frame of GUI ends

E1_ply, E2_ply, nu12_ply, G12_ply = Ply_properties(E1, E2, G12, nu12, Em, num, Vf)


def C_basicmatrix(E1_ply, E2_ply, nu12_ply, G12_ply):
	# Calculate compliance matrix for th eeffective properties
	# E1_ply, E2_ply, nu12_ply, G12_ply = float
    C = np.matrix([[1/E1_ply, -nu12_ply/E1_ply, 0],
                   [-nu12_ply/E1_ply, 1/E2_ply, 0],
                   [0, 0, 1/G12_ply]])
    return C


def Q_basicmatrix(C):
	# Calculate stiffness matrix
	# C = matrix
    return np.linalg.inv(C)

def Q_rotation(Q, angle):
	# Rotate stiffness matrix 
	# Q = matrix, angle = angle at which we need to rotate(in degrees)
    angle = angle * np.pi/180  # convert to radians
    m = np.cos(angle)
    n = np.sin(angle)
    T1 = np.matrix([[m**2, n**2, 2*m*n],
                    [n**2, m**2, -2*m*n],
                    [-m*n, m*n, m**2-n**2]])
    T2 = np.matrix([[m**2, n**2, m*n],
                    [n**2, m**2, -m*n],
                    [-2*m*n, 2*m*n, m**2-n**2]])
    Q_rotation = np.linalg.inv(T1) * Q * T2
    return Q_rotation

def C_rotation(C, angle):
	# Rotate stiffness matrix 
	# C = matrix, angle = angle at which we need to rotate(in degrees)
    angle = angle * np.pi/180  # convert to radians
    m = np.cos(angle)
    n = np.sin(angle)
    T1 = np.matrix([[m**2, n**2, 2*m*n],
                    [n**2, m**2, -2*m*n],
                    [-m*n, m*n, m**2-n**2]])
    T2 = np.matrix([[m**2, n**2, m*n],
                    [n**2, m**2, -m*n],
                    [-2*m*n, 2*m*n, m**2-n**2]])
    C_rotation = np.linalg.inv(T2) * C * T1
    return C_rotation

# Change effective properties from Gpa to Pa 
E1_plyd = E1_ply*10**9
E2_plyd = E2_ply*10**9
G12_plyd = G12_ply*10**9


C = C_basicmatrix(E1_plyd, E2_plyd, nu12_ply, G12_plyd)
Q = Q_basicmatrix(C)


def E1vtheta():
      
    # Code to plot E1 v/s theta in gui
    newWindow = Toplevel(root1)
    newWindow.wm_title("Graph of E1 vs theta")

    fig = Figure(figsize=(5, 4), dpi=100)
    x = np.linspace(0, 90, 100)
    angle = np.linspace(0, 90, 100)
    angle = angle*np.pi/180
    m = np.cos(angle)
    n = np.sin(angle)
    E1_plot = 1/(m**4/E1_ply + m**2*n**2*(1/G12_ply - 2*nu12_ply/E1_ply) + n**4/E2_ply)
    fig.add_subplot(111).plot(x,E1_plot)
    fig.suptitle('E1 v/s theta')
    fig.text(0.5, 0.04, 'theta (in degrees)', ha = 'center')
    fig.text(0.04, 0.5, 'E1 (GPa)', ha = 'center', rotation = 'vertical')

    canvas = FigureCanvasTkAgg(fig, master=newWindow)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, newWindow)
    toolbar.update()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    def on_key_press(event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)


    canvas.mpl_connect("key_press_event", on_key_press)
    def _quit():
        newWindow.quit()     # stops mainloop
        newWindow.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    button = Button(master=newWindow, text="Quit", command=_quit)
    button.pack(side=BOTTOM)

    mainloop()

def E2vtheta():

      
    # Code to plot E2 v/s theta in gui
    newWindow = Toplevel(root1)
    newWindow.wm_title("Graph of E2 vs theta")

    fig = Figure(figsize=(5, 4), dpi=100)
    x = np.linspace(0, 90, 100)
    angle = np.linspace(0, 90, 100)
    angle = angle*np.pi/180
    m = np.cos(angle)
    n = np.sin(angle)
    E2_plot = 1/(n**4/E1_ply + m**2*n**2*(1/G12_ply - 2*nu12_ply/E1_ply) + m**4/E2_ply)
    fig.add_subplot(111).plot(x,E2_plot)
    fig.suptitle('E2 vs theta')
    fig.text(0.5, 0.04, 'theta (in degrees)', ha = 'center')
    fig.text(0.04, 0.5, 'E2 (GPa)', ha = 'center', rotation = 'vertical')

    canvas = FigureCanvasTkAgg(fig, master=newWindow)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, newWindow)
    toolbar.update()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    def on_key_press(event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)


    canvas.mpl_connect("key_press_event", on_key_press)
    def _quit():
        newWindow.quit()     # stops mainloop
        newWindow.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    button = Button(master=newWindow, text="Quit", command=_quit)
    button.pack(side=BOTTOM)

    mainloop()


def G12vtheta():
      
    # # Code to plot G12 v/s theta in gui
    newWindow = Toplevel(root1)
    newWindow.wm_title("Graph of G_12 vs theta")

    fig = Figure(figsize=(5, 4), dpi=100)
    x = np.linspace(0, 90, 100)
    angle = np.linspace(0, 90, 100)
    angle = angle*np.pi/180
    m = np.cos(angle)
    n = np.sin(angle)
    G_xy_plot = 1/(4*m**2*n**2*((1/E1_ply)+(1/E2_ply)+(2*nu12_ply/E1_ply)) + (m**2-n**2)**2*(1/G12_ply))
    fig.add_subplot(111).plot(x,G_xy_plot)
    fig.suptitle('G_12 vs theta')
    fig.text(0.5, 0.04, 'theta (in degrees)', ha = 'center')
    fig.text(0.04, 0.5, 'G_12 (GPa)', ha = 'center', rotation = 'vertical')


    canvas = FigureCanvasTkAgg(fig, master=newWindow)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, newWindow)
    toolbar.update()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    def on_key_press(event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)


    canvas.mpl_connect("key_press_event", on_key_press)
    def _quit():
        newWindow.quit()     # stops mainloop
        newWindow.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    button = Button(master=newWindow, text="Quit", command=_quit)
    button.pack(side=BOTTOM)

    mainloop()
    
def Q16vtheta():
      
    # Code to plot Q16 v/s theta in gui
    newWindow = Toplevel(root1)
    newWindow.wm_title("Graph of Q_16 vs theta")

    fig = Figure(figsize=(5, 4), dpi=100)
    x = np.linspace(0, 90, 100)
    angle = np.linspace(0, 90, 100)
    angle = angle * np.pi/180  # convert to radians
    m = np.cos(angle)
    n = np.sin(angle)
    Q16 = m**3*n*(Q[0,0] - Q[0,1]) + m*n**3*(Q[0,1]-Q[1,1]) - 2*m*n*(m**2-n**2)*Q[2,2]
    fig.add_subplot(111).plot(x,Q16)
    fig.suptitle('Q16 vs theta')
    fig.text(0.5, 0.04, 'theta (in degrees)', ha = 'center')
    fig.text(0.04, 0.5, 'Q16 (GPa)', ha = 'center', rotation = 'vertical')


    canvas = FigureCanvasTkAgg(fig, master=newWindow)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, newWindow)
    toolbar.update()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    def on_key_press(event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)


    canvas.mpl_connect("key_press_event", on_key_press)
    def _quit():
        newWindow.quit()     # stops mainloop
        newWindow.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    button = Button(master=newWindow, text="Quit", command=_quit)
    button.pack(side=BOTTOM)

    mainloop()


# Code for second frame of gui
root1 = Tk() 

titlelbl = Label(root1, text="Mechanics Of Composite Term Project",anchor="center",fg="white",bg="red")
titlelbl.grid(column=0,row=0,columnspan=3,sticky='WE')

subtitlelbl = Label(root1, text="By Prayas Sambhare", anchor="center",fg="white",bg="black")
subtitlelbl.grid(column=0,row=1,columnspan=3,sticky='WE')

root1.geometry("740x700")  
root1.title("Composites GUI")

# output labels
E1olbl = Label(root1, text = "E1_ply (Gpa)  ",fg="white",bg="SlateGray4")

E2olbl = Label(root1, text = "E2_ply (Gpa) ",fg="white",bg="SlateGray4")

nu12olbl = Label(root1, text = "nu12_ply  ",fg="white",bg="SlateGray4")

G12olbl = Label(root1, text = "G12_ply(Gpa)  ",fg="white",bg="SlateGray4")

Layuplbl = Label(root1, text = "Enter Layup sequence(comma separated)  ",fg="white",bg="SlateGray4")

Thicknesslbl = Label(root1, text = "Enter thickness of laminate(in mm)  ",fg="white",bg="SlateGray4")

Loadlbl = Label(root1, text = "Enter load [Nx,Ny,Nz,Mx,My,Mz] comma separated (in N/mm)  ",fg="white",bg="SlateGray4")



E1olbl.grid(column=0,row=2,columnspan=1,sticky='EW',padx = 20, pady = 10)
E2olbl.grid(column=0,row=3,columnspan=1,sticky='EW',padx = 20, pady = 10)
nu12olbl.grid(column=0,row=4,columnspan=1,sticky='EW',padx = 20, pady = 10)
G12olbl.grid(column=0,row=5,columnspan=1,sticky='EW',padx = 20, pady = 10)
Layuplbl.grid(column=0,row=24,columnspan=1,sticky='EW',padx = 20, pady = 10)
Thicknesslbl.grid(column=0,row=27,columnspan=1,sticky='EW',padx = 20, pady = 10)
Loadlbl.grid(column=0,row=28,columnspan=1,sticky='EW',padx = 20, pady = 10)

E1oen = Entry(root1, width=50)
E2oen = Entry(root1)
nu12oen = Entry(root1)
G12oen = Entry(root1)
Layupen = Entry(root1)
Thicknessen = Entry(root1)
Loaden = Entry(root1)

E1oen.grid(row = 2, column = 1,sticky='Ew',padx = 20, pady = 10)
E2oen.grid(row = 3, column = 1,sticky='Ew',padx = 20, pady = 10)
nu12oen.grid(row = 4, column = 1,sticky='Ew',padx = 20, pady = 10)
G12oen.grid(row = 5, column = 1,sticky='Ew',padx = 20, pady = 10)
Layupen.grid(row = 24, column = 1,sticky='Ew',padx = 20, pady = 10)
Thicknessen.grid(row = 27, column = 1,sticky='Ew',padx = 20, pady = 10)
Loaden.grid(row=28,column = 1,sticky='Ew',padx = 20, pady = 10)

def display_effprop():
    #Displays the effective properties on the gui
    E1oen.insert(0, round(E1_ply,2))
    E2oen.insert(0, round(E2_ply,2))
    nu12oen.insert(0, round(nu12_ply,2))
    G12oen.insert(0, round(G12_ply,2))

def get_2ndinput():
	# Get the input of layup sequence and load from gui

    global layupseq, t, load
    
    layupstr =  Layupen.get()
    layupseq = [float(x) for x in layupstr.split(',') if x]
    t = float(Thicknessen.get())
    loadstr = Loaden.get()
    load = [float(x) for x in loadstr.split(',') if x]

def sel():
	# get the input of symmetric/unsymmetric from gui
    global selected
    selected = var.get()

var = IntVar()
R1 = Radiobutton(root1, text = "Symmetric matrix", variable = var, value = 1, command = sel)
R2 = Radiobutton(root1, text = "Unsymmetric matrix", variable = var, value = 2, command = sel)

R1.grid(row =25 , column = 1)
R2.grid(row =26 , column = 1)

Button(root1, text="Display ply effective properties", command = display_effprop).grid(row = 17, column = 0, pady = 10,padx=100)
btn1 = Button(root1,text ="E1 vs theta",command = E1vtheta, fg="red",bg="yellow").grid(row = 20, column = 0, pady = 5,padx=100) 
btn2 = Button(root1,text ="E2 vs theta",command = E2vtheta, fg="red",bg="yellow").grid(row = 21, column = 0, pady = 5,padx=100)
btn3 = Button(root1,text ="G12 vs theta",command = G12vtheta, fg="red",bg="yellow").grid(row = 22, column = 0, pady = 5,padx=100)
btn4 = Button(root1,text ="Q16 vs theta",command = Q16vtheta, fg="red",bg="yellow").grid(row = 23, column = 0, pady = 5,padx=100)
btn5 = Button(root1, text="S A V E  D A T A", command = get_2ndinput).grid(row = 29, column = 0, pady = 10,padx=100)

btn6 = Button(root1, text="N E X T", command= root1.destroy).grid(row = 30, column = 0, pady = 10,padx=100)
root1.mainloop()  


if (selected==1):
    k = layupseq
    for el in reversed(k):
        layupseq.append(el)

# Change effective properties from Gpa to Pa         
E1_ply = E1_ply * 10**9
E2_ply = E2_ply * 10**9
G12_ply = G12_ply * 10**9
t = t*10**(-3) #Change thickness of laminate from mm to m
t_lamina = t/len(layupseq)
thickness = [t_lamina]*len(layupseq)


def createABD(Q, angles, thickness):
	# Q = matrix(3,3)
	# angles = list of angles of the layup sequence (in degrees)
	# thickness = list of thickness of each ply
    h = sum(thickness) / 2

    # Create empty matricces for A B en D.
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))

    # Loop over all plies
    for i in range(len(angles)):
        # Calculate the z coordinates of the top and bottom of the ply.
        z_top = np.sum(thickness[:i]) - h
        z_bot = np.sum(thickness[:i+1]) - h

        # Rotate the local stiffenss matrix.
        Q_bar = Q_rotation(Q[i], angles[i])

        # Calculate the contribution to the A, B and D matrix of this layer.
        Ai = Q_bar * (z_bot - z_top)
        Bi = 1/2 * Q_bar * (z_bot**2 - z_top**2)
        Di = 1/3 * Q_bar * (z_bot**3 - z_top**3)

        # Summ this layer to the previous ones.
        A = A + Ai
        B = B + Bi
        D = D + Di

    # Compile the entirety of the ABD matrix.
    ABD = np.matrix([[A[0, 0], A[0, 1], A[0, 2], B[0, 0], B[0, 1], B[0, 2]],
                     [A[1, 0], A[1, 1], A[1, 2], B[1, 0], B[1, 1], B[1, 2]],
                     [A[2, 0], A[2, 1], A[2, 2], B[2, 0], B[2, 1], B[2, 2]],
                     [B[0, 0], B[0, 1], B[0, 2], D[0, 0], D[0, 1], D[0, 2]],
                     [B[1, 0], B[1, 1], B[1, 2], D[1, 0], D[1, 1], D[1, 2]],
                     [B[2, 0], B[2, 1], B[2, 2], D[2, 0], D[2, 1], D[2, 2]]])

    
    return ABD


Qlist = [Q]*len(layupseq)

ABD = createABD(Qlist, layupseq, thickness)

abd = np.linalg.inv(ABD)    # abd = inverse(ABD)

laminate_strain = np.dot(abd,load) # claculate laminate strain
laminate_strain = np.array(laminate_strain)


strain_membrane = [laminate_strain[0,0],laminate_strain[0,1],laminate_strain[0,2]] # mid plane strain
curvature = [laminate_strain[0,3],laminate_strain[0,4],laminate_strain[0,5]] # mid plane curvature

h = np.sum(thickness) / 2
z = []
for i in range(len(layupseq)):
        # Calculate the z coordinates of the top and bottom of the ply.
    z_top = np.sum(thickness[:i]) - h
    z_bot = np.sum(thickness[:i+1]) - h
    z.append((z_top, z_bot))


for i in range(len(z)):
	# Calculate the mean Z cordinate of all plies
    z[i] = (z[i][0]+z[i][1])/2


strain_ply = []
for i in range(len(thickness)):
	# calculate strain in every ply
    l = []
    for j in range(len(strain_membrane)):
        strain_plyz = strain_membrane[j] +  z[i] * curvature[j] # ep = ep0 + z*k
        l.append(strain_plyz)
    strain_ply.append(l)



stress_ply = []
for i in range(len(layupseq)):
	# Calculate stress in every ply [stress] = [Q] * [strain]
    Q_bar = Q_rotation(Q, layupseq[i])
    stress_plyz = np.dot(Q_bar,strain_ply[i])
    stress_plyz = list(stress_plyz)
    stress_ply.append(stress_plyz)

stress_plym = []
for i in range(len(stress_ply)):
    stress = [stress_ply[i][0][0,0], stress_ply[i][0][0,1], stress_ply[i][0][0,2]]
    stress_plym.append(stress)


def stress_rotation(stress, angle):
	# rotate stess
	# stress = list [sigma11, sigma22, sigma12]
	# angle = list of angles in degrees
    angle = angle * np.pi/180  # convert to radians
    m = np.cos(angle)
    n = np.sin(angle)
    T1_inv = np.matrix([[m**2, n**2, 2*m*n],
                        [n**2, m**2, -2*m*n],
                        [-m*n, m*n, m**2-n**2]])
    stress_rot = np.dot(T1_inv,stress)
    return stress_rot


stress_ply_rot = []
for i in range(len(layupseq)):
    stress_ply_rot.append(stress_rotation(stress_plym[i], layupseq[i]))
stress_ply_rotl = []
for i in range(len(stress_ply_rot)):
    stress = [stress_ply_rot[i][0][0,0], stress_ply_rot[i][0][0,1], stress_ply_rot[i][0][0,2]]
    stress_ply_rotl.append(stress)


def display_ABD():
	# Display ABD matrix on gui
    for i in range(6):
        for j in range(6):

            e = Entry()

            e.grid(row=i+3, column=j, sticky=NSEW)

            e.insert(END, round(ABD[i, j],2))

#             cols.append(e)

#     rows.append(cols)
def getinput3():
	# get input from gui 
    global X, Xd, Y, Yd, S12, S23
    X =  float(X_in.get())
    Xd = float(Xd_in.get())
    Y = float(Y_in.get())
    Yd = float(Yd_in.get())
    S12 = float(S12_in.get())
    S23 = float(S23_in.get())

# Code for 3rd frame of gui
root2 = Tk() 

titlelbl = Label(root2, text="Mechanics Of Composite Term Project",anchor="center",fg="white",bg="red")
titlelbl.grid(column=0,row=0,columnspan=6,sticky='WE')

subtitlelbl = Label(root2, text="By Prayas Sambhare", anchor="center",fg="white",bg="black")
subtitlelbl.grid(column=0,row=1,columnspan=6,sticky='WE')




root2.geometry("1050x600")  
root2.title("Composites GUI")

Button(root2, text="Display ABD matrix", command = display_ABD,fg="red",bg="yellow").grid(row = 3, column = 2, pady = 10,padx=100)

    

Xlbl = Label(root2, text = "Enter X(Mpa) : ",fg="white",bg="SlateGray4")  

Xdlbl = Label(root2, text = "Enter X_dash(Mpa) : ",fg="white",bg="SlateGray4") 

Ylbl = Label(root2, text = "Enter Y(Mpa) : ",fg="white",bg="SlateGray4")

Ydlbl = Label(root2, text = "Enter Y_dash(Mpa) : ",fg="white",bg="SlateGray4") 

S12lbl = Label(root2, text = "Enter S12(Mpa) : ",fg="white",bg="SlateGray4") 

S23lbl = Label(root2, text = "Enter S23(Mpa) : ",fg="white",bg="SlateGray4") 


Xlbl.grid(column=0,row=10,columnspan=1,sticky='EW',padx = 20, pady = 10)
Xdlbl.grid(column=0,row=11,columnspan=1,sticky='EW',padx = 20, pady = 10)
Ylbl.grid(column=0,row=12,columnspan=1,sticky='EW',padx = 20, pady = 10)
Ydlbl.grid(column=0,row=13,columnspan=1,sticky='EW',padx = 20, pady = 10)
S12lbl.grid(column=0,row=14,columnspan=1,sticky='EW',padx = 20, pady = 10)
S23lbl.grid(column=0,row=15,columnspan=1,sticky='EW',padx = 20, pady = 10)



X_in = Entry(root2)  
Xd_in = Entry(root2)  
Y_in = Entry(root2)
Yd_in = Entry(root2)
S12_in = Entry(root2)
S23_in = Entry(root2)

X_in.grid(row = 10, column = 1, sticky='Ew',padx = 20, pady = 10)  
Xd_in.grid(row = 11, column = 1, sticky='Ew',padx = 20, pady = 10)  
Y_in.grid(row = 12, column = 1, sticky='Ew',padx = 20, pady = 10)
Yd_in.grid(row = 13, column = 1, sticky='Ew',padx = 20, pady = 10)
S12_in.grid(row = 14, column = 1, sticky='Ew',padx = 20, pady = 10)
S23_in.grid(row = 15, column = 1, sticky='Ew',padx = 20, pady = 10)


Button(root2, text="S A V E  D A T A", command = getinput3).grid(row = 16, column = 2, pady = 10,padx=100)

Button(root2, text="Proceed to Hashin Failure", command = root2.destroy).grid(row = 17, column = 2, pady = 10,padx=100)

root2.mainloop()


# Convert to SI units
X = X*10**6
Xd = Xd*10**6
Y = Y*10**6
Yd = Yd*10**6
S12 = S12*10**6
S23 = S23*10**6


import sympy as sym
R = sym.symbols('R')
R_overall = []
for i in range(len(stress_ply_rotl)):
    if(stress_ply_rotl[i][1]>0):           # Matrix tensile failure mode 
        p = (stress_ply_rotl[i][1]*R/(Y*Y))**2 + (stress_ply_rotl[i][2]*R/S12)**2-1
    else:                                  # Matrix compreessive failure mode
        p = ((Yd/2*S23)**2-1)*(R*stress_ply_rotl[i][1])/Yd + ((R*stress_ply_rotl[i][1])/(4*S23))**2 + (R*stress_ply_rotl[i][2]/S12)**2-1
    R_mat = sym.solve(p, R)
    
    if(stress_ply_rotl[i][0]>0):           # Tensile fibre failure mode
        k = (stress_ply_rotl[i][0]*R/X)**2 + (stress_ply_rotl[i][2]*R/S12)**2 - 1
    else:                                  # Compressive fibre failure mode
        k = (stress_ply_rotl[i][0]*R/Xd)**2 - 1
    R_f = sym.solve(k, R)
    
    # print("R_mat ",R_mat)
    # print("R_f", R_f)
    Rfm = abs(R_f[0])
    Rmm=[]
    for i in R_mat:
        if i.is_real:
            Rmm.append(i)
    Rmf = 0
    for i in Rmm:
        if(i>=0):
            Rmf = i
    if(Rfm<Rmf):
        R_overall.append((Rfm,'fibre will fail first'))
    else:
        R_overall.append((Rmf,'Matrix will fail first'))

R_matrix = []
for i in range(len(layupseq)):
    k = []
    k.append(layupseq[i])
    k.append(R_overall[i])
    R_matrix.append(k)


def display_Rvalues():
	# Display Rvalues on gui
    for i in range(len(layupseq)):
        for j in range(2):

            e = Entry()

            e.grid(row=i+5, column=j+1, sticky=NSEW)

            e.insert(END, R_matrix[i][j])

#             cols.append(e)

#     rows.append(cols)

# Code for 4th frame of gui
root3 = Tk() 

titlelbl = Label(root3, text="Mechanics Of Composite Term Project",anchor="center",fg="white",bg="red")
titlelbl.grid(column=0,row=0,columnspan=6,sticky='WE')

subtitlelbl = Label(root3, text="By Prayas Sambhare", anchor="center",fg="white",bg="black")
subtitlelbl.grid(column=0,row=1,columnspan=6,sticky='WE')


Button(root3, text="Display R values", command = display_Rvalues).grid(row = 3, column = 2, pady = 10,padx=100)

root3.geometry("500x500")  
root3.title("Composites GUI")


root3.mainloop()

