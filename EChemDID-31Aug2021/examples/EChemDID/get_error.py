import math
import numpy as np
#from matplotlib.pylab import plt

# Inputs
filein = "traj.dump" # dump file
V0 = 8 # Voltage
lcontact = 3 # Size of the contact region
L = 120.-lcontact # Length of the electrode - contact (3*40 - 4)
K = 6. # Diffusion coefficient
nele = 10 # number of elctronic step
na = 10 # number of analytical curves

# Function
def u(x,t):
    s = 0.0
    for n in range(1,100):
        s += -(V0/np.pi)*(1./n)*math.exp(-t*K*(n*np.pi/L)**2.)*math.sin(n*np.pi*x/L)
    return V0*(1-x/L)/2.0+s

def read_dump(filein):
    print("Read dump file..")
    f = open(filein,"r")
    r = f.readlines()
    f.close()
    nat = int(r[3])
    nfr = int(len(r)/(nat+9))
    dfreq = int(r[nat+10])-int(r[1])
    print("Nat:",nat)
    print("Nframes:",nfr)
    print("Dump freq:",dfreq)
    start = 9
    x,y = [],[]
    for i in range(nfr):
        x_,y_ = [],[]
        for j in range(nat):
            s = r[start+j].split()
            x_.append(float(s[4]))
            y_.append(float(s[6]))
        start += nat+9
        x.append(x_)
        y.append(y_)
    return x,y,dfreq,nfr,nat 

# Get dump file electrochemical potential
x,y,dfreq,nfr,nat = read_dump(filein)

# Plot data
#plt.scatter(x[0],y[0],color='r',marker='+',s=20,linewidth=0.6,label="Electrochemical potential data")
#for i in range(1,na): plt.scatter(x[i],y[i],color='r',marker='+',s=20,linewidth=0.6)

# Compute error
err_tot = 0.0
for i in range(1,na): 
    ya = []
    nj = int(len(x[i])/2)
    err = 0.0
    for j in range(nj): 
        ya.append(u(x[i][j],(i)*nele))
        err += (y[i][j]-u(x[i][j],(i)*nele))**2
#    plt.scatter(np.array(x[i][:nj])+lcontact,ya,color='b',s=20,linewidth=0.6,facecolors='none')
    print("Step:",i,"Error:",err/nj)
    err_tot += err/nj
print("Total Error:",err_tot/na)
print("If Total Error < 5% your installation of EChemDID is most likely correct..")
#plt.scatter(np.array(x[i][:nj])+lcontact,ya,color='b',label="Analytical solution",s=20,linewidth=0.6,facecolors='none')

# Analytical function
#ya_list = []
#for i in range(nfr):
#    ya = []
#    for j in range(len(z)):
#        ya.append(u(z[j],time[i]))
#    ya_list.append(ya)
#    print("time = ",time[i])

#plt.plot(z+lcontact,ya_list[0],color='b',label="Analytical solution")
#for i in range(1,na): plt.plot(z+lcontact,ya_list[i],'b')

#plt.tick_params(axis='both', which='major', labelsize=22)
#plt.xlim(0,L)
#plt.ylim(-0.1,(V0/2+0.1))
#plt.xlabel("x [A]",fontsize=22)
#plt.ylabel("$\Phi$ [eV]",fontsize=22)
#plt.legend()
#plt.tight_layout()
#plt.show()
