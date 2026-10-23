
import numpy as np
import matplotlib.pyplot as plt

def compare_results_cylinder(points,u_num,sig_num,E,nu,b,a,Pressure):
    ne=points.shape[0]
    u_an=np.zeros(ne)
    sig_an=np.zeros((ne,2))
    for i in range(ne):
        P=Pressure
        x=points[i,0]
        y=points[i,1]
        r=np.sqrt(x**2+y**2)
        theta=np.arctan2(y,x)
        sigma_r=(a**2*P)/(b**2-a**2)*(1-b**2/r**2)
        sigma_theta=(a**2*P)/(b**2-a**2)*(1+b**2/r**2)
        u=(1+nu)*P*a**2/(E*(b**2-a**2))*((1-2*nu)*r+b**2/r)
        u_an[i]=u
        sig_an[i,:]=[sigma_r,sigma_theta]

    plt.plot(points[:,0],u_an,'ro-',label = 'u analitical')
    plt.plot(points[:,0],u_num,'bo-',label = 'u numerical')

    plt.legend()
    plt.grid()
    plt.title('Displacements')
    plt.figure()

    plt.plot(points[:,0],sig_an[:,0],'ro-',label = 'sig_r analitical');
    plt.plot(points[:,0],sig_num[:,1],'bo-',label = 'sig_r numerical');

    plt.plot(points[:,0],sig_an[:,1],'rs-',label = 'sig_theta analitical');
    plt.plot(points[:,0],sig_num[:,0],'bs-',label = 'sig_theta numerical');
    plt.legend()
    plt.grid()
    plt.title('Stresses')
    return