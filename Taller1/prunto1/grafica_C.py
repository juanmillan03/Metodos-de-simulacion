import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def Punto_C(data_path):
    datos=pd.read_csv(data_path,sep=" " ,names=['t','s',"i","r"])
    datos_rales=pd.read_csv("Datos_reales.txt",sep="\t" ,names=['t',"i"],skiprows=1)
    fig, ax = plt.subplots( )
    ax.scatter(datos_rales['t'],datos_rales['i'],label="Reales",s=1,c="r")
    ax.plot(datos['t'],datos['i'],label="Simulados", linestyle='-.')
    ax.set_xticks(np.linspace(0, 300,9))
    ax.set_xlabel("Tiempo(Díaz)")
    ax.set_ylabel("Infectados ")
    ax.grid()
    ax.legend()
    plt.savefig("Modelosir_C.pdf")


    
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python graficador.py <data_path>")
        sys.exit(1)
    data_C = sys.argv[1]
    Punto_C(data_C)