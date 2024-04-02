import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def Punto_B(data_path):
    data=pd.read_excel(data_path)
    datos_importantes = data[["dias","Fallecidos","Recuperados","Infectados activos","Infectados factor"]]
    datos_importantes =datos_importantes[datos_importantes["dias"].notna()]
    dia=37
    plt.scatter(datos_importantes["dias"],datos_importantes["Infectados factor"],s=10)
    plt.axvline(x=dia+1, color='r', linestyle='--', linewidth=1)
    plt.scatter(datos_importantes["dias"][dia],datos_importantes["Infectados factor"][dia])
    plt.xlabel("Días")
    plt.ylabel("Factor de infectados")
    plt.grid()
    plt.title("Factor de infectados contra los Días")
    plt.savefig("Modelosir_B.pdf")
    print("----------------------------------\n")

    print("Elegimos el día de inicio de la pandemia")
    print("Día=",datos_importantes["dias"][dia-1],"i=",datos_importantes["Infectados activos"][dia-1])
    print("Valores iniciales")
    i=datos_importantes["Infectados activos"][dia-1]/45000000
    print("t=0 r=0","i=",i,"s=",1-i,"a=",datos_importantes["Infectados factor"][dia-1])
    print("Promedio de los siguientes 120 días")
    i=datos_importantes["Infectados activos"][dia-3:dia+4].mean()/45000000
    print("t=0 r=0","i=",i,"s=",1-i,"a=",datos_importantes["Infectados factor"][dia:dia+120].mean())
    print("----------------------------------\n")

def Punto_A(data_path):
    datos=pd.read_csv(data_path,sep=" " ,names=['t','s',"i","r"])
    fig, ax = plt.subplots( )
    ax.plot(datos['t'],datos['s'],label="Susceptibles", linestyle='-.')
    ax.plot(datos['t'],datos['i'],label="Infectados", linestyle='-.')
    ax.plot(datos['t'],datos['r'],label="Recuperados", linestyle='-.')
    ax.set_xlabel("Tiempo(Díaz)")
    ax.set_ylabel("Susceptibles, infectados y recuperados")
    ax.grid()
    ax.legend()
    plt.title("SIR sin muertes ni nacimentos")
    plt.savefig("modelosir_A.pdf")




    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python graficador.py <data_path>")
        sys.exit(1)
        
    data_B = sys.argv[1]
    Punto_B(data_B)
    data_A = sys.argv[2]
    Punto_A(data_A)
 