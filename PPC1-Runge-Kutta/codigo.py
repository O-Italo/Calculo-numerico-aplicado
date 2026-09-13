import numpy as np
import matplotlib.pyplot as plt
import os

# Pasta onde os gráficos serão salvos
os.makedirs("resultados_ppc1", exist_ok=True)
h = 0.01
tmax = 8.0

# Equação para o regime de Stokes
def stokes(t, v, st):
    dvdt = (1 - v) / st
    return dvdt

# Equação de Oseen
def oseen(t, v, st, re):
    dvdt = (1 - v - (3 / 8) * re * v**2) / st
    return dvdt
# Runge-Kutta de 4a ordem
def rk4_stokes(t, v, h, st):
    k1 = stokes(t, v, st)
    k2 = stokes(t + h / 2, v + h * k1 / 2, st)
    k3 = stokes(t + h / 2, v + h * k2 / 2, st)
    k4 = stokes(t + h, v + h * k3, st)
    v_novo = v + h * (k1 + 2*k2 + 2*k3 + k4) / 6
    return v_novo
def rk4_oseen(t, v, h, st, re):
    k1 = oseen(t, v, st, re)
    k2 = oseen(t + h / 2, v + h * k1 / 2, st, re)
    k3 = oseen(t + h / 2, v + h * k2 / 2, st, re)
    k4 = oseen(t + h, v + h * k3, st, re)
    v_novo = v + h * (k1 + 2*k2 + 2*k3 + k4) / 6
    return v_novo

# 1 - Comparação entre RK4 e solução exata para Stokes
st_lista = [0.25, 0.5, 1.0, 2.0]
for st in st_lista:
    t = 0
    v = 0
    tempos = [t]
    velocidades = [v]
    while t < tmax:
        v = rk4_stokes(t, v, h, st)
        t = t + h
        tempos.append(t)
        velocidades.append(v)
        
    # Pontos calculados pelo método numérico
    plt.plot(tempos[::40], velocidades[::40], 'o')
    
    # Solução exata
    v_exato = []
    for tempo in tempos:
        v_exato.append(1 - np.exp(-tempo / st))
    plt.plot(tempos, v_exato, label="St = " + str(st))
plt.xlabel("t*")
plt.ylabel("v*")
plt.title("Regime de Stokes: RK4 e solução exata")
plt.legend()
plt.savefig("resultados_ppc1/tarefa1_stokes_validacao.png")
plt.clf()

# 2 - Erro do método de RK4 para diferentes passos
h_lista = [1.0, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001]
erros = []
for h_teste in h_lista:
    t = 0
    v = 0
    tempos = [t]
    velocidades = [v]
    while t < 5.0:
        v = rk4_stokes(t, v, h_teste, 1.0)
        t = t + h_teste
        tempos.append(t)
        velocidades.append(v)
    erro = 0
    for j in range(len(tempos)):
        v_exato = 1 - np.exp(-tempos[j])
        erro_atual = abs(velocidades[j] - v_exato)
        if erro_atual > erro:
            erro = erro_atual
    erros.append(erro)
plt.plot(h_lista, erros, 'o-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("h")
plt.ylabel("Erro máximo")
plt.title("Convergência do método de RK4")
plt.savefig("resultados_ppc1/tarefa2_convergencia_h.png")
plt.clf()

# 3 - Influência do número de Reynolds
re_lista = [0.0, 0.25, 0.5, 1.0, 2.0]
h = 0.005
st = 1.0
for re in re_lista:
    t = 0
    v = 0
    tempos = []
    velocidades = []
    while t < tmax:
        tempos.append(t)
        velocidades.append(v)
        if re == 0:
            v = rk4_stokes(t, v, h, st)
        else:
            v = rk4_oseen(t, v, h, st, re)
        t = t + h
        
    # Resultado numérico
    plt.plot(tempos[::30], velocidades[::30], 'o')
    # Solução analítica
    velocidade_exata = []
    for tempo in tempos:
        if re == 0:
            v_exato = 1 - np.exp(-tempo / st)
        else:
            raiz = np.sqrt(1 + 1.5 * re)
            vb = (4 / (3 * re)) * (raiz - 1)
            va = -(4 / (3 * re)) * (raiz + 1)
            exponencial = np.exp((va - vb) * (-3 * re / (8 * st)) * tempo)
            v_exato = (va * vb * (1 - exponencial)/ (vb - va * exponencial))
        velocidade_exata.append(v_exato)
    plt.plot(tempos, velocidade_exata, label="Re = " + str(re))
plt.xlabel("t*")
plt.ylabel("v*")
plt.title("Efeito do número de Reynolds")
plt.legend()
plt.savefig("resultados_ppc1/tarefa35_efeito_reynolds.png")
plt.clf()

