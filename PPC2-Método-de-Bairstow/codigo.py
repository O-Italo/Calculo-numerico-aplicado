"""
================================================================================
                               MÉTODO DE BAIRSTOW
================================================================================
* Programa: PPC2 - MÉTODO DE BAIRSTOW
* Atividade: Programa Para Casa 02 | Tema: Raízes de Polinômios
*
* Repositório: Cálculo Numérico Aplicado
* Curso: Engenharia Mecatrônica - Universidade de Brasília (UnB)
* Autor: Ítalo B Tavares | Prof. Rafael Gabler Gontijo
*
*       https://github.com/O-Italo/Calculo-numerico-aplicado
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# CABEÇALHO DO TERMINAL
# =============================================================================
def imprimir_cabecalho():
    print("=" * 80)
    print("MÉTODO DE BAIRSTOW")
    print("=" * 80)
    print("* Programa: PPC2 - MÉTODO DE BAIRSTOW")
    print("* Atividade: Programa Para Casa 02 | Tema: Raízes de Polinômios")
    print("*")
    print("* Repositório: Cálculo Numérico Aplicado")
    print("* Curso: Engenharia Mecatrônica - Universidade de Brasília (UnB)")
    print("* Autor: Ítalo B Tavares | Prof. Rafael Gabler Gontijo")
    print("*")
    print("*       https://github.com/O-Italo/Calculo-numerico-aplicado")
    print("=" * 80)
    print()

# =============================================================================
# FUNÇÕES DO MÉTODO DE BAIRSTOW
# =============================================================================

def calcula_b(coeficientes, r, s):
    # Realiza a divisão sintética do polinómio pelo divisor quadrático
    grau = len(coeficientes) - 1
    b = np.zeros(grau + 1)
    
    b[0] = coeficientes[0]
    if grau >= 1:
        b[1] = coeficientes[1] + r * b[0]
        
    for i in range(2, grau + 1):
        b[i] = coeficientes[i] + r * b[i - 1] + s * b[i - 2]
        
    return b

def calcula_c(b, r, s):
    # Segunda divisão sintética para encontrar as derivadas parciais
    grau = len(b) - 1
    c = np.zeros(grau)
    
    c[0] = b[0]
    if grau >= 2:
        c[1] = b[1] + r * c[0]
        
    for i in range(2, grau):
        c[i] = b[i] + r * c[i - 1] + s * c[i - 2]
        
    return c

def iteracao_bairstow(coeficientes, r, s):
    # Resolve o sistema linear 2x2 para encontrar os incrementos delta_r e delta_s
    grau = len(coeficientes) - 1
    b = calcula_b(coeficientes, r, s)
    c = calcula_c(b, r, s)

    b0 = b[grau]
    b1 = b[grau - 1]

    c1 = c[grau - 1]
    c2 = c[grau - 2]
    c3 = c[grau - 3]

    determinante = c2**2 - c1*c3

    # Evitar divisão por zero se o sistema for singular
    if abs(determinante) < 1e-15:
        return 0.0, 0.0, b1, b0

    delta_r = (b0 * c3 - b1 * c2) / determinante
    delta_s = (b1 * c1 - b0 * c2) / determinante

    return delta_r, delta_s, b1, b0

def bairstow(coeficientes, r_inicial=0.5, s_inicial=0.5, tolerancia=1e-10, limite_iteracoes=1000):
    # Normaliza o polinómio (a_n = 1)
    polinomio = np.array(coeficientes, dtype=float)
    polinomio = polinomio / polinomio[0]

    raizes_encontradas = []
    historico_iteracoes = []

    # Loop principal de deflação
    while len(polinomio) > 3:
        grau_atual = len(polinomio) - 1
        r = r_inicial
        s = s_inicial
        convergiu = False

        for iteracao in range(limite_iteracoes):
            delta_r, delta_s, b1, b0 = iteracao_bairstow(polinomio, r, s)

            # Se o delta for nulo por erro de divisão por zero, dá um pequeno salto para sair do ponto
            if delta_r == 0.0 and delta_s == 0.0 and abs(b0) > tolerancia:
                r += 0.01
                s += 0.01
                continue

            r += delta_r
            s += delta_s

            # Calcula o erro relativo
            erro_r = abs(delta_r / r) if r != 0 else abs(delta_r)
            erro_s = abs(delta_s / s) if s != 0 else abs(delta_s)

            if erro_r < tolerancia and erro_s < tolerancia:
                convergiu = True
                break

        if not convergiu:
            print(f"Aviso: Não convergiu para r={r:.4f}, s={s:.4f}.")

        historico_iteracoes.append(iteracao + 1)

        # Calcula as duas raízes do fator quadrático x^2 - rx - s = 0
        delta_eq = r**2 + 4 * s
        if delta_eq >= 0:
            raiz1 = (r + np.sqrt(delta_eq)) / 2
            raiz2 = (r - np.sqrt(delta_eq)) / 2
        else:
            real = r / 2
            imag = np.sqrt(-delta_eq) / 2
            raiz1 = complex(real, imag)
            raiz2 = complex(real, -imag)

        raizes_encontradas.extend([raiz1, raiz2])

        # Reduz o grau do polinómio usando os coeficientes b (deflação)
        b = calcula_b(polinomio, r, s)
        polinomio = b[:grau_atual - 1]

    # Trata o que sobrou do polinómio (grau 1 ou 2)
    if len(polinomio) == 3:
        a2, a1, a0 = polinomio
        delta_eq = a1**2 - 4 * a2 * a0
        if delta_eq >= 0:
            raizes_encontradas.append((-a1 + np.sqrt(delta_eq)) / (2 * a2))
            raizes_encontradas.append((-a1 - np.sqrt(delta_eq)) / (2 * a2))
        else:
            raizes_encontradas.append((-a1 + 1j * np.sqrt(-delta_eq)) / (2 * a2))
            raizes_encontradas.append((-a1 - 1j * np.sqrt(-delta_eq)) / (2 * a2))
        historico_iteracoes.append(0)

    elif len(polinomio) == 2:
        raizes_encontradas.append(-polinomio[1] / polinomio[0])
        historico_iteracoes.append(0)

    return np.array(raizes_encontradas, dtype=complex), historico_iteracoes


# =============================================================================
# ROTINAS DE ANÁLISE SOLICITADAS NO GUIÃO
# =============================================================================

def analise_validacao_grau7():
    print("--- 1. Validação com Polinómio de Grau 7 ---")
    
    # Raízes arbitrárias conhecidas para teste
    raizes_reais = [1.0, -2.0, 3.0]
    raizes_complexas = [complex(1, 2), complex(1, -2), complex(-1, 1), complex(-1, -1)]
    raizes_esperadas = np.array(raizes_reais + raizes_complexas)
    
    # Gera os coeficientes a partir das raízes
    coeficientes_teste = np.poly(raizes_esperadas)
    
    raizes_calculadas, iteracoes = bairstow(coeficientes_teste, r_inicial=0.5, s_inicial=-0.5)
    
    print(f"Iterações gastas por raiz isolada: {iteracoes}")
    print("\nComparação de Raízes:")
    
    # Ordena para facilitar a visualização no terminal
    esperadas_ord = np.sort_complex(raizes_esperadas)
    calculadas_ord = np.sort_complex(raizes_calculadas)
    
    for resp, calc in zip(esperadas_ord, calculadas_ord):
        erro = abs(resp - calc)
        print(f"Esperada: {np.round(resp, 4)} | Calculada: {np.round(calc, 4)} | Erro: {erro:.2e}")

    # Gráfico simples
    plt.figure(figsize=(7, 5))
    plt.scatter(raizes_esperadas.real, raizes_esperadas.imag, s=100, label="Esperadas", marker="o")
    plt.scatter(raizes_calculadas.real, raizes_calculadas.imag, s=50, label="Bairstow", marker="x")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.xlabel("Eixo Real")
    plt.ylabel("Eixo Imaginário")
    plt.title("Validação do Método de Bairstow (Grau 7)")
    plt.legend()
    plt.grid(True)
    plt.show()

def analise_sistema_apc2():
    print("\n--- 2. Análise do Polinómio da APC2 ---")
    
    # P(lambda) = 2x^4 + 5x^3 + 12x^2 + 8x + 8
    coeficientes_apc2 = [2.0, 5.0, 12.0, 8.0, 8.0]
    
    raizes_apc2, iteracoes = bairstow(coeficientes_apc2, r_inicial=0.5, s_inicial=0.5)
    
    print("Autovalores encontrados (Lambda):")
    for i, raiz in enumerate(raizes_apc2):
        parte_real = raiz.real
        if parte_real < -1e-6:
            status = "Estável"
        else:
            status = "Instável"
            
        print(f"Lambda {i+1}: {np.round(raiz, 4)} -> {status}")

    plt.figure(figsize=(6, 5))
    for raiz in raizes_apc2:
        plt.plot(raiz.real, raiz.imag, 'ro', markersize=8)
    
    plt.axhline(0, color='black', linewidth=0.8)
    plt.axvline(0, color='black', linewidth=0.8)
    plt.xlabel("Parte Real")
    plt.ylabel("Parte Imaginária")
    plt.title("Autovalores do Sistema 2-GDL")
    plt.grid(True)
    plt.show()

def analise_sensibilidade():
    print("\n--- 3. Sensibilidade aos chutes iniciais (r0, s0) ---")
    coeficientes_apc2 = [2.0, 5.0, 12.0, 8.0, 8.0]
    
    valores_teste_r = [-1.0, 0.0, 1.0, 2.0]
    valores_teste_s = [-1.0, 0.0, 1.0, 2.0]
    
    print("r0\ts0\tIterações Totais")
    print("-" * 35)
    
    for r0 in valores_teste_r:
        for s0 in valores_teste_s:
            _, iteracoes = bairstow(coeficientes_apc2, r_inicial=r0, s_inicial=s0, limite_iteracoes=100)
            total_it = sum(iteracoes)
            print(f"{r0}\t{s0}\t{total_it}")

def gerar_fractal_bairstow():
    print("\n--- 4. Gerando Fractal de Bairstow ---")
    print("Isso pode demorar um pouco...")
    
    coeficientes = [2.0, 5.0, 12.0, 8.0, 8.0]
    resolucao = 200
    max_it = 50
    
    r_vetor = np.linspace(-3, 3, resolucao)
    s_vetor = np.linspace(-3, 3, resolucao)
    matriz_iteracoes = np.zeros((resolucao, resolucao))
    
    for i in range(resolucao):
        for j in range(resolucao):
            r0 = r_vetor[j]
            s0 = s_vetor[i]
            
            polinomio = np.array(coeficientes, dtype=float)
            polinomio = polinomio / polinomio[0]
            
            r, s = r0, s0
            iteracoes_gastas = max_it
            
            for k in range(max_it):
                if abs(r) > 100 or abs(s) > 100:
                    break
                    
                dr, ds, _, _ = iteracao_bairstow(polinomio, r, s)
                
                if dr == 0.0 and ds == 0.0:
                    break
                    
                r += dr
                s += ds
                
                erro_r = abs(dr/r) if r != 0 else abs(dr)
                erro_s = abs(ds/s) if s != 0 else abs(ds)
                
                if erro_r < 1e-6 and erro_s < 1e-6:
                    iteracoes_gastas = k
                    break
                    
            matriz_iteracoes[i, j] = iteracoes_gastas

    plt.figure(figsize=(7, 6))
    plt.imshow(matriz_iteracoes, origin='lower', extent=[-3, 3, -3, 3], cmap='hot')
    plt.colorbar(label="Iterações")
    plt.xlabel("Valor inicial de r0")
    plt.ylabel("Valor inicial de s0")
    plt.title("Fractal de Bairstow (Convergência)")
    plt.show()

# =============================================================================
# EXECUÇÃO DO SCRIPT
# =============================================================================
if __name__ == "__main__":
    imprimir_cabecalho()
    analise_validacao_grau7()
    analise_sistema_apc2()
    analise_sensibilidade()
    gerar_fractal_bairstow()

