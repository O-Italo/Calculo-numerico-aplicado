
# 🔢 Raízes de Polinômios — Método de Bairstow

> **Programa para Casa #2 · Cálculo Numérico Aplicado**


---

## 📌 Sobre o projeto

Este repositório contém a implementação do **método de Bairstow** para a determinação de todas as raízes (reais e complexas) de polinômios de grau arbitrário. O método é aplicado à análise do sistema dinâmico massa-mola-amortecedor com 2 graus de liberdade (GDL) estudado na APC2, cujo polinômio característico é:

$$P(\lambda) = 2\lambda^4 + 5\lambda^3 + 12\lambda^2 + 8\lambda + 8$$

Adicionalmente, o código gera o **fractal de Bairstow**: um mapa de convergência no plano $(r_0, s_0)$ que revela a estrutura geométrica das bacias de atração do método.

**Nenhuma biblioteca de fatoração de polinômios é utilizada** — o método é implementado do zero a partir dos conceitos matemáticos de divisão sintética e do método de Newton-Raphson.

---

## ⚙️ O método de Bairstow

O algoritmo divide iterativamente o polinômio $f_n(x)$ por um fator quadrático real $D(x) = x^2 - rx - s$, ajustando iterativamente $(r, s)$ para anular o resto da divisão.

### 1. Divisão sintética (coeficientes $b_i$)
Para $f_n(x) = a_n x^n + a_{n-1}x^{n-1} + \cdots + a_0$, a divisão por $D(x)$ produz o resto associado a $b_0$ e $b_1$:
$$b_n = a_n, \quad b_{n-1} = a_{n-1} + r\,b_n, \quad b_k = a_k + r\,b_{k+1} + s\,b_{k+2}, \; k = n{-}2, \ldots, 0$$

### 2. Derivadas parciais (coeficientes $c_i$)
Uma segunda divisão sintética fornece as derivadas parciais do resto em relação a $r$ e $s$:
$$c_n = b_n, \quad c_{n-1} = b_{n-1} + r\,c_n, \quad c_k = b_k + r\,c_{k+1} + s\,c_{k+2}, \; k = n{-}2, \ldots, 1$$

### 3. Passo de Newton-Raphson
Resolve-se o sistema $2 \times 2$ para os incrementos $(\Delta r, \Delta s)$ usando a Regra de Cramer:
$$\begin{bmatrix} c_2 & c_3 \\ c_1 & c_2 \end{bmatrix} \begin{bmatrix} \Delta r \\ \Delta s \end{bmatrix} = \begin{bmatrix} -b_1 \\ -b_0 \end{bmatrix}$$

### 4. Extração das raízes
Quando os incrementos são menores que a tolerância estipulada, extraem-se as raízes pela fórmula quadrática:
$$x_{r} = \frac{r \pm \sqrt{r^2 + 4s}}{2}$$
O polinômio é deflacionado (substituído pelo quociente da divisão) e o ciclo recomeça até encontrar todas as raízes.
---

##  Como executar

### Pré-requisitos

- Python 3.8 ou superior
- Bibliotecas: `numpy` e `matplotlib`

```bash
pip install numpy matplotlib
```

### Rodando o programa completo

```bash
PPC2.py
```

O programa executa as quatro análises automaticamente e salva todos os resultados em `resultados/`. A geração do fractal pode levar alguns minutos dependendo da resolução.

---

## Análises realizadas

###  Validação de grau 7

Constrói um polinômio com raízes conhecidas (3 reais + 2 pares complexos conjugados) e verifica se o algoritmo as recupera perfeitamente, plotando um comparativo no plano complexo
---

### Sistema dinâmico 2-GDL (APC2)

Extrai os autovalores do polinômio de 4º grau que rege o sistema 2-GDL. O algoritmo encontra raízes com parte real negativa, confirmando fisicamente a estabilidade do sistema amortecido
---

### Sensibilidade ao chute inicial

Varre uma grade de pontos $(r_0, s_0)$ e imprime o número de iterações gastas, comprovando a forte dependência que o método de Newton-Raphson tem do ponto de partida.
---

### Análise 4 — Fractal de Bairstow
Gera um mapa de calor de $200 \times 200$ pontos no domínio $[-3, 3]^2$. Cada coordenada recebe uma cor baseada na quantidade de iterações necessárias para convergir. inferno (claro = poucos iters, escuro = muitos iters, preto = divergiu)
```


---

## 📐 Tabela de funções principais

| Função | Descrição |
|---|---|
| `calcula_b(coeficientes, r, s)` | Divisão sintética → coeficientes $b_i$ |
| `calcula_c(b, r, s)` | Segunda divisão sintética → coeficientes $c_i$ |
| `iteracao_bairstow(coeficientes, r, s)` | Um passo Newton-Raphson → $(\Delta r, \Delta s)$ |
| `bairstow(coeffs, r0, s0, tol, max_iter)` | Método completo → todas as raízes |
| `analise_validacao_grau7()` | Validação com polinômio grau 7 |
| `analise_sistema_apc2()` | Autovalores do sistema 2-GDL |
| `analise_sensibilidade()` | Grade de chutes iniciais |
| `gerar_fractal_bairstow())` | Mapa fractal no plano $(r_0, s_0)$ |
| `exportar_raizes(...)` | Salva raízes em `.csv` |

---

##  ratamento de Erros e Limitações

- Caso o determinante do sistema linear de Newton-Raphson se aproxime de zero (matriz singular), o algoritmo aplica uma pequena perturbação nos valores atuais (+ 0.01) e reinicia a iteração, evitando quebra do código por erro de divisão por zero.
- A convergência tende a ser rápida (quadrática), mas em áreas escuras do fractal pode falhar em atingir a tolerância no limite de passos.
---

## 📚 Referências

1. Chapra, S. C., Canale, R. P. — *Métodos Numéricos para Engenharia*, McGraw-Hill, 5ª ed. (2008) — Capítulo 9.
2. Bairstow, L. — *Applied Aerodynamics*, Longmans (1920) — Apêndice.
3. Gontijo, R. G. — *Notas de Aula: Cálculo Numérico Aplicado*, UnB (2026).

---
