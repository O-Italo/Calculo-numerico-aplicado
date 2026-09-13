# 📘 PPC1 — Método de Runge-Kutta

## 📌 Sobre o projeto

Este projeto foi desenvolvido para estudar a aplicação de métodos numéricos na solução de um **problema de valor inicial (PVI)**.

O problema escolhido descreve o movimento de uma esfera sedimentando em um fluido viscoso. A partir da equação física do movimento, foi obtida uma equação diferencial adimensional que pode ser resolvida numericamente.

Para resolver o problema, foi implementado o método de **Runge-Kutta de 4ª ordem (RK4)** em Python.

Além de resolver a equação, o programa realiza algumas análises para verificar o comportamento e a precisão do método numérico.

---

# 🎯 Objetivos

Os principais objetivos deste trabalho são:

- Implementar o método RK4 em Python;
- Resolver numericamente uma equação diferencial;
- Comparar a solução numérica com uma solução analítica conhecida;
- Verificar a convergência do método conforme o passo `h` é alterado;
- Analisar a influência do número de Reynolds na solução;
- Gerar gráficos para visualizar os resultados;
- Salvar automaticamente os resultados obtidos.

---

# 🧰 Tecnologias utilizadas

O programa foi desenvolvido utilizando **Python 3.8+**.

As principais bibliotecas utilizadas foram:

### NumPy

Utilizada para operações numéricas e manipulação de vetores e arrays.

```python
import numpy as np
```

### Matplotlib

Utilizada para gerar os gráficos das simulações e das análises.

```python
import matplotlib.pyplot as plt
```

### Bibliotecas da linguagem

Também são utilizados recursos próprios do Python, como:

```python
import os
```

O módulo `os` é utilizado para criar automaticamente a pasta onde os resultados serão armazenados.

---

# 📚 O problema

O problema estudado consiste em analisar a velocidade de uma esfera que se movimenta em um fluido sob ação da gravidade.

Para facilitar a análise, as variáveis físicas são transformadas em variáveis **adimensionais**.

A velocidade característica utilizada é a velocidade de Stokes:

```math
U_s =
\frac{2a^2\Delta\rho\,g}{9\eta}
```

onde:

- `a` — raio da partícula;
- `Δρ` — diferença entre a densidade da partícula e a densidade do fluido;
- `g` — aceleração da gravidade;
- `η` — viscosidade dinâmica do fluido.

A velocidade da partícula é então representada pela variável adimensional:

```math
v^{*} =
\frac{v_z}{U_s}
```

e o tempo adimensional por:

```math
t^{*} =
\frac{tU_s}{a}
```

---

# 🔢 Parâmetros adimensionais

## Número de Stokes

O número de Stokes utilizado no problema é:

```math
St =
\frac{2\rho_s aU_s}{9\eta}
```

Ele está relacionado aos efeitos de inércia da partícula em comparação com o arrasto viscoso.

## Número de Reynolds

O número de Reynolds da partícula é:

```math
Re_s =
\frac{\rho_fU_sa}{\eta}
```

Esse parâmetro permite analisar a influência dos efeitos inerciais do fluido.

---

# 🧮 Equação utilizada

## Regime de Stokes

Quando o número de Reynolds tende a zero, o problema pode ser representado pela equação:

```math
\frac{dv^{*}}{dt^{*}}
=
\frac{1-v^{*}}{St}
```

A condição inicial utilizada é:

```math
v^{*}(0)=0
```

Nesse caso, existe uma solução analítica que pode ser utilizada para verificar se o método numérico está funcionando corretamente:

```math
v^{*}(t^{*})
=
1-e^{-t^{*}/St}
```

Essa solução é especialmente importante porque permite comparar diretamente o resultado obtido pelo RK4 com uma solução conhecida.

---

# 🌊 Regime de Oseen

Para considerar a influência de um número de Reynolds diferente de zero, é utilizada a correção de Oseen:

```math
\frac{dv^{*}}{dt^{*}}
=
\frac{
1-v^{*}-\frac{3}{8}Re_s(v^{*})^2
}{St}
```

com condição inicial:

```math
v^{*}(0)=0
```

Podemos perceber que, em relação à equação de Stokes, aparece um novo termo:

```math
\frac{3}{8}Re_s(v^{*})^2
```

Esse termo representa a correção associada ao efeito do número de Reynolds.

Quando:

```text
Re_s = 0
```

a equação volta a ser a equação do regime de Stokes.

---

# ⚙️ Como o RK4 foi implementado

O método utilizado para resolver numericamente a equação foi o **Runge-Kutta de 4ª ordem**.

Para uma equação geral:

```math
\frac{dv}{dt}=f(t,v)
```

o método calcula quatro valores intermediários.

### Primeiro passo

```math
k_1=f(t_i,v_i)
```

### Segundo passo

```math
k_2=
f\left(
t_i+\frac{h}{2},
v_i+\frac{hk_1}{2}
\right)
```

### Terceiro passo

```math
k_3=
f\left(
t_i+\frac{h}{2},
v_i+\frac{hk_2}{2}
\right)
```

### Quarto passo

```math
k_4=
f(t_i+h,v_i+hk_3)
```

Depois dos quatro cálculos, a velocidade é atualizada por:

```math
v_{i+1}
=
v_i+
\frac{h}{6}
\left(
k_1+2k_2+2k_3+k_4
\right)
```

O método possui erro global da ordem de:

```math
O(h^4)
```

Por isso, uma das análises do trabalho consiste justamente em verificar se o comportamento do erro acompanha essa ordem.

---

# 💻 Organização do código

O código principal está no arquivo:

```text
PPC1_222005878.py
```

Dentro dele, o problema é dividido em funções para facilitar a organização.

Entre as principais funções utilizadas estão:

```python
f()
```

Função responsável por representar a equação diferencial.

```python
rk4()
```

Função responsável por realizar a integração numérica utilizando o método RK4.

```python
solucao_analitica_stokes()
```

Função utilizada para calcular a solução analítica no regime de Stokes.

Dessa forma, o código pode ser organizado em três partes principais:

```text
Equação diferencial
        ↓
      RK4
        ↓
Resultados e gráficos
```

---

# 📁 Estrutura do repositório

A organização do projeto é:

```text
.
├── PPC1_222005878.py
│
├── resultados_ppc1/
│   ├── tarefa1_stokes_validacao.png
│   ├── tarefa2_convergencia_h.png
│   └── tarefa35_efeito_reynolds.png
│
└── README.md
```

O arquivo Python contém o programa principal.

A pasta `resultados_ppc1` contém os gráficos gerados durante a execução.

O `README.md` contém a explicação do problema, dos métodos utilizados e dos resultados.

---

# 🚀 Como executar

## 1. Instalar o Python

Primeiramente, é necessário possuir o **Python 3.8 ou superior** instalado.

Para verificar a versão instalada, execute:

```bash
python --version
```

---

## 2. Instalar as bibliotecas

O projeto utiliza principalmente `numpy` e `matplotlib`.

Para instalar as duas bibliotecas:

```bash
pip install numpy matplotlib
```

---

## 3. Baixar ou clonar o repositório

Depois de obter o projeto, entre na pasta onde está o arquivo:

```text
PPC1_222005878.py
```

---

## 4. Executar o programa

No terminal, execute:

```bash
python PPC1_222005878.py
```

O programa realizará automaticamente as análises e criará a pasta:

```text
resultados_ppc1/
```

com os gráficos gerados.

---

# 📊 Análises realizadas

O programa foi utilizado para realizar três análises principais.

---

## 1️⃣ Validação do RK4 no regime de Stokes

Primeiramente, é necessário verificar se o método RK4 está produzindo resultados confiáveis.

Para isso, o programa resolve numericamente a equação de Stokes e compara o resultado com a solução analítica:

```math
v^{*}(t^{*})
=
1-e^{-t^{*}/St}
```

Foram utilizados os valores:

```text
St = [0.25, 0.5, 1.0, 2.0]
```

e:

```text
h = 0.01
t_max = 8.0
Re_s = 0
```

Como `Re_s = 0`, estamos trabalhando no regime de Stokes.

O resultado dessa análise é salvo em:

```text
resultados_ppc1/tarefa1_stokes_validacao.png
```

### O que essa análise mostra?

O gráfico permite visualizar duas curvas:

- solução numérica obtida pelo RK4;
- solução analítica.

Se as duas curvas estiverem praticamente sobrepostas, isso indica que a implementação do RK4 está apresentando um bom resultado para esse problema.

---

# 2️⃣ Análise de convergência

Depois de verificar a solução, é interessante analisar o comportamento do método quando o tamanho do passo `h` é alterado.

Foram testados diferentes valores de `h`:

```text
h = [1.0, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001]
```

Para essa análise foram utilizados:

```text
St = 1.0
Re_s = 0
t_max = 5.0
```

Para cada valor de `h`, o programa:

1. resolve a equação utilizando RK4;
2. calcula a solução analítica;
3. compara as duas soluções;
4. calcula o erro máximo;
5. armazena o resultado.

O erro máximo utilizado é:

```math
E_{\max}
=
\max_j
\left|
v_{\mathrm{RK4},j}
-
v_{\mathrm{exato},j}
\right|
```

O resultado é apresentado em um gráfico log-log.

O gráfico é salvo em:

```text
resultados_ppc1/tarefa2_convergencia_h.png
```

### O que esperamos observar?

Como o RK4 possui erro global de ordem:

```math
O(h^4)
```

espera-se que, ao diminuir o valor de `h`, o erro também diminua aproximadamente seguindo uma relação de quarta ordem.

Essa análise serve para verificar numericamente uma das principais características do método RK4.

---

# 3️⃣ Influência do número de Reynolds

Na terceira análise, o objetivo é observar o que acontece quando o número de Reynolds aumenta.

Foram utilizados:

```text
Re_s = [0.0, 0.25, 0.5, 1.0, 2.0]
```

com:

```text
St = 1.0
h = 0.005
t_max = 8.0
```

Para cada valor de `Re_s`, o programa resolve numericamente a equação correspondente.

Quando:

```text
Re_s = 0
```

é utilizado o modelo de Stokes.

Para valores maiores de `Re_s`, é utilizado o modelo de Oseen:

```math
\frac{dv^{*}}{dt^{*}}
=
\frac{
1-v^{*}-\frac{3}{8}Re_s(v^{*})^2
}{St}
```

O gráfico dessa análise é salvo em:

```text
resultados_ppc1/tarefa35_efeito_reynolds.png
```

### O que podemos observar?

Aumentando `Re_s`, o termo:

```math
\frac{3}{8}Re_s(v^{*})^2
```

passa a ter maior influência na equação.

Com isso, a solução se afasta progressivamente daquela obtida no limite de Stokes.

---

# 📈 Resumindo o que foi feito

O trabalho pode ser entendido através do seguinte fluxo:

```text
Problema físico
      ↓
Adimensionalização
      ↓
Equação diferencial
      ↓
Implementação do RK4
      ↓
Solução numérica
      ↓
Comparação com solução analítica
      ↓
Análise de convergência
      ↓
Análise do efeito de Reynolds
      ↓
Geração dos gráficos
```

Dessa maneira, o projeto não utiliza o RK4 apenas para obter uma solução numérica. Ele também verifica o comportamento do método e analisa como os parâmetros do problema influenciam a solução.

---

# 📐 Principais parâmetros utilizados

| Parâmetro | Descrição |
|---|---|
| `St` | Número de Stokes |
| `Re_s` | Número de Reynolds da partícula |
| `h` | Passo de integração |
| `t_max` | Tempo máximo da simulação |
| `U_s` | Velocidade de Stokes |
| `v*` | Velocidade adimensional |
| `t*` | Tempo adimensional |

---

# 📂 Resultados

Após executar o programa, os resultados ficam organizados na pasta:

```text
resultados_ppc1/
```

Os arquivos gerados são:

```text
tarefa1_stokes_validacao.png
```

Validação da solução numérica com a solução analítica.

```text
tarefa2_convergencia_h.png
```

Análise da convergência do método conforme o passo `h` é reduzido.

```text
tarefa35_efeito_reynolds.png
```

Análise da influência do número de Reynolds na solução.

---

# 📚 Referências

1. Sobral, Y. D., Oliveira, T. F., Cunha, F. R. — *On the unsteady forces during the motion of a sedimenting particle*, Powder Technology, 178 (2007), 129–141.

2. Chapra, S. C.; Canale, R. P. — *Métodos Numéricos para Engenharia*, McGraw-Hill, 5ª ed., 2008.
