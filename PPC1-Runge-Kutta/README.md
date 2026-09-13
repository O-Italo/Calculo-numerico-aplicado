## Programa de casa numero 1
o código tem o objetivo de resolver numericamente pelo método de runge-kutta de 4 ordem
implementando do 0, sem bibliotecas, apenas as necessárias para criar gráficos e pastas
equação do movimento de uma esfera sedimentando em um fluido viscoso foi dado no rotéiro
em duas formulações adimensionais:

Regime de Stokes (Re → 0): St·dv/dt = 1 − v
Regime de Oseen (Re ≠ 0): St·dv/dt = 1 − v − (3/8)·Re·v²

## Requisitos
- **Python 3.8+**
- `numpy` 
- `matplotlib` 

## Como rodar: escreva no terminal python3 PPC1_222005878.py
O script cria automaticamente a pasta resultados_ppc1/ com três gráficos:

			
	
	
	
| Arquivo	| Conteúdo | 
|---|---|
| tarefa1_stokes_validacao.png |RK4 vs. solução exata no regime de Stokes, para vários números de Stokes (St)|
| tarefa2_convergencia_h.png | Erro máximo do RK4 em função do passo de tempo h (verifica a convergência de 4ª ordem) |
| tarefa35_efeito_reynolds.png | RK4 vs. solução exata no regime de Oseen, para vários números de Reynolds (Re) |
