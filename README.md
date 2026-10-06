# REMD-GUI.jl

Interface gráfica interativa para estudar **Replica Exchange Molecular Dynamics (T-REMD)** num sistema simples: um fluido de Lennard-Jones em 2D. Feita para ensino: cada parâmetro de controle do método pode ser alterado e o efeito na amostragem aparece nos gráficos.

Inspirada na interface do [FundamentosDMC.jl](https://github.com/m3g/FundamentosDMC.jl).

## O que a GUI mostra

| Painel | Pergunta que responde |
| --- | --- |
| Painéis de réplicas (direita) | Onde está cada configuração agora? A última troca foi aceita? |
| Random walk | Cada configuração percorre todos os kT? |
| Energia potencial por réplica | O sistema equilibrou? As réplicas se separam por temperatura? |
| P(U) por réplica | As distribuições de energia de kT vizinhos se sobrepõem? |
| Réplica de referência | Qual configuração está no menor kT? |
| Mapa: ocupação | Fração do tempo de cada walker em cada kT |
| Mapa: menor energia | Em que kT está a configuração de menor U |
| Mapa: taxa de troca (kT) | Trocas aceitas / tentadas por par de kT vizinhos |
| Mapa: trocas entre walkers | Há grupos de configurações isolados? |

## Instalação

Requer [Julia](https://julialang.org/downloads/) e suporte a OpenGL 3.3 (exigência do GLMakie). Testado em Ubuntu 20.04.

```bash
git clone https://github.com/LP-IQ/REMD-GUI.jl.git
cd REMD-GUI.jl
julia -e 'using Pkg; Pkg.add("GLMakie")'
```

## Uso

```bash
julia -t auto remd_gui.jl
```

`-t auto` usa uma thread por núcleo; cada réplica é propagada numa thread. A primeira execução demora alguns minutos, porque o Julia compila o GLMakie.

No REPL:

```julia
include("remd_gui.jl")
remd_gui()
```

Sequência básica:

1. Ajuste o número de réplicas e os limites de kT e clique em **Gerar escada geométrica** (ou edite a lista `kTs` à mão).
2. Clique em **1. Minimizar + Equilibrar**.
3. Clique em **2. Run**. Um novo clique continua de onde parou.
4. Use **Salvar figuras** para gravar um PNG de cada gráfico e **Salvar trajetória do menor kT** para gerar um `.xyz` e um script `.tcl` (abrir com `vmd -e arquivo.tcl`).

## Modelo

- Partículas de Lennard-Jones em 2D, caixa quadrada periódica, sem raio de corte (todos os pares, imagem mínima).
- Unidades reduzidas: massa 1, energia em unidades de epsilon, temperatura como kT.
- Dinâmica de Langevin (integrador BAOAB) em cada réplica.
- Tentativas de troca entre kT vizinhos, alternando pares pares e ímpares, com o critério

```math
P = \min\left\{1,\ \exp\left[\left(\frac{1}{kT_i}-\frac{1}{kT_j}\right)\left(U_i-U_j\right)\right]\right\}
```

- Após uma troca aceita, as velocidades são reescaladas por √(kT novo / kT antigo).

## Documentação

[docs/guia_parametros.md](docs/guia_parametros.md) explica cada parâmetro, sua importância no REMD, como ler cada painel e um roteiro de experimentos.

## Limitações

- O cálculo de forças é O(N²); acima de algumas centenas de partículas a simulação fica lenta.
- É uma ferramenta didática, não um código de produção.

## Autores

Lucas Pinheiro (IQ-UNICAMP).
