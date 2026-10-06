# REMD na GUI: guia de parâmetros

## O que a GUI simula

A GUI roda um T-REMD: várias cópias do mesmo sistema, cada uma numa temperatura fixa, que de tempos em tempos tentam trocar de configuração com a vizinha.

O sistema é o mesmo do FundamentosDMC: partículas de Lennard-Jones em 2D, numa caixa quadrada periódica. As unidades são reduzidas: massa 1, energia em unidades de epsilon (a profundidade do poço de Lennard-Jones), temperatura dada como kT.

Dois termos aparecem o tempo todo e não são a mesma coisa:

- **Estado**: uma temperatura da escada. Cada painel da direita é um estado, e o kT dele nunca muda.
- **Walker** (configuração): um conjunto de posições e velocidades. É o que tem cor. Quando uma troca é aceita, dois walkers trocam de painel e levam a cor junto.

O objetivo do método é amostrar bem o estado mais frio. Os estados quentes cruzam barreiras com facilidade e, pelas trocas, entregam configurações novas ao frio. Cada estado continua amostrando o ensemble canônico da sua temperatura, porque o critério de troca respeita o balanço detalhado.

Neste sistema, o comportamento muda com kT. Em kT baixo as partículas se juntam em aglomerados e a energia potencial é bem negativa. Em kT alto elas se espalham como um gás e a energia fica perto de zero. A região intermediária, onde os aglomerados se formam e se desfazem, é a mais difícil de atravessar na escada.

## A escada de temperaturas

A escada é a principal variável de controle do REMD: ela decide se as trocas acontecem e se as réplicas quentes são quentes o bastante para ajudar.

### kT mín

É a temperatura do estado de interesse, aquele que você quer amostrar. Não é um parâmetro de ajuste: ele vem da pergunta física. Quanto mais baixo, mais o sistema fica preso em aglomerados e mais o REMD tem a oferecer.

### kT máx

É a temperatura do estado mais quente. Ele precisa ser alto o bastante para o sistema esquecer rapidamente onde estava, ou seja, para os aglomerados se desfazerem e se refazerem em poucos passos.

- **Baixo demais**: as trocas funcionam, mas nenhuma réplica cruza as barreiras. O random walk parece saudável e a amostragem do estado frio não melhora.
- **Alto demais**: a faixa a cobrir cresce, e são necessárias mais réplicas para manter a aceitação.

### Réplicas

É o número de estados entre kT mín e kT máx. Com a faixa fixa, mais réplicas significam vizinhos mais próximos e aceitação maior, mas também mais custo e mais degraus para um walker percorrer.

O número necessário cresce com o tamanho do sistema. A diferença de energia média entre dois estados cresce com n, enquanto a largura das flutuações cresce só com a raiz de n. Para manter a sobreposição, o espaçamento tem de encolher como 1/√n.

### kTs e o botão "Gerar escada geométrica"

O campo kTs é a lista que a simulação usa de fato. O botão preenche essa lista com uma progressão geométrica entre kT mín e kT máx:

```math
kT_i = kT_{\min}\left(\frac{kT_{\max}}{kT_{\min}}\right)^{\frac{i-1}{R-1}}
```

A progressão geométrica mantém constante a razão entre temperaturas vizinhas. Ela dá aceitação uniforme quando a capacidade calorífica não muda ao longo da escada.

Neste sistema ela muda: perto da condensação a energia varia muito com kT, e os pares frios aceitam menos que os quentes. A lista pode ser editada à mão para aproximar os valores nessa região e afastá-los na região quente. Depois de mudar a lista, é preciso clicar em Reset ou no botão 1.

### Como julgar a escada

A faixa-alvo usual de aceitação por par é de 0,2 a 0,3; até 0,4 é aceitável. Abaixo de 0,1 em algum par, a escada tem um gargalo e os walkers não o atravessam. Acima de 0,6 em todos, há réplicas sobrando: dá para tirar algumas ou esticar o kT máx.

Nos testes com kT de 0,55 a 0,86, 6 réplicas e n = 100, a aceitação ficou entre 0,5 e 0,8, crescendo do par mais frio para o mais quente. Essa escada está mais densa do que precisa.

## O sistema: n, L, epsilon e sigma

Esses quatro parâmetros definem o que está sendo simulado. No REMD eles importam porque fixam a escala de energia e o tamanho das flutuações, e portanto a aceitação.

### n partículas

É o parâmetro do sistema que mais afeta o REMD. Com mais partículas, as distribuições de energia de estados vizinhos se afastam mais rápido do que se alargam, e a aceitação cai para a mesma escada. Dobrar n pede cerca de √2 vezes mais réplicas na mesma faixa.

O custo por passo cresce com n², porque a GUI calcula todos os pares.

### L (caixa)

O lado da caixa quadrada. Junto com n, define a densidade, n/L². Com n = 100 e L = 100 o sistema é diluído: as partículas passam a maior parte do tempo longe umas das outras, e os aglomerados só se formam em kT baixo.

Uma caixa menor aumenta a densidade, torna a energia mais negativa em todas as temperaturas e desloca a região de condensação para kT mais alto.

### epsilon (poço LJ)

A profundidade do poço de Lennard-Jones, ou seja, a energia ganha quando duas partículas se encostam. O que importa fisicamente é a razão kT/epsilon. Dobrar epsilon tem o mesmo efeito que dividir por dois todas as temperaturas da escada.

### sigma (diâmetro LJ)

O diâmetro das partículas. O mínimo do potencial fica em 1,12 × sigma. Define o tamanho dos aglomerados na tela e, com n e L, a fração da caixa ocupada.

O potencial usado é:

```math
U(r) = 4\,\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]
```

A GUI não usa raio de corte: todos os pares entram, com a imagem periódica mais próxima. Isso mantém a energia usada no critério de troca exata.

## Integração e termostato: dt e gamma

Esses dois parâmetros controlam como cada réplica anda entre uma tentativa de troca e a seguinte. Cada réplica é propagada com dinâmica de Langevin, que é o que a mantém no seu kT.

### dt

O passo de integração. Um dt maior cobre mais tempo simulado por passo, mas, acima de um limite, a integração fica instável e a energia explode.

No REMD há um cuidado a mais: o limite é dado pela réplica mais quente, onde as partículas são mais rápidas e as colisões mais fortes. Um dt que funciona em kT = 0,35 pode falhar em kT = 1,5. Se a energia de um painel quente disparar, reduza dt.

Nos testes com a GUI, dt = 0,05 ficou estável até kT = 0,86. Acima disso, confira antes de confiar.

### gamma (atrito)

O coeficiente de atrito do termostato de Langevin. Ele define a força do acoplamento com o banho térmico: o tempo de relaxação da temperatura é da ordem de 1/gamma.

- **Pequeno** (0,01): a dinâmica é quase newtoniana e a temperatura demora a se ajustar.
- **Grande** (1 ou mais): a temperatura se ajusta rápido, mas o movimento vira difusivo e o sistema explora o espaço mais devagar.

No REMD isso importa logo depois de uma troca. O walker chega a uma temperatura nova com as velocidades reescaladas, mas as posições ainda são típicas da temperatura anterior. Ele precisa de um tempo da ordem de 1/gamma para relaxar. Com gamma = 0,1 e dt = 0,05, isso dá cerca de 200 passos.

## As trocas

### O critério

A cada tentativa, dois estados vizinhos i e j comparam suas energias potenciais. A troca é aceita com probabilidade:

```math
P = \min\left\{1,\ \exp\left[\left(\frac{1}{kT_i} - \frac{1}{kT_j}\right)\left(U_i - U_j\right)\right]\right\}
```

Se a configuração que está no estado quente tem energia menor que a do estado frio, a troca é sempre aceita. No caso contrário, a probabilidade cai exponencialmente com a diferença de energia e com a distância entre as temperaturas.

É por isso que a aceitação depende da sobreposição dos histogramas P(U): só há troca quando as duas réplicas podem ter energias parecidas.

Quando a troca é aceita, as velocidades são multiplicadas por √(kT novo / kT antigo), para que a energia cinética já corresponda à nova temperatura.

As tentativas alternam entre dois conjuntos de pares: numa rodada 1-2, 3-4, 5-6; na seguinte 2-3, 4-5. Cada par é testado uma vez a cada duas rodadas.

### Troca a cada (passos)

O número de passos de dinâmica entre duas rodadas de tentativas.

- **Pequeno**: mais tentativas por tempo simulado, e os walkers percorrem a escada mais rápido. Tentar com mais frequência não viola o balanço detalhado.
- **Grande**: cada réplica tem mais tempo para relaxar e explorar na temperatura nova antes da próxima tentativa, mas a escada é percorrida devagar.

O valor padrão, 50 passos, é menor que o tempo de relaxação do termostato (cerca de 200 passos). Isso é aceitável para amostragem, mas vale testar 200 e comparar o número de round trips pelo mesmo número total de passos.

### Pausa se trocou (s)

Só afeta a visualização. Quando alguma troca é aceita, a tela congela por esse tempo para dar para ver as cores mudando de painel. Não altera a simulação. Use 0 para rodar na velocidade máxima.

## Preparação e duração

### Botão 1: Minimizar + Equilibrar

A configuração inicial tem posições sorteadas, com partículas sobrepostas e energia enorme. O botão faz duas coisas em sequência:

1. **Minimização**: desloca as partículas na direção das forças até remover as sobreposições. Sem isso, o primeiro passo de dinâmica explode.
2. **Equilibração**: copia a configuração minimizada para todas as réplicas, sorteia velocidades na temperatura de cada uma e roda a dinâmica sem trocas.

### Equilibração (passos)

O número de passos dessa etapa sem trocas. Todas as réplicas partem da mesma configuração, então no início nenhuma está em equilíbrio com seu kT. Se as trocas começarem cedo demais, as primeiras refletem essa relaxação inicial e não o comportamento da escada.

No gráfico de energia, a equilibração está boa quando as curvas param de derivar e se separam por temperatura. As réplicas frias demoram mais, porque formar aglomerados é lento. Para kT mín abaixo de 0,4, 2000 passos provavelmente é pouco. A equilibração aparece no gráfico em passos negativos, e a linha tracejada em 0 marca o início das trocas.

### nsteps

O número de passos de cada clique em Run. Um novo clique continua de onde parou, então dá para acumular estatística aos poucos.

A duração deve ser julgada pelo número de round trips (idas e voltas), não pelo número de passos. Com poucos round trips, o estado frio recebeu poucas configurações independentes e as médias ainda não são confiáveis.

### Semente

A semente do gerador de números aleatórios, que define as posições iniciais, as velocidades, o ruído do termostato e os sorteios das trocas. Com 0, uma semente nova é sorteada a cada Reset. Com qualquer outro valor, a simulação se repete.

Para comparar duas escadas de forma justa, use a mesma semente nas duas.

## Como ler a janela

| Painel | O que mostra | Sinal de que está bom | Sinal de problema |
| --- | --- | --- | --- |
| Random walk | Em que kT cada walker está, a cada tentativa | Cada linha sobe até o topo e desce até a base várias vezes | Linhas presas numa faixa; um degrau que ninguém cruza |
| Energia potencial por réplica | Média móvel de U de cada estado, de azul (menor kT) a vermelho (maior kT) | Curvas ordenadas por temperatura, sem deriva | Curva fria ainda caindo: falta equilibração |
| P(U) por réplica | Histograma da energia de cada estado | Vizinhos se sobrepõem em parte | Vizinhos sem sobreposição (sem trocas) ou quase idênticos (réplicas sobrando) |
| Réplica de referência | A configuração que está no menor kT | A cor muda com frequência | Sempre a mesma cor: o frio não recebe configurações novas |
| Painéis da direita | Cada estado, com kT, U e a probabilidade da última tentativa | Títulos verdes ("Trocou") aparecem em todos os pares | Um par sempre vermelho ("Falhou") |
| Taxa de troca acumulada | Trocas aceitas / tentativas, por par | Entre 0,2 e 0,3, parecida em todos os pares | Um par abaixo de 0,1 |
| Round trips (idas e voltas) | Quantas vezes cada walker foi do menor ao maior kT e voltou | Cresce de forma parecida para todos os walkers | Zero ou muito desigual entre walkers |

A aceitação sozinha não basta. Uma escada com kT máx baixo tem aceitação ótima e muitos round trips, e mesmo assim não ajuda, porque o topo não é quente o bastante para mudar as configurações. Olhe também o painel do kT mais quente: os aglomerados precisam se desfazer ali.

### O seletor Mapa

O painel inferior direito mostra um de quatro mapas, escolhido no menu "Mapa". Cada um responde a uma pergunta diferente.

| Mapa | O que cada célula mostra | Pergunta que responde | Sinal de problema |
| --- | --- | --- | --- |
| Ocupação | Fração do tempo de cada walker em cada kT | Os walkers percorrem a escada toda? | Blocos; células longe de 1/R |
| Menor energia | Fração do tempo em que a configuração de menor U estava naquele kT e era aquele walker | Em que kT fica a menor energia? | Sempre o mesmo walker: falta mistura |
| Taxa de troca (kT) | Trocas aceitas / tentadas entre dois kT vizinhos | A escada está bem espaçada? | Um par perto de 0, ou todos acima de 0,5 |
| Trocas entre walkers | Fração do total de trocas aceitas que ocorreu entre dois walkers | Há grupos isolados? | Blocos: walkers que só trocam entre si |

No mapa "Taxa de troca (kT)" só as células vizinhas à diagonal ficam preenchidas, porque apenas kT vizinhos tentam trocar.

A menor energia aparece às vezes num kT mais alto. Isso é esperado: as distribuições P(U) se sobrepõem, e é essa sobreposição que permite as trocas. O REMD não procura o mínimo de energia; ele amostra o ensemble canônico em cada kT.

### Marcações e cores

- **Moldura magenta**: nos painéis da direita, marca o estado que contém a configuração de menor U naquele instante.
- **Chave "menor U no random walk"**: quando ligada, pontos magenta no random walk mostram em que kT estava a menor U a cada tentativa.
- **Menu "Cores"**: "Destaque" pinta um walker de vermelho e os outros de cinza, para seguir uma configuração. "Colorido" dá uma cor a cada walker.

## Salvar figuras e trajetória

- **Salvar figuras (PNG por gráfico)**: cria uma pasta com um PNG de cada gráfico, os quatro mapas e a janela inteira. O nome da pasta traz réplicas, kT, n, intervalo de troca e passo. A resolução é a da janela, então maximize antes de salvar.
- **Salvar trajetória do menor kT (XYZ)**: grava um frame a cada 10 passos do Run, sempre do estado de menor kT, e um script `.tcl` para abrir no VMD com `vmd -e arquivo.tcl`.

A trajetória é contínua no estado, não na configuração: a cada troca aceita no par mais frio as posições saltam. A linha de comentário de cada frame diz qual walker estava ali.

## Roteiro de experimentos para o tutorial

Cada experimento muda um parâmetro e mantém o resto. Use a mesma semente (por exemplo, 42) em todos, para que a diferença venha do parâmetro. Os resultados esperados abaixo são previsões; só o primeiro foi observado.

1. **Referência.** 6 réplicas, kT de 0,55 a 0,86. Observado: aceitação de 0,5 a 0,8 e cerca de 40 round trips em 10 000 passos. Serve para apresentar os painéis.
2. **Menos réplicas na mesma faixa.** 3 réplicas, kT de 0,55 a 0,86. Esperado: aceitação menor, histogramas com menos sobreposição. Mostra que o número de réplicas controla a aceitação.
3. **Faixa larga com poucas réplicas.** 3 réplicas, kT de 0,35 a 1,2. Esperado: histogramas separados, aceitação perto de zero, mapa de ocupação em blocos. É o caso em que o REMD vira simulações independentes.
4. **A mesma faixa larga com mais réplicas.** 8 a 10 réplicas, kT de 0,35 a 1,2. Esperado: as trocas voltam, com os pares frios aceitando menos que os quentes.
5. **Escada ajustada à mão.** Partindo do experimento 4, aproxime os kTs da região fria e afaste os da quente, editando o campo kTs. Esperado: aceitação mais uniforme e mais round trips com o mesmo número de réplicas.
6. **Frequência de troca.** Na escada do experimento 5, compare "troca a cada" 10, 50 e 500 passos, com o mesmo nsteps. Esperado: menos round trips quanto maior o intervalo.
7. **Tamanho do sistema.** Na escada do experimento 5, mude n de 100 para 200, com L = 141 para manter a densidade. Esperado: a aceitação cai em todos os pares.
8. **Gargalo (opcional).** Uma escada com um buraco no meio, por exemplo 0,55, 0,58, 0,61, 0,80, 0,83, 0,86. Esperado: taxa de troca baixa no par do meio e blocos no mapa "Trocas entre walkers".

Para cada experimento, use o botão "Salvar figuras" ao final da corrida e anote a taxa de troca por par e o total de round trips.

## Resumo dos parâmetros

| Parâmetro | Padrão | O que controla | Efeito no REMD |
| --- | --- | --- | --- |
| réplicas | 6 | Número de temperaturas | Mais réplicas: aceitação maior, custo maior |
| kT mín | 0,35 | Temperatura do estado de interesse | Define o problema; quanto mais frio, mais difícil |
| kT máx | 0,86 | Temperatura do estado mais quente | Tem de ser alto o bastante para desfazer os aglomerados |
| kTs | lista | As temperaturas usadas de fato | Espaçamento entre vizinhos define a aceitação de cada par |
| n partículas | 100 | Tamanho do sistema | Mais partículas: aceitação menor na mesma escada |
| L (caixa) | 100 | Densidade, com n | Muda onde ocorre a condensação |
| epsilon | 1,0 | Profundidade do poço | Só a razão kT/epsilon importa |
| sigma | 2,0 | Diâmetro das partículas | Tamanho dos aglomerados e fração ocupada |
| dt | 0,05 | Passo de integração | Limitado pela réplica mais quente |
| gamma (atrito) | 0,1 | Acoplamento do termostato | Tempo de relaxação após uma troca, cerca de 1/gamma |
| troca a cada | 50 | Passos entre tentativas | Menor: escada percorrida mais rápido |
| pausa se trocou | 0,3 s | Só a visualização | Nenhum |
| equilibração | 2000 | Passos sem trocas no início | Pouca: primeiras trocas enviesadas |
| nsteps | 20 000 | Duração de cada Run | Julgar pelo número de round trips |
| semente | 0 | Números aleatórios | Valor fixo torna a simulação repetível |
