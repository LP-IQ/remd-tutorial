#!/usr/bin/env julia
# REMD interativo (T-REMD) em um fluido de Lennard-Jones 2D
#
# Instalação (uma vez):   julia -e 'using Pkg; Pkg.add("GLMakie")'
# Execução:               julia -t auto remd_gui.jl
#   ou, no REPL:          include("remd_gui.jl"); remd_gui()
#
# Unidades reduzidas: massa = 1, energia em unidades de epsilon (poço de Lennard-Jones), temperatura como kT.
# Cada painel da direita é um ESTADO (kT fixo); a cor é a identidade da CONFIGURAÇÃO (walker).
# Troca aceita = duas configurações trocam de painel.
# Documentação: README.md e docs/tutorial.md

using GLMakie, Random, Printf

const MAXREP = 12

# ----------------------------------------------------------------------------------------
# Núcleo da simulação
# ----------------------------------------------------------------------------------------
mutable struct Replica
    x::Vector{Float64}; y::Vector{Float64}       # posições
    vx::Vector{Float64}; vy::Vector{Float64}     # velocidades
    fx::Vector{Float64}; fy::Vector{Float64}     # forças
    U::Float64                                   # energia potencial (unidades de epsilon)
    walker::Int                                  # identidade da configuração (cor)
end
Replica(n::Int, w::Int) = Replica(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), 0.0, w)

Base.@kwdef mutable struct Sim
    reps::Vector{Replica}
    kTs::Vector{Float64}
    L::Float64 = 100.0
    eps::Float64 = 1.0
    sig::Float64 = 2.0
    dt::Float64 = 0.05
    gam::Float64 = 0.1
    step::Int = 0
    nattempt::Int = 0
    att::Vector{Int} = Int[]              # tentativas por par
    acc::Vector{Int} = Int[]              # trocas aceitas por par
    last_prob::Vector{Float64} = Float64[]  # P da última tentativa de cada par (NaN = nunca tentado)
    rngs::Vector{MersenneTwister} = MersenneTwister[]
    prepared::Bool = false
end

function forces!(r::Replica, L, eps, sig)
    n = length(r.x); fill!(r.fx, 0.0); fill!(r.fy, 0.0); U = 0.0; s2 = sig^2
    @inbounds for i in 1:n-1
        xi = r.x[i]; yi = r.y[i]
        for j in i+1:n
            dx = xi - r.x[j]; dy = yi - r.y[j]
            dx -= L * round(dx / L); dy -= L * round(dy / L)
            r2 = dx * dx + dy * dy
            s6 = (s2 / r2)^3
            U += 4 * eps * (s6 * s6 - s6)
            fm = 24 * eps * (2 * s6 * s6 - s6) / r2
            r.fx[i] += fm * dx; r.fy[i] += fm * dy
            r.fx[j] -= fm * dx; r.fy[j] -= fm * dy
        end
    end
    r.U = U
    return U
end

function minimize!(r::Replica, L, eps, sig; niter=500)
    n = length(r.x)
    for _ in 1:niter
        forces!(r, L, eps, sig)
        @inbounds for i in 1:n
            f = sqrt(r.fx[i]^2 + r.fy[i]^2)
            f < 1e-12 && continue
            d = min(0.05, 0.01 * f)
            r.x[i] = mod(r.x[i] + d * r.fx[i] / f, L)
            r.y[i] = mod(r.y[i] + d * r.fy[i] / f, L)
        end
    end
    return forces!(r, L, eps, sig)
end

# Langevin BAOAB na temperatura kT
function langevin_step!(r::Replica, kT, dt, gam, L, eps, sig, rng)
    n = length(r.x); c1 = exp(-gam * dt); c2 = sqrt(kT * (1 - c1^2)); h = 0.5 * dt
    @inbounds for i in 1:n
        r.vx[i] += h * r.fx[i]; r.vy[i] += h * r.fy[i]
        r.x[i] += h * r.vx[i]; r.y[i] += h * r.vy[i]
        r.vx[i] = c1 * r.vx[i] + c2 * randn(rng)
        r.vy[i] = c1 * r.vy[i] + c2 * randn(rng)
        r.x[i] = mod(r.x[i] + h * r.vx[i], L)
        r.y[i] = mod(r.y[i] + h * r.vy[i], L)
    end
    forces!(r, L, eps, sig)
    @inbounds for i in 1:n
        r.vx[i] += h * r.fx[i]; r.vy[i] += h * r.fy[i]
    end
    return nothing
end

function step_all!(s::Sim)
    Threads.@threads for i in eachindex(s.reps)
        langevin_step!(s.reps[i], s.kTs[i], s.dt, s.gam, s.L, s.eps, s.sig, s.rngs[i])
    end
    s.step += 1
    return nothing
end

# P = min{1, exp[(β_i − β_j)(U_i − U_j)]}; pares pares e ímpares alternados
function attempt_exchange!(s::Sim, rng)
    R = length(s.reps); out = Tuple{Int,Bool,Float64}[]
    for i in (1 + s.nattempt % 2):2:(R - 1)
        a = s.reps[i]; b = s.reps[i+1]
        delta = (1 / s.kTs[i] - 1 / s.kTs[i+1]) * (a.U - b.U)
        P = min(1.0, exp(delta))
        ok = rand(rng) < P
        s.att[i] += 1
        s.last_prob[i] = P
        if ok
            s.acc[i] += 1
            sc = sqrt(s.kTs[i+1] / s.kTs[i])       # v ∝ sqrt(T_novo/T_antigo)
            a.vx .*= sc; a.vy .*= sc
            b.vx ./= sc; b.vy ./= sc
            s.reps[i], s.reps[i+1] = b, a
        end
        push!(out, (i, ok, P))
    end
    s.nattempt += 1
    return out
end

function new_sim(kTs, n, L, eps, sig, dt, gam; seed=rand(1:10^6))
    R = length(kTs); rng = MersenneTwister(seed)
    x0 = L .* rand(rng, n); y0 = L .* rand(rng, n)
    reps = Replica[]
    for w in 1:R
        r = Replica(n, w); r.x .= x0; r.y .= y0; push!(reps, r)
    end
    return Sim(; reps, kTs=collect(Float64, kTs), L, eps, sig, dt, gam,
               att=zeros(Int, R - 1), acc=zeros(Int, R - 1), last_prob=fill(NaN, R - 1),
               rngs=[MersenneTwister(seed + 17 * i) for i in 1:R])
end

# Minimiza uma configuração e copia para todas as réplicas
function prepare!(s::Sim)
    R = length(s.reps); r1 = s.reps[1]
    minimize!(r1, s.L, s.eps, s.sig)
    for i in 1:R
        r = s.reps[i]; kT = s.kTs[i]
        r.x .= r1.x; r.y .= r1.y
        r.vx .= sqrt(kT) .* randn(s.rngs[i], length(r.x))
        r.vy .= sqrt(kT) .* randn(s.rngs[i], length(r.x))
        forces!(r, s.L, s.eps, s.sig)
    end
    s.prepared = true
    return nothing
end

# histograma normalizado -> pontos (centro do bin, densidade)
function hist_points(v::Vector{Float64}, lo, hi, nb)
    c = zeros(Int, nb); w = (hi - lo) / nb
    for u in v
        c[clamp(floor(Int, (u - lo) / w) + 1, 1, nb)] += 1
    end
    tot = max(1, length(v))
    return [Point2f(lo + (b - 0.5) * w, c[b] / (tot * w)) for b in 1:nb]
end

# ----------------------------------------------------------------------------------------
# Interface gráfica
# ----------------------------------------------------------------------------------------
tbtext(tb) = (d = tb.displayed_string[]; d === nothing ? "" : strip(String(d)))
getf(tb) = parse(Float64, tbtext(tb))
geti(tb) = round(Int, getf(tb))
# Testes do tutorial (docs/tutorial.md, seção 9): nome no menu => valores dos campos.
# Campos omitidos voltam ao padrão (PADRAO). A lista kTs é gerada pela escada geométrica,
# exceto quando "kts" é dado explicitamente.
const PADRAO = Dict("nrep" => "6", "kmin" => "0.55", "kmax" => "0.86", "n" => "100", "L" => "100.0",
                    "dt" => "0.05", "nsteps" => "20000", "nex" => "50", "pause" => "0", "nequil" => "2000",
                    "gam" => "0.1", "eps" => "1.0", "sig" => "2.0", "seed" => "42")
const TESTES = [
    "(escolher)"          => Dict{String,String}(),
    "1: referência"       => Dict{String,String}(),
    "2a: 3 réplicas"      => Dict("nrep" => "3"),
    "2b: 2 réplicas"      => Dict("nrep" => "2"),
    "3: larga, 3 rép."    => Dict("nrep" => "3", "kmin" => "0.35", "kmax" => "1.2", "nequil" => "5000"),
    "4a: larga, 8 rép."   => Dict("nrep" => "8", "kmin" => "0.35", "kmax" => "1.2", "nequil" => "5000"),
    "4b: larga, 10 rép."  => Dict("nrep" => "10", "kmin" => "0.35", "kmax" => "1.2", "nequil" => "5000"),
    "6a: troca a cada 10" => Dict("nex" => "10"),
    "6b: troca a cada 50" => Dict("nex" => "50"),
    "6c: troca a cada 500"=> Dict("nex" => "500"),
    "7a: n = 200"         => Dict("n" => "200", "L" => "141.0"),
    "7b: n = 400"         => Dict("n" => "400", "L" => "200.0"),
    "8: gargalo"          => Dict("kts" => "0.55,0.58,0.61,0.80,0.83,0.86"),
    "9: faixa fria"       => Dict("kmin" => "0.35", "kmax" => "0.45", "nequil" => "5000"),
]

function settext!(tb, s::String)
    tb.displayed_string[] = s; tb.stored_string[] = s
    return nothing
end

const PALETTES = Dict(
    "Destaque" => [RGBf(0.9, 0.1, 0.1); fill(RGBf(0.65, 0.65, 0.65), MAXREP - 1)],
    "Colorido" => [
        RGBf(0.95, 0.1, 0.1), RGBf(0.1, 0.5, 0.8), RGBf(0.2, 0.7, 0.2), RGBf(0.9, 0.6, 0.1),
        RGBf(0.6, 0.3, 0.8), RGBf(0.2, 0.8, 0.8), RGBf(0.9, 0.3, 0.6), RGBf(0.4, 0.4, 0.4),
        RGBf(0.7, 0.8, 0.2), RGBf(0.6, 0.4, 0.2), RGBf(0.1, 0.4, 0.4), RGBf(0.8, 0.8, 0.3)],
)
# cor por temperatura: mapa 'turbo' (azul = menor kT -> verde -> amarelo -> vermelho = maior kT)
const TURBO = GLMakie.Makie.to_colormap(:turbo)
function tempcolor(k, R)
    t = R > 1 ? (k - 1) / (R - 1) : 0.0
    c = TURBO[clamp(round(Int, 1 + (0.08 + 0.84 * t) * (length(TURBO) - 1)), 1, length(TURBO))]
    return RGBf(c.r, c.g, c.b)
end

function remd_gui()   # se ainda cortar: redimensione a janela ou edite o tamanho em Figure(size=...)
    monitor = GLMakie.GLFW.GetPrimaryMonitor()
    vidmode = GLMakie.GLFW.GetVideoMode(monitor)
    # os gráficos preenchem a janela e acompanham o redimensionamento; altura mínima útil ~650 px
    fig = Figure(size=(vidmode.width * 0.92, vidmode.height * 0.80), fontsize=12)

    left = fig[1, 1] = GridLayout(tellheight=false, valign=:top)   # não fixa a altura da janela
    center = fig[1, 2] = GridLayout()
    right_grid = fig[1, 3] = GridLayout()
    colsize!(fig.layout, 1, Fixed(300))
    colsize!(fig.layout, 3, Relative(0.24))

    # ---------------- controles ----------------
    row = Ref(0)
    function field(label, default)
        row[] += 1
        Label(left[row[], 1], label; halign=:right, fontsize=13.5)
        return Textbox(left[row[], 2]; stored_string=default, width=125, height=25, fontsize=13.5, textpadding=(7, 7, 3, 3))
    end
    tb_nrep = field("réplicas (máx $MAXREP)", "6")
    tb_kmin = field("kT mín (em epsilon)", "0.55")
    tb_kmax = field("kT máx (em epsilon)", "0.86")
    row[] += 1
    btn_ladder = Button(left[row[], 1:2]; label="Gerar escada geométrica", height=27, fontsize=13.5)
    tb_kts = field("kTs", "0.550,0.601,0.658,0.719,0.786,0.860")
    tb_n = field("n partículas", "100")
    tb_L = field("L (caixa)", "100.0")
    tb_dt = field("dt", "0.05")
    tb_nsteps = field("nsteps", "20000")
    tb_nex = field("troca a cada (passos)", "50")
    tb_pause = field("pausa se trocou (s)", "0.3")
    tb_nequil = field("equilibração", "2000")
    tb_gam = field("gamma (atrito)", "0.1")
    tb_eps = field("epsilon (poço LJ)", "1.0")
    tb_sig = field("sigma (diâmetro LJ)", "2.0")
    tb_seed = field("semente (0 = aleatória)", "0")
    row[] += 1
    Label(left[row[], 1], "Teste"; halign=:right, fontsize=13.5)
    menu_teste = Menu(left[row[], 2]; options=first.(TESTES), default="(escolher)", height=27, fontsize=13.5, width=145)
    row[] += 1
    Label(left[row[], 1], "Cores"; halign=:right, fontsize=13.5)
    menu_cores = Menu(left[row[], 2]; options=["Destaque", "Colorido"], default="Destaque", height=27, fontsize=13.5, width=145)
    row[] += 1
    Label(left[row[], 1], "Mapa"; halign=:right, fontsize=13.5)
    menu_map = Menu(left[row[], 2]; options=["Ocupação", "Menor energia", "Taxa de troca (kT)", "Trocas entre walkers"], default="Ocupação", height=27, fontsize=13.5, width=145)
    row[] += 1
    Label(left[row[], 1], "menor U no random walk"; halign=:right, fontsize=13.5)
    tg_min = Toggle(left[row[], 2]; active=false, halign=:left)
    row[] += 1
    bl = left[row[], 1:2] = GridLayout()
    bk = (; height=29, fontsize=13.5)
    btn_min = Button(bl[1, 1:3]; label="1. Minimizar + Equilibrar", buttoncolor=RGBf(0.55, 0.80, 0.55), bk...)
    btn_run = Button(bl[2, 1]; label="2. Run", buttoncolor=RGBf(1.0, 0.85, 0.2), bk...)
    btn_stop = Button(bl[2, 2]; label="Stop", buttoncolor=RGBf(0.9, 0.45, 0.45), bk...)
    btn_reset = Button(bl[2, 3]; label="Reset", buttoncolor=RGBf(0.85, 0.75, 0.9), bk...)
    btn_save = Button(bl[3, 1:3]; label="Salvar figuras (PNG por gráfico)", buttoncolor=RGBf(0.75, 0.85, 0.95), bk...)
    btn_traj = Button(bl[4, 1:3]; label="Salvar trajetória do menor kT (XYZ)", buttoncolor=RGBf(0.75, 0.85, 0.95), bk...)
    rowgap!(bl, 4); colgap!(bl, 4)
    row[] += 1
    stats = Observable("Taxa de troca: —")
    Label(left[row[], 1:2], stats; halign=:left, valign=:top, justification=:left, tellwidth=false, tellheight=false, fontsize=12.5)   # última linha absorve a sobra
    rowgap!(left, 4)

    # ---------------- centro: random walk, energia, histogramas, réplica de referência e mapa ----------------
    ax_W = Axis(center[1, 1:2]; title="Random walk das configurações", xlabel="rodada de tentativas de troca", ylabel="índice do kT (1 = menor)")
    ax_E = Axis(center[2, 1]; title="Energia potencial por réplica (média móvel)", xlabel="passo (< 0 = equilibração)", ylabel="U (unidades de epsilon)")
    ax_H = Axis(center[2, 2]; title="P(U) por réplica (azul = menor kT, vermelho = maior kT)", xlabel="U (unidades de epsilon)", ylabel="P(U)")
    ax_main = Axis(center[3, 1]; aspect=DataAspect(), title="Réplica de referência (menor kT)", xticklabelsvisible=false, yticklabelsvisible=false)
    g_O = center[3, 2] = GridLayout()
    ax_O = Axis(g_O[1, 1]; aspect=DataAspect(), title="Ocupação", xlabel="índice do kT (1 = menor)", ylabel="walker")
    rowsize!(center, 1, Relative(0.27)); rowsize!(center, 2, Relative(0.27))

    pts_main = Observable(Point2f[]); col_main = Observable(RGBf(0, 0, 0))
    scatter!(ax_main, pts_main; color=col_main, markersize=9, strokewidth=1.0, strokecolor=:black)

    epts = [Observable(Point2f[]) for _ in 1:MAXREP]      # energia por estado (kT)
    wpts = [Observable(Point2f[]) for _ in 1:MAXREP]      # random walk por walker
    hpts = [Observable(Point2f[]) for _ in 1:MAXREP]      # histograma por estado (kT)
    tcols = [Observable(tempcolor(k, MAXREP)) for k in 1:MAXREP]   # cor do estado k
    lw_e = [Observable(k == 1 ? 2.8 : 1.2) for k in 1:MAXREP]
    line_cols = [Observable(RGBf(0, 0, 0)) for _ in 1:MAXREP]
    lw_lines = [Observable(1.5) for _ in 1:MAXREP]
    alpha_W = Observable(0.85)                 # opacidade das linhas do random walk
    for k in MAXREP:-1:1
        lines!(ax_E, epts[k]; color=tcols[k], linewidth=lw_e[k])   # menor kT desenhado por último (fica por cima)
        lines!(ax_W, wpts[k]; color=line_cols[k], linewidth=lw_lines[k], alpha=alpha_W)
        lines!(ax_H, hpts[k]; color=tcols[k], linewidth=lw_e[k])
    end

    # rastreio da configuração de menor U: em qual kT ela está a cada tentativa de troca (pontos magenta, opcionais)
    minpts = Observable(Point2f[])
    COR_MENOR_U = RGBf(0.85, 0.0, 0.75)   # magenta: cor que não aparece na paleta dos walkers
    scatter!(ax_W, minpts; color=COR_MENOR_U, markersize=8, strokewidth=0, visible=tg_min.active)
    # com a chave ligada, as linhas dos walkers ficam claras para destacar os pontos da menor U
    on(tg_min.active) do ativo
        alpha_W[] = ativo ? 0.2 : 0.85
    end

    # mapa de ocupação: occ[estado, walker] = fração do tempo; NaN = célula não usada
    vlines!(ax_E, [0.0]; color=:gray40, linestyle=:dash, linewidth=1.5)   # fim da equilibração

    occ = Observable(fill(NaN32, MAXREP, MAXREP))
    hm = heatmap!(ax_O, occ; colormap=:Blues, colorrange=(0.0, 0.34))
    occ_pos = [Point2f(k, w) for w in 1:MAXREP for k in 1:MAXREP]
    occ_txt = Observable(fill("", MAXREP * MAXREP))
    text!(ax_O, occ_pos; text=occ_txt, align=(:center, :center), fontsize=12, color=:black)
    cb_O = Colorbar(g_O[1, 2], hm; label="fração", width=14)

    pts = [Observable(Point2f[]) for _ in 1:MAXREP]
    cols = [Observable(RGBf(0, 0, 0)) for _ in 1:MAXREP]
    axs = Axis[]

    logtext = Observable("Pronto. Configure e clique em '1. Minimizar + Equilibrar'.")
    Label(fig[2, 1:3], logtext; halign=:left, tellwidth=false)
    rowsize!(fig.layout, 2, Fixed(30))

    # ---------------- estado ----------------
    sim = Ref{Union{Nothing,Sim}}(nothing)
    running = Ref(false)        # pedido de continuar (Stop põe em false)
    busy = Ref(false)           # há uma tarefa de simulação ativa
    lastflag = fill("", MAXREP)
    xrng = MersenneTwister(1234)
    current_colors = Observable(PALETTES["Destaque"])
    Uh = [Float64[] for _ in 1:MAXREP]     # amostras de U por estado (para P(U))
    occ_count = zeros(Int, MAXREP, MAXREP)  # [estado, walker]: nº de registros
    swap_count = zeros(Int, MAXREP, MAXREP) # [walker a, walker b]: trocas aceitas entre os dois
    minU_count = zeros(Int, MAXREP, MAXREP) # [estado, walker]: nº de registros em que a menor U estava ali
    # frames do estado de menor kT durante o REMD: (passo, walker, U, x, y)
    traj = Tuple{Int,Int,Float64,Vector{Float32},Vector{Float32}}[]
    map_mode = Ref("Ocupação")
    Usm = fill(NaN, MAXREP)                 # média móvel exponencial de U por estado (só para o gráfico)
    rt_stage = zeros(Int, MAXREP)          # 0: ainda não tocou o menor kT; 1: a caminho do maior kT; 2: voltando ao menor kT
    rt_count = zeros(Int, MAXREP)          # round trips (menor kT -> maior kT -> menor kT) por walker

    mean_or_nan(v) = isempty(v) ? NaN : sum(v) / length(v)

    function update_hist()
        s = sim[]; R = length(s.reps)
        isempty(Uh[1]) && return
        lo = Inf; hi = -Inf
        for k in 1:R
            isempty(Uh[k]) && continue
            lo = min(lo, minimum(Uh[k])); hi = max(hi, maximum(Uh[k]))
        end
        hi - lo < 1e-9 && return
        for k in 1:R
            hpts[k][] = hist_points(Uh[k], lo, hi, 28)
        end
        autolimits!(ax_H)
        return nothing
    end

    function update_occ()
        s = sim[]; R = length(s.reps)
        m = fill(NaN32, MAXREP, MAXREP); txt = fill("", MAXREP * MAXREP)
        if map_mode[] == "Ocupação"
            # m[estado, walker] = fração do tempo do walker naquele kT (cada linha soma 1)
            for w in 1:R
                tot = max(1, sum(@view occ_count[1:R, w]))
                for k in 1:R
                    f = occ_count[k, w] / tot
                    m[k, w] = f
                    txt[(w - 1) * MAXREP + k] = @sprintf("%.2f", f)
                end
            end
            hm.colorrange = (0.0, 2.0 / R)
            ax_O.xlabel = "índice do kT (1 = menor)"; ax_O.ylabel = "walker"
            ax_O.title = @sprintf("Ocupação: fração do tempo de cada walker em cada kT (ideal = %.2f)", 1 / R)
        elseif map_mode[] == "Menor energia"
            # m[estado, walker] = fração dos registros em que a configuração de MENOR U estava naquele kT e era aquele walker
            tot = max(1, sum(@view minU_count[1:R, 1:R])); mx = 0.0
            for w in 1:R, k in 1:R
                f = minU_count[k, w] / tot
                m[k, w] = f; mx = max(mx, f)
                txt[(w - 1) * MAXREP + k] = @sprintf("%.2f", f)
            end
            pk = [sum(@view minU_count[k, 1:R]) / tot for k in 1:R]
            hm.colorrange = (0.0, max(mx, 1e-3))
            ax_O.xlabel = "índice do kT (1 = menor)"; ax_O.ylabel = "walker"
            ax_O.title = "Onde está a configuração de menor U?\nP por kT: " * join([@sprintf("%.2f", q) for q in pk], "  ")
        elseif map_mode[] == "Taxa de troca (kT)"
            # m[i, j] = P_acc = aceitas/tentadas entre os estados kT_i e kT_j (só vizinhos tentam trocar)
            for i in 1:R-1
                s.att[i] == 0 && continue
                f = s.acc[i] / s.att[i]
                m[i, i+1] = f; m[i+1, i] = f
                txt[i * MAXREP + i] = @sprintf("%.2f", f)             # célula (i, i+1)
                txt[(i - 1) * MAXREP + i + 1] = @sprintf("%.2f", f)   # célula (i+1, i)
            end
            hm.colorrange = (0.0, 1.0)
            ax_O.xlabel = "índice do kT (1 = menor)"; ax_O.ylabel = "índice do kT"
            ax_O.title = "Taxa de troca: P_acc = aceitas/tentadas por par de kT (alvo 0.2-0.3)"
        else
            # m[a, b] = fração de TODAS as trocas aceitas que foram entre os walkers a e b (simétrico)
            tot = max(1, sum(@view swap_count[1:R, 1:R]) ÷ 2)
            npairs = max(1, R * (R - 1) ÷ 2)
            for b in 1:R, a in 1:R
                a == b && continue
                f = swap_count[a, b] / tot
                m[a, b] = f
                txt[(b - 1) * MAXREP + a] = @sprintf("%.3f", f)
            end
            hm.colorrange = (0.0, 2.0 / npairs)
            ax_O.xlabel = "walker"; ax_O.ylabel = "walker"
            ax_O.title = @sprintf("Trocas entre walkers: fração do total de trocas aceitas (ideal = %.3f)", 1 / npairs)
        end
        occ[] = m; occ_txt[] = txt
        return nothing
    end

    function update_plots()
        s = sim[]; s === nothing && return
        R = length(s.reps); pal = current_colors[]
        for k in 1:R
            r = s.reps[k]
            pts[k][] = [Point2f(r.x[i], r.y[i]) for i in eachindex(r.x)]
            cols[k][] = pal[r.walker]
            pstr = (k < R && !isnan(s.last_prob[k])) ? @sprintf(" | P(->)=%.2f", s.last_prob[k]) : ""
            axs[k].title = @sprintf("kT=%.2f | U=%.0f%s\n%s", s.kTs[k], r.U, pstr, lastflag[k])
        end
        kmin = argmin([r.U for r in s.reps])            # painel com a configuração de menor U: moldura magenta
        for k in 1:R
            c = k == kmin ? COR_MENOR_U : RGBf(0, 0, 0); lw = k == kmin ? 4.0 : 1.0
            ax = axs[k]
            ax.leftspinecolor = c; ax.rightspinecolor = c; ax.topspinecolor = c; ax.bottomspinecolor = c
            ax.spinewidth = lw
        end
        notify(minpts)
        pts_main[] = pts[1][]
        col_main[] = pal[s.reps[1].walker]
        ax_main.title = @sprintf("Réplica de referência (menor kT) | kT=%.2f | walker %d | U=%.0f", s.kTs[1], s.reps[1].walker, s.reps[1].U)
        ax_E.title = @sprintf("Energia potencial (média móvel) | <U>: menor kT = %.1f, maior kT = %.1f", mean_or_nan(Uh[1]), mean_or_nan(Uh[R]))
        for w in 1:MAXREP
            if w <= R
                line_cols[w][] = pal[w]
                lw_lines[w][] = (w == 1) ? 3.5 : 1.2
            end
            notify(epts[w]); notify(wpts[w])
        end
        autolimits!(ax_E); autolimits!(ax_W)
        update_occ()
        if sum(s.att) > 0
            lines_acc = [@sprintf("%.2f-%.2f: %.2f (%d/%d)", s.kTs[i], s.kTs[i+1], s.acc[i] / max(1, s.att[i]), s.acc[i], s.att[i]) for i in 1:R-1]
            stats[] = "Taxa de troca acumulada por par de kT:\n" * join(lines_acc, "\n") *
                      "\nRound trips (kT mín-máx-mín): " * string(sum(rt_count[1:R])) *
                      "\n  por walker: " * join(rt_count[1:R], " ")
        end
        return nothing
    end

    function trim!(v, n)
        length(v) > n && deleteat!(v, 1:(length(v) - n))
        return nothing
    end

    # energia (sempre) e amostras para P(U) (só durante o REMD, collect = true)
    function record_energy(collect::Bool)
        s = sim[]
        for k in eachindex(s.reps)
            U = s.reps[k].U
            Usm[k] = isnan(Usm[k]) ? U : 0.9 * Usm[k] + 0.1 * U     # janela de ~10 registros (100 passos no Run; 200 na equilibração)
            v = epts[k][]
            push!(v, Point2f(s.step, Usm[k])); trim!(v, 4000)
            if collect
                push!(Uh[k], s.reps[k].U); trim!(Uh[k], 20000)
            end
        end
        if collect
            kmin = argmin([r.U for r in s.reps])
            minU_count[kmin, s.reps[kmin].walker] += 1
            r = s.reps[1]                                   # réplica no menor kT (a cada 10 passos)
            length(traj) < 50000 && push!(traj, (s.step, r.walker, r.U, Float32.(r.x), Float32.(r.y)))
        end
        return nothing
    end

    function record_walk()
        s = sim[]; R = length(s.reps)
        push!(minpts[], Point2f(s.nattempt, argmin([r.U for r in s.reps])))
        for k in 1:R
            w = s.reps[k].walker
            push!(wpts[w][], Point2f(s.nattempt, k))
            occ_count[k, w] += 1
            if k == 1
                rt_stage[w] == 2 && (rt_count[w] += 1)
                rt_stage[w] = 1
            elseif k == R && rt_stage[w] == 1
                rt_stage[w] = 2
            end
        end
        return nothing
    end

    function clear_series!()
        for k in 1:MAXREP
            pts[k][] = Point2f[]; epts[k][] = Point2f[]; wpts[k][] = Point2f[]; hpts[k][] = Point2f[]
            empty!(Uh[k])
        end
        minpts[] = Point2f[]; empty!(traj); fill!(occ_count, 0); fill!(swap_count, 0); fill!(minU_count, 0); fill!(rt_stage, 0); fill!(rt_count, 0); fill!(lastflag, ""); fill!(Usm, NaN)
        return nothing
    end

    function reset!()
        v = [parse(Float64, strip(a)) for a in split(tbtext(tb_kts), ',') if !isempty(strip(a))]
        sort!(v); R = length(v)
        R < 2 && error("informe ao menos 2 temperaturas")
        R > MAXREP && error("máximo de $MAXREP réplicas")
        L = getf(tb_L); seed = geti(tb_seed)
        seed <= 0 && (seed = rand(1:10^6))
        sim[] = new_sim(v, geti(tb_n), L, getf(tb_eps), getf(tb_sig), getf(tb_dt), getf(tb_gam); seed)
        Random.seed!(xrng, seed + 99)
        clear_series!()
        if length(axs) != R
            for ax in axs
                delete!(ax)
            end
            empty!(axs)
            for k in 1:R
                ax = Axis(right_grid[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1]; aspect=DataAspect(), titlesize=12)
                hidedecorations!(ax)
                scatter!(ax, pts[k]; color=cols[k], markersize=6, strokewidth=0.5, strokecolor=:black)
                push!(axs, ax)
            end
        end
        for k in 1:R
            axs[k].title = ""; axs[k].titlecolor = :black
            limits!(axs[k], 0, L, 0, L)
        end
        limits!(ax_main, 0, L, 0, L)
        limits!(ax_O, 0.5, R + 0.5, 0.5, R + 0.5)
        ax_O.xticks = 1:R; ax_O.yticks = 1:R
        ax_W.yticks = 1:R                          # índices inteiros no random walk
        hm.colorrange = (0.0, 2.0 / R)            # ideal (1/R) no meio da escala
        for k in 1:MAXREP
            tcols[k][] = tempcolor(min(k, R), R)
            lw_e[k][] = (k == 1 || k == R) ? 3.0 : 1.3      # extremos (menor e maior kT) em destaque
        end
        stats[] = "Taxa de troca: —"
        record_walk()
        update_plots()
        logtext[] = "Reset: $R réplicas criadas (semente $seed). Clique em '1. Minimizar + Equilibrar'."
        return nothing
    end

    function equilibrate()
        s = sim[]
        logtext[] = "Minimizando..."; yield()
        prepare!(s)
        neq = geti(tb_nequil)
        logtext[] = "Equilibrando cada réplica no seu kT (sem trocas)..."
        running[] = true
        for it in 1:neq
            running[] || break
            step_all!(s)
            if it % 20 == 0
                record_energy(false); update_plots(); sleep(0.001)
            end
        end
        stopped = !running[]
        running[] = false
        for k in 1:MAXREP                       # equilibração fica em passos negativos; produção começa em 0
            epts[k][] = [Point2f(q[1] - s.step, q[2]) for q in epts[k][]]
        end
        s.step = 0
        minpts[] = Point2f[]; empty!(traj); fill!(occ_count, 0); fill!(swap_count, 0); fill!(minU_count, 0); fill!(rt_stage, 0); fill!(rt_count, 0)
        for k in 1:MAXREP
            wpts[k][] = Point2f[]
        end
        record_walk()
        update_plots()
        logtext[] = stopped ? "Equilibração interrompida (Stop). O sistema pode não estar equilibrado." :
                              "Equilibração concluída. Clique em '2. Run'."
        return nothing
    end

    function run_remd()
        s = sim[]
        s.prepared || error("minimize e equilibre primeiro (botão 1)")
        nsteps = geti(tb_nsteps); nex = max(1, geti(tb_nex)); ptime = getf(tb_pause)
        s.dt = getf(tb_dt); s.gam = getf(tb_gam)
        running[] = true
        for it in 1:nsteps
            running[] || break
            step_all!(s)
            if s.step % nex == 0
                res = attempt_exchange!(s, xrng)
                fill!(lastflag, "")
                for k in eachindex(s.reps)
                    axs[k].titlecolor = :black
                end
                swapped = false
                for (i, ok, P) in res
                    lastflag[i] = ok ? "Trocou" : "Falhou"; lastflag[i+1] = lastflag[i]
                    c = ok ? RGBf(0.1, 0.6, 0.1) : RGBf(0.8, 0.2, 0.2)
                    axs[i].titlecolor = c; axs[i+1].titlecolor = c
                    swapped |= ok
                    if ok
                        wa = s.reps[i].walker; wb = s.reps[i+1].walker
                        swap_count[wa, wb] += 1; swap_count[wb, wa] += 1
                    end
                end
                record_walk()
                update_plots(); update_hist()
                if swapped && ptime > 0
                    sleep(ptime)                 # pausa só quando houve troca aceita
                end
            end
            if s.step % 10 == 0
                record_energy(true); update_plots(); sleep(0.001)
                logtext[] = "Rodando REMD: passo $(s.step) ($it/$nsteps)"
            end
        end
        running[] = false
        update_plots(); update_hist()
        logtext[] = "Parado no passo $(s.step). '2. Run' continua daqui; Reset recomeça."
        return nothing
    end

    # Executa f em tarefa assíncrona: uma por vez, com erro mostrado no rodapé
    function launch(f)
        if busy[]
            logtext[] = "Já há uma simulação em andamento. Use Stop antes."
            return nothing
        end
        busy[] = true
        @async try
            f()
        catch err
            logtext[] = "ERRO: " * sprint(showerror, err)
        finally
            running[] = false; busy[] = false
        end
        return nothing
    end

    # carrega os valores de um teste do tutorial nos campos (não roda nada)
    on(menu_teste.selection) do sel
        (sel === nothing || sel == "(escolher)") && return
        v = merge(PADRAO, Dict(TESTES)[sel])
        for (k, tb) in (("nrep", tb_nrep), ("kmin", tb_kmin), ("kmax", tb_kmax), ("n", tb_n), ("L", tb_L),
                        ("dt", tb_dt), ("nsteps", tb_nsteps), ("nex", tb_nex), ("pause", tb_pause),
                        ("nequil", tb_nequil), ("gam", tb_gam), ("eps", tb_eps), ("sig", tb_sig), ("seed", tb_seed))
            settext!(tb, v[k])
        end
        if haskey(v, "kts")
            settext!(tb_kts, v["kts"])
        else
            R = parse(Int, v["nrep"]); a = parse(Float64, v["kmin"]); b = parse(Float64, v["kmax"])
            settext!(tb_kts, join([@sprintf("%.3f", a * (b / a)^((i - 1) / (R - 1))) for i in 1:R], ","))
        end
        logtext[] = "Teste $sel carregado (semente 42). Clique em '1. Minimizar + Equilibrar' e depois em '2. Run'."
    end
    on(menu_cores.selection) do pal
        pal === nothing && return
        current_colors[] = PALETTES[pal]
        update_plots()
    end
    on(menu_map.selection) do sel
        sel === nothing && return
        map_mode[] = sel
        sim[] === nothing || update_occ()
    end
    on(btn_ladder.clicks) do _
        try
            R = clamp(geti(tb_nrep), 2, MAXREP); a = getf(tb_kmin); b = getf(tb_kmax)
            kTs = [a * (b / a)^((i - 1) / (R - 1)) for i in 1:R]
            settext!(tb_kts, join([@sprintf("%.3f", k) for k in kTs], ","))
            logtext[] = "Escada gerada. Clique em Reset (ou no botão 1) para aplicar."
        catch err
            logtext[] = "ERRO: " * sprint(showerror, err)
        end
    end
    on(btn_reset.clicks) do _
        if busy[]
            running[] = false
            logtext[] = "Simulação parada. Clique em Reset de novo para recriar o sistema."
        else
            try
                reset!()
            catch err
                logtext[] = "ERRO: " * sprint(showerror, err)
            end
        end
    end
    on(btn_min.clicks) do _
        launch() do
            reset!(); equilibrate()
        end
    end
    on(btn_run.clicks) do _
        launch(run_remd)
    end
    # recorta da imagem da janela a região ocupada pelos blocos 'objs' (eixo + título + rótulos) e salva
    function crop_save(fn, img, objs; pad=10)
        H, W = size(img); sc = W / size(fig.scene)[1]          # pixels da imagem por unidade da figura
        x0 = Inf; x1 = -Inf; y0 = Inf; y1 = -Inf
        for o in objs
            bb = o.layoutobservables.computedbbox[]; pr = o.layoutobservables.protrusions[]
            x0 = min(x0, bb.origin[1] - pr.left - pad);   x1 = max(x1, bb.origin[1] + bb.widths[1] + pr.right + pad)
            y0 = min(y0, bb.origin[2] - pr.bottom - pad); y1 = max(y1, bb.origin[2] + bb.widths[2] + pr.top + pad)
        end
        cols = clamp(floor(Int, x0 * sc) + 1, 1, W):clamp(ceil(Int, x1 * sc), 1, W)
        rows = clamp(H - ceil(Int, y1 * sc) + 1, 1, H):clamp(H - floor(Int, y0 * sc), 1, H)   # linha 1 = topo
        save(fn, img[rows, cols])
        return nothing
    end

    on(btn_save.clicks) do _
        try
            s = sim[]; R = length(s.reps)
            dir = @sprintf("remd_figs_R%d_kT%.2f-%.2f_n%d_troca%d_passo%d", R, s.kTs[1], s.kTs[end],
                           length(s.reps[1].x), geti(tb_nex), s.step)
            mkpath(dir)
            img = GLMakie.Makie.colorbuffer(fig)
            save(joinpath(dir, "janela_completa.png"), img)
            crop_save(joinpath(dir, "random_walk.png"), img, [ax_W])
            crop_save(joinpath(dir, "energia.png"), img, [ax_E])
            crop_save(joinpath(dir, "P_U.png"), img, [ax_H])
            crop_save(joinpath(dir, "replica_referencia.png"), img, [ax_main])
            crop_save(joinpath(dir, "replicas.png"), img, axs)
            old = map_mode[]                                   # um PNG para cada opção do seletor "Mapa"
            for (mode, tag) in (("Ocupação", "ocupacao"), ("Menor energia", "menor_energia"),
                                ("Taxa de troca (kT)", "taxa_de_troca"), ("Trocas entre walkers", "trocas_walkers"))
                map_mode[] = mode; update_occ()
                crop_save(joinpath(dir, "mapa_" * tag * ".png"), GLMakie.Makie.colorbuffer(fig), [ax_O, cb_O])
            end
            map_mode[] = old; update_occ()
            logtext[] = "Figuras salvas em: " * joinpath(pwd(), dir)
        catch err
            logtext[] = "ERRO ao salvar: " * sprint(showerror, err)
        end
    end
    on(btn_traj.clicks) do _
        try
            isempty(traj) && error("sem frames: rode o REMD primeiro (botão 2)")
            s = sim[]; R = length(s.reps); n = length(traj[1][4])
            fn = @sprintf("remd_traj_kT%.2f_R%d_n%d_%dframes.xyz", s.kTs[1], R, n, length(traj))
            open(fn, "w") do io
                for (st, w, U, x, y) in traj
                    println(io, n)
                    @printf(io, "passo=%d walker=%d U=%.4f kT=%.3f L=%.2f\n", st, w, U, s.kTs[1], s.L)
                    for i in 1:n
                        @printf(io, "H %.4f %.4f 0.0000\n", x[i], y[i])   # "H": o VMD nunca liga H com H
                    end
                end
            end
            # script do VMD: sem ligações automáticas, esferas de diâmetro ~sig, caixa periódica, vista ortográfica
            tcl = replace(fn, r"\.xyz$" => ".tcl")
            open(tcl, "w") do io
                println(io, "mol new {$fn} type xyz autobonds off waitfor all")
                println(io, "mol delrep 0 top")
                println(io, "mol representation VDW $(round(s.sig / 2 / 1.0, digits=3)) 20")   # raio = sig/2 (H: 1.0 A)
                println(io, "mol color ColorID 1")
                println(io, "mol addrep top")
                println(io, "pbc set {$(s.L) $(s.L) 10.0} -all")
                println(io, "pbc box -center origin -shiftcenter {$(s.L / 2) $(s.L / 2) 0}")
                println(io, "display projection Orthographic")
                println(io, "display resetview")
            end
            logtext[] = "Trajetória salva (.xyz + .tcl, $(length(traj)) frames) em $(pwd()). Abrir com: vmd -e " * tcl
        catch err
            logtext[] = "ERRO ao salvar: " * sprint(showerror, err)
        end
    end
    on(btn_stop.clicks) do _
        running[] = false
    end

    reset!()
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    fig = remd_gui()
    wait(display(fig))
end
