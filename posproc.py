import numpy as np
import matplotlib.pyplot as plt

def compute_analytical_solution(n, L, rho, E, nu):
    """
    Computes the analytical resonance frequencies for the 1D P-wave propagation problem.

    Args:
        n (int or array-like): The mode number(s) (0, 1, 2, ...).
        L (float): Length of the domain.
        rho (float): Density of the material.
        E (float): Elastic modulus.
        nu (float): Poisson's ratio.

    Returns:
        omega (float or array-like): The analytical resonance frequency(ies) in s^-1.
    """
    # Constrained modulus for P-wave in 1D (plane strain assumption usually, but text uses c1)
    # Assuming 1D bar velocity c1 = sqrt(E / rho) based on context, or constrained P-wave velocity.
    # Let's use the constrained P-wave velocity c1 = sqrt((E * (1 - nu)) / (rho * (1 + nu) * (1 - 2*nu)))
    # For a simple 1D rod, c = sqrt(E/rho). We'll use the constrained one which is more common in these plane problems.

    M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    c1 = np.sqrt(M / rho)

    omega = (2 * np.array(n) + 1) * np.pi * c1 / (2 * L)
    return omega

# Example usage to verify the first frequency:
# Using parameters from the text: L=6, rho=100 (from text, though input_data has 7850), E needs to match mu=10^6
# Let's recalculate based on text to see if it matches 45.345
# mu = E / (2*(1+nu)) => E = mu * 2 * (1+nu) = 10^6 * 2 * 1.25 = 2.5 * 10^6
# L = 6
# rho = 100
# nu = 0.25

L_text = 6.0
rho_text = 100.0
nu_text = 0.25
mu_text = 1.0e6
E_text = mu_text * 2 * (1 + nu_text)

w_n = compute_analytical_solution(np.array([0, 1, 2]), L_text, rho_text, E_text, nu_text)
print(f"Analytical frequencies (n=0,1,2): {w_n}")

import numpy as np
import matplotlib.pyplot as plt

# 1. Parâmetros do problema
L_text = 6.0
rho_text = 100.0
nu_text = 0.25
mu_text = 1.0e6
E_text = mu_text * 2 * (1 + nu_text)
p = 100.0  # Tração aplicada

# Módulo de deformação confinada e velocidade da onda P
M = E_text * (1 - nu_text) / ((1 + nu_text) * (1 - 2 * nu_text))
c1 = np.sqrt(M / rho_text)

# Frequência selecionada para comparação (combina com os dados numéricos)
w_eval = 70.0

# 2. Solução Analítica (Deslocamento 1D)
def analytical_displacement(y, w, p, c1, M, L):
    # u(y) = [p * c1 / (M * w * cos(w * L / c1))] * sin(w * y / c1)
    return (p * c1 / (M * w * np.cos(w * L / c1))) * np.sin(w * y / c1)


# Gerando pontos para a curva analítica
y_analyt = np.linspace(0, L_text, 100)
u_analyt = analytical_displacement(y_analyt, w_eval, p, c1, M, L_text)

# 3. Dados Numéricos (Extraídos anteriormente de squa4.out para w=70)
y_num = np.array([1., 2., 3., 4., 5.])
u_num_real = np.array([-0.4265E-04, -0.7852E-04, -0.1018E-03, -0.1085E-03, -0.9723E-04])

# 4. Plot da Comparação
plt.figure(figsize=(10, 6))
plt.plot(y_analyt, u_analyt * 1e4, 'r-', linewidth=2, \
         label=f'Solução Analítica ($\\omega$ = {w_eval} rad/s)')
plt.plot(y_num, u_num_real * 1e4, 'bs', markersize=8, \
         label='Resultados Numéricos (QUADPLEH)')

plt.xlabel('Posição $y$ (m)')
plt.ylabel('Parte Real do Deslocamento Vertical $u_y \\times 10^4$ (m)')
plt.title('Comparação: Solução Analítica vs Numérica ao longo do eixo central')
plt.legend()
plt.grid(True)
plt.xlim(0, L_text)
plt.savefig("result_plot.png")

frequencies = np.linspace(10.,150.,30)

# 1. Carrega os dados do arquivo disp.dat gerado pelo solver Fortran
freq_num = []
uy_num_mod = []

try:
    with open('disp_trac.dat', 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

        for i, line in enumerate(lines):
            if i >= len(frequencies):
                break

            parts = line.split()
            if len(parts) < 2:
                continue

            # A estrutura gerada é: (ux_re,ux_im) (uy_re,uy_im)
            # O deslocamento 'y' é o segundo item (índice 1)
            uy_str = parts[1].strip('()')

            re_uy, im_uy = uy_str.split(',')
            uy_mod = np.abs(complex(float(re_uy), float(im_uy)))

            freq_num.append(frequencies[i])
            uy_num_mod.append(uy_mod)

    freq_num = np.array(freq_num)
    uy_num_mod = np.array(uy_num_mod)
    print(f"Foram lidos {len(freq_num)} registros de disp_trac.dat.")

    # 2. Parâmetros analíticos
    L_text = 6.0
    rho_text = 100.0
    nu_text = 0.25
    mu_text = 1.0e6
    E_text = mu_text * 2 * (1 + nu_text)
    p = 100.0

    M = E_text * (1 - nu_text) / ((1 + nu_text) * (1 - 2 * nu_text))
    c1 = np.sqrt(M / rho_text)
    w1 = np.pi * c1 / (2 * L_text) # 1ª Frequência natural

    # 3. Curva Analítica contínua para u(L)
    w_analyt = np.linspace(0.1, np.max(freq_num) + 10, 1000)
    u_analyt = (p * c1 / (M * w_analyt)) * np.tan(w_analyt * L_text / c1)
    u_analyt_mod = np.abs(u_analyt)

    # 4. Plot da Comparação
    plt.figure(figsize=(10, 6))
    plt.plot(w_analyt / w1, u_analyt_mod * 1e4, 'r-', label='Solução Analítica (Contínua)')
    plt.plot(freq_num / w1, uy_num_mod * 1e4, 'bs', markersize=8, label='Resultados Numéricos (disp.dat)')

    # Linha vertical para a primeira ressonância
    plt.axvline(1.0, color='gray', linestyle='--', alpha=0.5, label='1ª Freq. Ressonância')

    plt.xlabel('Frequência Normalizada ($\\omega / \\omega_1$)')
    plt.ylabel('Módulo do Deslocamento Vertical $|u_y| \\times 10^4$ (m)')
    plt.title("Comparação FRF: Analítico vs Numérico (Extraído de disp.dat)")

    # Limitando o eixo y para focar nos resultados sem que as assíntotas estraguem a escala
    if len(uy_num_mod) > 0:
        plt.ylim(0, np.max(uy_num_mod * 1e4) + 2)
    plt.xlim(0, np.max(w_analyt / w1))
    plt.legend()
    plt.grid(True)
    plt.savefig("result_plot.png")

except Exception as e:
    print(f"Erro ao processar disp.dat: {e}")
