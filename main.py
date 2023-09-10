import numpy as np
import matplotlib.pyplot as plt

pi = np.pi

# Dados coletados experimentalmente
y0 = 0.113
y3 = 0.05843
delta = 0.2198
E = 2.78 * 1000
I = 2.08
L_viga = 0.25
L_total = 0.30
m = 0.012
A = y0

personalizar_dados = False

# Entradas pelo usuário
m2 = float(
    input("Insira a massa adicionada a extremidade da viga (kg): ")) + m

duracao = int(
    input("Insira o tempo desejado para a simulação (s): "))


# Validação de alteração dos dados base do problema
while (personalizar_dados != 's' and personalizar_dados != "n"):
    personalizar_dados = input(
        "Deseja inserir seus próprios dados para as demais variáveis do problema? (s/n): ")

    if (personalizar_dados != 's' and personalizar_dados != "n"):
        print("Responda com 's' para sim ou 'n' para não. \n")


# Alteração de dados base do problema
if (personalizar_dados == 's'):
    y0 = float(
        input("Insira a posição em t=0 (m): "))

    y3 = float(
        input("Insira a posição do terceiro pico (m): "))

    E = float(
        input("Insira o módulo de elasticidade (GPa): ")) * 1000

    I = float(
        input("Insira o momento de inércia (mm^4): "))

    L_total = float(
        input("Insira o comprimento total (m): "))

    L_viga = float(
        input("Insira o comprimento livre da viga (m): "))

    m = float(
        input("Insira a massa total da viga (kg): "))


def calcular_decremento_logaritmico():
    return np.log(y0 / y3) / 3


def calcular_razao_amortecimento(delta):
    return delta / np.sqrt((2*pi)**2 + delta**2)


def calcular_rigidez_equivalente():
    return ((3 * E * I) / ((L_viga * 1000) ** 3)) * 1000


def calcular_massa_equivalente(m):
    m_viga = (m * L_viga) / L_total
    return (33/140) * m_viga


def calcular_frequencia_natural(k_eq, m_eq):
    omega_n = np.sqrt(k_eq / m_eq)
    return omega_n


def calcular_frequencia_amortecida(zeta, omega_n):
    return omega_n * np.sqrt((1 - zeta**2))


def calcular_periodo_amortecido(omega_d):
    return (2 * pi) / omega_d


def sistema_amortecido_subcritico(t, A, zeta, omega_n, omega_d, phi=0):
    return A * np.exp(-zeta * omega_n * t) * np.cos(omega_d * t + phi)


# Cálculo das variáveis para a descrição da posição sem massa adicional
delta = calcular_decremento_logaritmico()
zeta = calcular_razao_amortecimento(delta)
k_eq = calcular_rigidez_equivalente()
m_eq = calcular_massa_equivalente(m)
omega_n = calcular_frequencia_natural(k_eq, m_eq)
omega_d = calcular_frequencia_amortecida(zeta, omega_n)
tau = calcular_periodo_amortecido(omega_d)

# Cálculo das variáveis para a descrição da posição com massa adicional
m_eq2 = calcular_massa_equivalente(m2)
omega_n2 = calcular_frequencia_natural(k_eq, m_eq2)
omega_d2 = calcular_frequencia_amortecida(zeta, omega_n2)
tau2 = calcular_periodo_amortecido(omega_d2)

# Criação do vetor de tempo e cálculo dos deslocamentos
tempo = np.linspace(0, duracao, 1000)

deslocamento_sem_massa = sistema_amortecido_subcritico(
    tempo, A, zeta, omega_n, omega_d)
deslocamento_massa_adicional = sistema_amortecido_subcritico(
    tempo, A, zeta, omega_n2, omega_d2)

# Exibição dos valores das variáveis calculadas
print(
    f"\n \n VARIÁVEIS CALCULADAS: \n delta = {delta} \n zeta = {zeta} \n k_eq = {k_eq}")
print(
    f"\n SEM MASSA ADICIONAL: \n m_eq = {m_eq} \n omega_n = {omega_n} \n omega_d = {omega_d} \n tau = {tau}")
print(
    f"\n COM MASSA ADICIONAL: \n m_eq = {m_eq2} \n omega_n = {omega_n2} \n omega_d = {omega_d2} \n tau = {tau2} \n")

# Plotagem do gráfico
plt.figure(figsize=(10, 6))
plt.plot(tempo, deslocamento_sem_massa,
         label='Deslocamento sem massa adicional')
plt.plot(tempo, deslocamento_massa_adicional,
         label='Deslocamento com massa adicional')
plt.xlabel('Tempo (s)')
plt.ylabel('Deslocamento (m)')
plt.title('Vibração de viga engastada')
plt.grid(True)
plt.legend()
plt.show()
