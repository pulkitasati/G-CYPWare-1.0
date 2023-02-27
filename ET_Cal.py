import math
import sys

lamda = float(sys.argv[1])
delta_G  = float(sys.argv[2])
R = float(sys.argv[3])
file_name = sys.argv[4]
print(f'{file_name}.txt') 


A = 0.022  # pre-exponential factor (eV)
b = 1.39  
Ra = 3.6  # 

Hab = A * math.exp(-b * (R - Ra) / 2)
hab2= Hab * Hab

print("Hab = {:.10f} eV".format(Hab)) 
print("hab2 = {:.10f} eV".format(hab2))  # output Hab value with 4 decimal places


p = 6.582119569e-16   # momentum in eV s / m
T = 300               # temperature in K
kb = 8.617333262e-5   # Boltzmann constant in eV / K
term1 = 2 * math.pi / p
term2 = hab2
term3 = 1 / math.sqrt(4 * math.pi * lamda * kb * T)
term4 = math.exp(-(delta_G + lamda) ** 2 / (4 * lamda * kb * T))
K = term1 * term2 * term3 * term4
num = float(K)
log_num = math.log(num)

print(f"lamda = {lamda}")
print(f"delta_G = {delta_G}")
print(f"term1 = {term1}")
print(f"term2 = {term2}")
print(f"term3 = {term3}")
print(f"term4 = {term4}")
print(f"K = {K}")
print("loge({}) = {}".format(num, log_num))

with open(f'{file_name}.txt', 'w', encoding='utf-8') as f:

    # Write the variable names and values to the file
    f.write(f"lamda = {lamda}\n")
    f.write(f"delta_G = {delta_G}\n")
    f.write(f"term1 = {term1}\n")
    f.write(f"term2 = {term2}\n")
    f.write(f"term3 = {term3}\n")
    f.write(f"term4 = {term4}\n")
    f.write(f"K = {K}\n")
    f.write(f"log_num = {log_num}\n")
