from math import pi, sin
import matplotlib.pyplot as plt

# Newmark methods ;linear systems

# initializations
# Constant average acceleration method


gamma = 0.5
beta = 0.25
# initial displacement and velocity are zero
initials = [0, 0]
disp, vel = initials
K = int(input("Give value for initial stiffness(kN/m):  "))

T_n = 1
# td/Tn = 0.75
t_d = 0.75 * 1
omega = 2 * pi / t_d
mass = K / omega ** 2

# undamped SDOF => c = 0
c = 0

# initial calculations

forcingFunction = 10
acc = ((forcingFunction * sin(0)) - c * vel - K * disp) / mass

delta_t = 0.01
a_1 = (mass / (beta * delta_t ** 2)) + ((gamma * c) / (beta * delta_t))
a_2 = (mass / (beta * delta_t)) + ((gamma / beta) - 1) * c
a_3 = ((1 / (2 * beta)) - 1) * mass + ((gamma / (2 * beta)) - 1) * (c * delta_t)
print(a_1, a_2, a_3)
k_prime = K + a_1

# bucket initializations of velocity, acc, and displacement
displacement, velocity, acceleration, timer = [], [], [], []
# initial append
displacement.append(disp)
velocity.append(vel)
acceleration.append(acc)
timer.append(0)
t = delta_t
while t <= 2 * t_d:
    p = 0
    if t <= t_d:
        p = forcingFunction * sin(omega * t)
    if t >= t_d:
        p = 0
    p_prime = p + (a_1 * disp) + (a_2 * vel) + (a_3 * acc)
    disp_new = p_prime / k_prime
    vel_new = (2 * (disp_new - disp) / delta_t) - vel
    acc = (4 * (disp_new - disp) / delta_t ** 2) - (4 * vel / delta_t) - acc
    timer.append(t)
    t = t + delta_t
    displacement.append(disp_new)
    velocity.append(vel_new)
    acceleration.append(acc)
    disp = disp_new
    vel = vel_new

# for i,j,k in zip(displacement,velocity,acceleration):
#	print(i ,"!" ,j,"!",k)

plt.plot(timer, displacement, c='r')
plt.title("Undamped SDOF system subjected to full cycle sine forcing function")
plt.show()
