from math import pi, sin

import matplotlib.pyplot as plt
# Newmark methods ;Non-linear systems
# You will have to input values for stiffness K

# initializations
# Average acceleration method
gamma = 0.5
beta = 0.25
ultimateResidualForce = 0.001
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
# initial state determination
kT = K
forcingFunction = K
fs = 0

yieldStrength = 2.0
acc = (fs - c * vel - K * disp) / mass

delta_t = 0.01
a_1 = (mass / (beta * delta_t ** 2)) + ((gamma * c) / (beta * delta_t))
a_2 = (mass / (beta * delta_t)) + ((gamma / beta) - 1) * c
a_3 = ((1 / (2 * beta)) - 1) * mass + ((gamma / (2 * beta)) - 1) * (c * delta_t)
print(a_1, a_2, a_3)

# bucket initializations
displacement, velocity, acceleration, timer = [], [], [], []
kism = []
# initial append
displacement.append(disp)
velocity.append(vel)
acceleration.append(acc)
timer.append(0)
kism.append(fs)
t = delta_t
while t <= 2 * t_d:
    disp_plus = disp
    fs_plus = fs
    kT_plus = kT
    p = 0
    if t <= t_d:
        p = forcingFunction * sin(omega * t)
    if t >= t_d:
        p = 0
    p_prime = p + (a_1 * disp) + (a_2 * vel) + (a_3 * acc)
    residualForce = p_prime - fs_plus - a_1 * disp_plus
    kism.append(fs_plus)
    while abs(residualForce) > ultimateResidualForce:
        kT_plusPrime = kT_plus + a_1
        delta_disp = residualForce / kT_plusPrime
        disp_plus = disp_plus + delta_disp
        fs_plus = fs + K * (disp_plus - disp)

        if fs_plus > yieldStrength:
            fs_plus = yieldStrength
            kT_plus = 0
        else:
            kT_plus = K
        residualForce = p_prime - fs_plus - a_1 * disp_plus

    vel_plus = (2 * (disp_plus - disp) / delta_t) - vel
    acc = (4 * (disp_plus - disp) / delta_t ** 2) - (4 * vel / delta_t) - acc
    timer.append(t)
    t = t + delta_t
    displacement.append(disp_plus)
    velocity.append(vel_plus)
    acceleration.append(acc)
    disp = disp_plus
    vel = vel_plus
    fs = fs_plus
    kT = kT_plus
# for i,j,k in zip(displacement,velocity,acceleration):
#   print(i ,"!" ,j,"!",k)
plt.subplot(2, 1, 1)
plt.title("Displacement response for yielded SDOF")
plt.plot(timer, displacement, c='b')
plt.subplot(2, 1, 2)
plt.plot(displacement, kism)
plt.title("Force_deformation relation")
plt.subplot_tool()  # adjust subplot
plt.show()
