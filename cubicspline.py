# TFY4106/4107/4125 Fysikk våren 2020.
#
# Programmet bestemmer høyden til de 8 festepunktene ved å trekke
# tilfeldige heltall mellom 50 og 300 (mm). Når festepunktene er
# akseptable, beregnes baneformen med mm som enhet. 
# Etter dette regnes høydeverdiene om til meter som enhet.
# Hele banens form y(x) beregnes
# ved hjelp av 7 ulike tredjegradspolynomer, på en slik måte
# at både banen y, dens stigningstall dy/dx og dens andrederiverte
# d2y/dx2 er kontinuerlige i alle 6 indre festepunkter.
# I tillegg velges null krumning (andrederivert) 
# i banens to ytterste festepunkter (med bc_type='natural' nedenfor).
# Dette gir i alt 28 ligninger som fastlegger de 28 koeffisientene
# i de i alt 7 tredjegradspolynomene.

# Programmet aksepterer de 8 festepunktene først når
# følgende betingelser er oppfylt:
#
# 1. Starthøyden er stor nok og en del høyere enn i de øvrige 7 festepunktene.
# 2. Helningsvinkelen i startposisjonen er ikke for liten.
# 3. Banens maksimale helningsvinkel er ikke for stor.
#
# Med disse betingelsene oppfylt vil 
# (1) objektet (kula/skiva/ringen) fullføre hele banen selv om det taper noe 
#     mekanisk energi underveis;
# (2) objektet få en fin start, uten å bruke for lang tid i nærheten av
#     startposisjonen; 
# (3) objektet forhåpentlig rulle rent, uten å gli/slure.

# Vi importerer nødvendige biblioteker:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.constants import g as GRAVITY, pi as PI
import math
from velocity import csv_to_list


# Horisontal avstand mellom festepunktene er 200 mm
h = 200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])
# Vi begrenser starthøyden (og samtidig den maksimale høyden) til
# å ligge mellom 250 og 300 mm
# yfast: tabell med 8 heltall mellom 50 og 300 (mm); representerer
# høyden i de 8 festepunktene
# inttan: tabell med 7 verdier for (yfast[n+1]-yfast[n])/h (n=0..7); dvs
# banens stigningstall beregnet med utgangspunkt i de 8 festepunktene.
# while-løkken sjekker om en eller flere av de 3 betingelsene ovenfor
# ikke er tilfredsstilt; i så fall velges nye festepunkter inntil
# de 3 betingelsene er oppfylt


# Når programmet her har avsluttet while-løkka, betyr det at
# tallverdiene i tabellen yfast vil resultere i en tilfredsstillende bane. 

# Omregning fra mm til m:
xfast = xfast/1000
yfast = np.asarray([0.297, 0.249, 0.235, 0.161, 0.175, 0.122, 0.193, 0.159])
# Programmet beregner deretter de 7 tredjegradspolynomene, et
# for hvert intervall mellom to nabofestepunkter.
EXPERIMENTAL = False
BALL = False
data = csv_to_list(csv="ballxy" if BALL else "circlexy", velocity=False)

x_experimental = np.asarray([data[i][1] for i in range(1, len(data))])[2:]
y_experimental = np.asarray([data[i][2] for i in range(1, len(data))])[2:]
time_data = np.asarray([data[i][0] for i in range(1, len(data))])[2:]

# Med scipy.interpolate-funksjonen CubicSpline:
time_data -= 0.15
c = 2/5
MASS_BALL = 30 / 1000 # kilograms
MASS_FILLED_CIRCLE = 27 / 1000# kilograms
MASS_EMPTY_CIRCLE = 13 / 1000 # kilograms

MASS = MASS_BALL if BALL else MASS_EMPTY_CIRCLE


def derivatives(experimental=False):
    cs = CubicSpline(xfast, yfast, bc_type='natural')
    xmin = 0.000
    xmax = 1.401
    data_points = len(x_experimental) if experimental else 1401
    dx = xmax / data_points

    if experimental:
        x = x_experimental
        Nx = len(x)
        y = y_experimental
    else:
        x = np.arange(xmin, xmax, dx)
        Nx = len(x)
        y = cs(x)
    dy = cs(x, 1) if not EXPERIMENTAL or BALL else cs(x, 1)[:-1]
    d2y = cs(x, 2) if not EXPERIMENTAL or BALL else cs(x, 2)[:-1]
    return x, Nx, y, dy, d2y


x_e, Nx_e, y_e, dy_e, d2y_e = derivatives(True)
x_n, Nx_n, y_n, dy_n, d2y_n = derivatives(False)


# Plotting
def track(x_data, y_data, x_fast, y_fast):
    track_plot = plt.figure('y(x)',figsize=(12,3))
    plt.plot(x_data, y_data, x_fast, y_fast, '*')
    plt.title('Banens form')
    plt.xlabel('$x$ (m)',fontsize=20)
    plt.ylabel('$y(x)$ (m)',fontsize=20)
    plt.ylim(0,0.350)
    plt.grid()
    plt.show()
    track_plot.savefig("track_plot.pdf", bbox_inches='tight')


# Helningsvinkel
def beta(dy):
    return np.asarray([math.atan(derivative)/(2 * PI) * 180 for derivative in dy])


# BETA FØRSTE PLOT
def beta_plot():
    y_numeric, y_experimental = beta(dy_n), beta(dy_e)
    x_label = r'$x$ (m)'
    y_label = r'$\beta$ (grader)'
    name = 'beta_plot'
    plot_data(x_n, y_numeric, x_e, y_experimental, x_label, y_label, name)
# END BETA FØRSTE PLOT


# START KRUMNINGSPLOT
def curve(d2y, dy):
    return d2y/(1 + dy ** 2) ** (3/2)


def curve_plot():
    y_numeric, y_experimental = curve(d2y_n, dy_n), curve(d2y_e, dy_e)
    x_label = r'$x$(m)'
    y_label = r'$\kappa$(x) (1/m)'
    name = "curve_plot"

    plot_data(x_n, y_numeric, x_e, y_experimental, x_label, y_label, name)
# END KRUMNINGSPLOTs


def v(y_data, c):
    result = []
    for instant_y in y_data:
        result.append(instant_v(y_data, instant_y, c))
    return np.asarray(result)


def instant_v(y, y_x, c):
    return math.sqrt(2 * GRAVITY * abs(y[0] - y_x) / (1 + c))


def v_plot(c, v_func):
    y_numeric = v_func(y_n, x_n, c)
    y_experimental = v_func(y_e, x_e, c)
    x_label = r'$x$ (m)'
    y_label = r'$v(x)$ m/s'
    plot_data(x_n, y_numeric, x_e, y_experimental, x_label, y_label)


def centripetal_acceleration(velocity, curve):
    result = []
    for i in range(len(velocity)):
        result.append(curve[i] * (velocity[i]**2))
    return np.asarray(result)


def normal_force(mass, acceleration, beta_data):
    result = []
    for i in range(len(acceleration)):
        result.append(mass * (GRAVITY * math.cos(beta_data[i] * 2 * PI / 180) + acceleration[i]))
    return np.asarray(result)


# BEGIN NORMAL FORCE PLOT
def normal_force_plot(mass_object, c):
    velocity_n = v(y_n, c)
    velocity_e = v(y_e, c)
    curve_n = curve(d2y_n, dy_n)
    curve_e = curve(d2y_e, dy_e)
    acceleration_n = centripetal_acceleration(velocity_n, curve_n)
    acceleration_e = centripetal_acceleration(velocity_e, curve_e)
    normal_n = normal_force(mass=mass_object, acceleration=acceleration_n, beta_data=beta(dy_n))
    normal_e = normal_force(mass=mass_object, acceleration=acceleration_e, beta_data=beta(dy_e))

    x_label = r'$x$ (m)'
    y_label = r'$N(X)$, N'
    name = "normal_plot"

    plot_data(x_n, normal_n, x_e, normal_e, x_label, y_label, name)


def friction(c, mass, beta_data):
    result = []
    for value in beta_data:
        result.append(c * mass * GRAVITY * math.sin(value * 2 * PI / 180) / (1 + c))
    return np.asarray(result)


def horizontal_velocity(velocity, beta_data):
    v_avg = [0]
    for i in range(1, len(beta_data)):
        v_avg.append(1/2 * (v_avg[i - 1] + instant_horizontal_velocity(beta_data[i], velocity[i])))
    return np.asarray(v_avg)


def time(velocity):
    time = [0]
    v_avg = horizontal_velocity(velocity, beta(dy_n))
    for i in range(1, len(velocity)):
        time.append(((1/1000) / v_avg[i]) + time[i - 1])
    return np.asarray(time)


def instant_horizontal_velocity(beta_data, velocity):
    return velocity * math.cos(beta_data * 2 * PI / 180)


def friction_plot(mass_object, c):
    y_numeric = friction(c=c, mass=mass_object, beta_data=beta(dy_n))
    y_experimental = friction(c=c, mass=mass_object, beta_data=beta(dy_e))
    x_label = r'$x$ (m)'
    y_label = r'$friction$'
    name = "friction_plot"
    plot_data(x_n, y_numeric, x_e, y_experimental, x_label, y_label, name)


def friction_normal_plot(mass_object, c):
    velocity_n = v(y_n, c)
    velocity_e = v(y_e, c)
    curve_n = curve(d2y_n, dy_n)
    curve_e = curve(d2y_e, dy_e)
    acceleration_n = centripetal_acceleration(velocity_n, curve_n)
    acceleration_e = centripetal_acceleration(velocity_e, curve_e)
    norm_n = normal_force(mass_object, acceleration_n, beta(dy_n))
    norm_e = normal_force(mass_object, acceleration_e, beta(dy_e))
    friction_n = friction(c=c, mass=mass_object, beta_data=beta(dy_n))
    friction_e = friction(c=c, mass=mass_object, beta_data=beta(dy_e))

    y_numeric = np.asarray([abs(friction_n[i] / norm_n[i]) for i in range(len(norm_n))])
    y_experimental = np.asarray([abs(friction_e[i] / norm_e[i]) for i in range(len(norm_e))])

    x_label = r'$x$ (m)'
    y_label = r'$|f/N|$'
    name = "friction_normal_plot"
    y_lim = [0, 4]
    plot_data(x_n, y_numeric, x_e, y_experimental, x_label, y_label, name, y_lim=y_lim)


def velocity_horizontal_displacement_plot(c):
    time_n = time(velocity=v(y_n, c))
    time_e = time_data
    x_numeric = time_n
    x_experimental = time_e
    y_numeric = x_n
    y_experimental = x_e
    x_label = r'$t$ (s)'
    y_label = r'$x$ (m)'
    name = "velcity_horizontal_displacement_plot"
    plot_data(x_numeric, y_numeric, x_experimental, y_experimental, x_label, y_label, name)


def velocity_time_plot(c):
    velocity_n = v(y_n, c)
    velocity_e = v(y_e, c)
    time_n = time(velocity=v(y_n, c))
    time_e = time_data
    x_numeric = time_n
    x_experimental = time_e
    y_numeric = horizontal_velocity(velocity_n, beta(dy_n))
    y_experimental = horizontal_velocity(velocity_e, beta(dy_e))

    x_label = r'$t$ (s)'
    y_label = r'$v(t)$ (m/s)'
    name = "velocity_time_plot"

    plot_data(x_numeric, y_numeric, x_experimental, y_experimental, x_label, y_label, name)


def plot_data(x_numeric, y_numeric, x_experimental, y_experimental, x_label, y_label, name="", y_lim=None, compare_with=0):
    plot = plt.figure(x_label, figsize=(12, 3))
    plt.plot(x_numeric, y_numeric, x_experimental, y_experimental)
    plt.grid()
    plt.xlabel(x_label, fontsize=20)
    plt.ylabel(y_label, fontsize=20)

    if y_lim is not None:
        plt.gca().set_ylim(y_lim)

    plt.show()

    if BALL:
        name = "ball_" + name
    else:
        name = "circle_" + name
    name += "_comparison"
    if name != "":
        plot.savefig(f"{name}.pdf", bbox_inches="tight")


track(x_n, y_n, xfast, yfast)
beta_plot()
curve_plot()
normal_force_plot(mass_object=MASS, c=c)
friction_plot(c=c, mass_object=MASS)
friction_normal_plot(mass_object=MASS, c=c)
velocity_horizontal_displacement_plot(c=c)
velocity_time_plot(c=c)



