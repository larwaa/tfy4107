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

# Horisontal avstand mellom festepunktene er 200 mm
h = 200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])
# Vi begrenser starthøyden (og samtidig den maksimale høyden) til
# å ligge mellom 250 og 300 mm
ymax = 300
# yfast: tabell med 8 heltall mellom 50 og 300 (mm); representerer
# høyden i de 8 festepunktene
yfast=np.asarray(np.random.randint(50, ymax, size=8))
# inttan: tabell med 7 verdier for (yfast[n+1]-yfast[n])/h (n=0..7); dvs
# banens stigningstall beregnet med utgangspunkt i de 8 festepunktene.
inttan = np.diff(yfast)/h
attempts=1
# while-løkken sjekker om en eller flere av de 3 betingelsene ovenfor
# ikke er tilfredsstilt; i så fall velges nye festepunkter inntil
# de 3 betingelsene er oppfylt


# Når programmet her har avsluttet while-løkka, betyr det at
# tallverdiene i tabellen yfast vil resultere i en tilfredsstillende bane. 

# Omregning fra mm til m:
xfast = xfast/1000
yfast = yfast/1000

yfast = np.asarray([0.297, 0.249, 0.235, 0.161, 0.175, 0.122, 0.193, 0.159])
#Programmet beregner deretter de 7 tredjegradspolynomene, et
#for hvert intervall mellom to nabofestepunkter.

#Med scipy.interpolate-funksjonen CubicSpline:
cs = CubicSpline(xfast, yfast, bc_type='natural')
xmin = 0.000
xmax = 1.401
dx = 0.001
x = np.arange(xmin, xmax, dx)
Nx = len(x)
y = cs(x)
dy = cs(x,1)
d2y = cs(x,2)

#Plotting

baneform = plt.figure('y(x)',figsize=(12,3))
plt.plot(x,y,xfast,yfast,'*')
plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$y(x)$ (m)',fontsize=20)
plt.ylim(0,0.350)
plt.grid()
plt.show()
#baneform.savefig("baneform.pdf", bbox_inches='tight')
#baneform.savefig("baneform.png", bbox_inches='tight')


def beta():
	return np.asarray([math.atan(derivative)/(2 * PI) * 180 for derivative in dy])

# def v_x(x, c):
	# return math.sqrt(2 * GRAVITY * (yfast[0] - y(x)) / (1 + c))


# BETA FØRSTE PLOT
plt.plot(x, beta())
plt.grid()
plt.xlabel('$x$ (m)', fontsize=20)
plt.ylabel(r'$\beta$ (grader)', fontsize=20)
plt.show()
# END BETA FØRSTE PLOT


# START KRUMNINGSPLOT
def curve(d2y, dy):
	return d2y/(1 + dy ** 2) ** (3/2)


plt.plot(x, curve(d2y, dy))
plt.grid()
plt.xlabel(r'$x$(m)', fontsize=20)
plt.ylabel(r'$\kappa$(x) (1/m)')
plt.show()
# END KRUMNINGSPLOT


def v(y, x, c):
	a = []
	for instant_y in y:
		a.append(instant_v(y, instant_y, c))
	return np.asarray(a)


def instant_v(y, y_x, c):
	return math.sqrt(2 * GRAVITY * (y[0] - y_x) / (1 + c))


print(GRAVITY)
c = 2/5
plt.plot(x, v(y, x, c))
plt.grid()
plt.xlabel(r'$x$ (m)', fontsize=20)
plt.ylabel(r'$v(x)$ m/s', fontsize=20)
plt.show()

def t(x, y, c):
	t = [0]
	for n in range(1,1401):
		t.append(instant_t(x, y, c, n))
	return np.asarray(t)

def instant_t(x, y, c, n):
	v_xn_prev = v(y, x, c)[n-1] * np.cos(beta()[n-1]*PI/180)
	v_xn = v(y, x, c)[n] * np.cos(beta()[n]*PI/180)
	v_xn_avg = (1/2) * (v_xn_prev + v_xn)
	return dx/v_xn_avg

plt.plot(t(x, y, c), x)
plt.grid()
plt.xlabel(r'$t$ (s)', fontsize=20)
plt.ylabel(r'$x$ m', fontsize=20)
plt.show()