import numpy
import math


def f(x):
    return pow(x, 3)-7*x-6


def split_method(eps, a_k=0, b_k=5):
    k = 0
    x_k = (a_k + b_k) / 2
    print("Expected numbers of iteration >="+str(int(math.log((b_k-a_k)/eps, 2))))
    print("         x_k", "  step", "        f(x_k)", sep="                   ")

    while not (-eps <= f(x_k) <= eps):
        if numpy.sign(f(a_k)) == numpy.sign(f(x_k)):
            a_k = x_k
        if numpy.sign(f(b_k)) == numpy.sign(f(x_k)):
            b_k = x_k
        k += 1
        x_k = (a_k + b_k) / 2
        print(str(x_k).center(20)+str(k).center(30)+str(f(x_k)).center(35))
    return x_k


def f_1(x):
    return pow(x, 3)-6*pow(x, 2)+5*x+12


def f_1_diff(x):
    return 3*pow(x, 2)-12*x+5


def f_1_2_diff(x):
    return 6*x-12


def newton(eps, x_k=-1.5):
    k = 0
    a = -2.5
    b = -0.5
    delta = min(x_k - a, b - x_k)
    if x_k > b or x_k < a or (f_1(x_k)*f_1_2_diff(x_k) <= 0):
        print("Wrong first approximation!")
        return 0
    m2 = abs(f_1_2_diff(x_k-delta))
    m1 = abs(f_1_diff(x_k+delta))
    q = (m2 * abs(x_k + delta)) / (2*m1)
    t = numpy.log(((abs(x_k + delta) / eps) / numpy.log(1/q))+1)
    number = math.log(t, 2)
    print("Expected numbers of iteration >=" + str(int(number+1)))
    print("         x_k", "  step", "        f(x_k)", sep="                   ")
    while not (-eps <= f_1(x_k) <= eps):
        x_k -= (f_1(x_k)/f_1_diff(x_k))
        k += 1
        print(str(x_k).center(20) + str(k).center(30) + str(f_1(x_k)).center(35))
    return x_k


def f_2(x):
    return pow(x, 3)+3*pow(x, 2)-x-3


def phi(x):
    return (2*x*x*x+3*x*x+3)/(3*x*x+6*x-1)


def phi_diff(x):
    return 6*(pow(x, 4)+4*pow(x, 3)+2*x*x-4*x-3)/((3*x*x+6*x-1)*(3*x*x+6*x-1))


def simple_iteration(eps, x_k=1.2, a=0.8, b=1.6):
    k = 0

    if x_k > b or x_k < a:
        print("Wrong first approximation!")
        return 0
    delta = min(x_k - a, b - x_k)
    q = max(abs(phi_diff(a)), abs(phi_diff(b)))
    print("q = "+str(q))
    if q > 1 or q < 0:
        print("Theorem doesnt work because q is is negative or more than 1!")
        return 0
    if abs(phi(x_k)-x_k) > (1-q)*delta:
        print("Theorem doesnt work!")
        return 0
    print("Expected numbers of iteration >=" + str(int((numpy.log((b - a) / eps*(1-q)))/numpy.log(1/q))+1)+" >= "+str(int((numpy.log((abs(phi(x_k) - x_k)) / eps*(1-q)))/numpy.log(1/q))+1))
    print("         x_k", "  step", "        f(x_k)", sep="                   ")
    while not (-eps <= x_k-phi(x_k) <= eps):
        x_k = phi(x_k)
        k += 1
        print(str(x_k).center(20) + str(k).center(30) + str(f_2(x_k)).center(35))
    return x_k


print("Split method")
split_method(0.001)
print("Newton")
newton(0.001)
print("Simple iteration method")
simple_iteration(0.001)

